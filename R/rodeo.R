###---------------------------------------------###
### rodeo (regularised ODE optimisation) object ###
###---------------------------------------------###

# Sanity check generic
#
# @description Performs sanity checks on the objects of the package. Highly irrelevant for most users.
#
# @param x The object to check.
# @param ... Additional arguments passed to \code{check}.
# @return check returns nothing, only throws errors and warnings.
# @keywords internal
.check <- function (x, ...) {
  UseMethod(".check", x)
}

# Sanity check of rodeo object
#
# @param x Object of class \code{\link{rodeo}}.
# @param ... Additional arguments passed to \code{.check.rodeo}.
# @inherit .check description return
# @keywords internal
.check.rodeo <- function(x, ...) {
  with(x, {
    .o_op_check(o, op)
    s <- op$s
    d <- o$d
    if (is.array(x0s) || is.matrix(x0s)) {
      if (dim(x0s)[1] != d * s) stop(paste0("Number of rows of x0s (", nrow(x0s), ") does not equal d * s (", d*s, "), extracted from 'ode'-object and y."))
    } else {
      stop("x0s is neither matrix nor array.")
    }
  })
}

.exclude <- function (x, ...) {
  UseMethod(".exclude", x)
}


#' Regularised Ordinary Differential Equation Optimisation (RODEO) generic
#'
#' @description Fit the parameters (and optionally initial states) for an Ordinary Differential Equation model with data sampled across contexts. Two implementations exist:
#' \describe{
#' \item{\code{\link{rodeo.ode}}}{The raw form of rodeo, based on \code{\link{ode}} (created via \code{\link{mak}}, \code{\link{plk}}, etc.) and \code{\link{opt}}-objects. If needed, everything can be user-specified through these two objects. However, parameter initialisations are always required and no action towards adjusting weights or scales of the parameters are taken.}
#' \item{\code{\link{rodeo.aim}}}{The more automatic form of rodeo, where automatic parameter initialisations and (optional) adjustments of scales and weights are available through \code{\link{aim}}.}
#' }
#' @param x Object that specifies ode system.
#' @param ... Additional arguments passed to \code{rodeo}.
#' @seealso rodeo.aim, rodeo.ode
#' @return An object with S3 class "rodeo":
#' \item{o}{Original \code{ode}-object.}
#' \item{op}{Original \code{opt}-object with default values for \code{lambda_min_ratio}, \code{lambda} and (if needed) \code{a} inserted, if these were originally \code{NULL}.}
#' \item{params}{Parameter estimates, stored as list of sparse column format matrices, "dgCMatrix" (or a list of those if multiple initialisations). Rows represent coordinates and columns represent the \code{lambda} value.}
#' \item{x0s}{Initial state estimates stored in a matrix (or array). Rows represent coordinates, columns represent the \code{lambda} value and (if multiple initialisations) slices represent initialisations.}
#' \item{dfs}{A matrix (or array, if multiple initialisations) of degrees of freedom. Row represents a parameter (the first is always the initial state parameter), columns represent lambda, slices represent initialisation, if multiple are provided.}
#' \item{codes}{A matrix (or array) of convergence codes organised as \code{dfs}.
#' \describe{
#' \item{0:}{The convergence criteria is met (see details in \code{\link{opt}}). Current estimate is probably a local minimum.}
#' \item{1:}{Backtracking in the last iteration yields no numerical improvement, but no unsual behavior observed. Current estimate is probably a local minimum. However, if \code{exact_gradient} = FALSE in the \code{reg}-object in the \code{ode}-object, changing this may improve the code. Alternatively one can adjust backtracking via \code{backtrack_max} and \code{tau_min} in \code{reg} objects in \code{ode} object.}
#' \item{2:}{The optimisation procedure exceeded maximal number of steps (\code{step_max} in \code{reg} objects).}
#' \item{3:}{The last gradient was unsually large. Either the tolerances in \code{reg} objects are off or the ODE systems is very sensitive and runs over long time spans. In the latter case, initialisation(s) may have inappropriate zeros (change initialisation and/or make sure they start at smaller lambda value).}
#' \item{4:}{The numeric ODE solver exceeded maximal number of steps. Check if supplied initial states were out of bounds, if not increase \code{step_max} (or \code{tol}) in \code{reg}-objects in \code{ode}-object.}
#' }
#'}
#' \item{steps}{A matrix (or array) holding number of steps used in optimisation procedure. Organised as \code{dfs}.}
#' \item{losses}{A vector (or matrix) of unpenalised losses at optimum for each lambda value (stored row-wise if  multiple are provided).}
#' \item{penalties}{A matrix (or array) of penalties for each parameter, organised as \code{dfs}.}
#' \item{jerr}{A matrix (or array) of summary codes (for internal debugging), organised as \code{dfs}.}
#'
#' @details
#' For details on the loss function, optimisation, etc. See documentation of \code{\link{opt}}.
#'
#' @examples
#' set.seed(123)
#' # Michaelis-Menten system with two 0-rate reactions
#' A <- matrix(c(1, 1, 0, 0,
#'  0, 0, 1, 0,
#'  0, 0, 1, 0,
#'  0, 0, 0, 1,
#'  1, 0, 0, 0), ncol = 4, byrow = TRUE)
#' B <- matrix(c(0, 0, 1, 0,
#'  1, 1, 0, 0,
#'  1, 0, 0, 1,
#'  0, 0, 1, 0,
#'  1, 0, 1, 0), ncol = 4, byrow = TRUE)
#' k <- c(1.1, 2.5, 0.25, 0, 0); x0 <- c(E = 2, S = 2, ES = 8.5, P = 2.5)
#' Time <- seq(0, 10, by = 1)
#'
#' # Simulate data, in second context the catalytic rate has been doubled
#' contexts <- cbind(1, c(1, 1, 2, 1, 1))
#' m <- mak(A, B, r = reg(contexts = contexts))
#' y <- numsolve(m, c(Time, Time), cbind(x0, x0 + c(2, -2, 0, 0)), contexts * k)
#' y[, -1] <- y[, -1] + matrix(rnorm(prod(dim(y[, -1])), sd = .5), nrow = nrow(y))
#'
#' # Example: fit data using rodeo on mak-object
#' op <- opt(y, nlambda = 10)
#' fit <- rodeo(m, op, x0 = NULL, params = NULL)
#' fit$params$rate
#'
#' # Example: fit dat using rodeo on aim-object
#' a <- aim(m, op)
#' a$params$rate
#' fit <- rodeo(a)
#' fit$params$rate
#'
#' @export
rodeo <- function (x, ...) {
  UseMethod("rodeo", x)
}



#' Regularised Ordinary Differential Equation Optimisation (RODEO)
#'
#' @description Fit the parameters (and optionally initial states) for a ordinary differential equation model with data sampled across different contexts.
#' @param x \code{ode}-object created via \code{\link{mak}}, \code{\link{plk}} etc.
#' @param op \code{opt}-object created via \code{\link{opt}}-function (compatibility check with \code{x} is conducted).
#' @param x0 A vector (or matrix) of non-negative initialisation(s) for the initial state in the optimisation problem, see details.
#' @param params A list of initialisations for the parameter arguments in the optimisation problem, see details. A list of matrices if multiple initialisations are desired. Sparse matrix "dgCMatrix" is allowed.
#' @param trace Logical indicating if status messages should be printed during \code{rodeo}.
#' @param ... Additional arguments passed to \code{rodeo}.
#' @seealso rodeo, rodeo.aim
#' @inherit rodeo return
#' @details
#' For running a single initialisation of the optimisation procedure, supply \code{x0} as vector (where the initial states for the contexts are concatinated) and \code{params} as a list with a vector entry for each parameter. The resulting estimates are stored in matrices, with each column representing a lambda-value. The convergence codes, steps and losses are stored as vectors, one entry for each value of lambda.
#'
#' For running multiple initialisations, supply \code{x0} as matrix with the individual initialisations as columns (each column vector as described above) and \code{params} as list of matrices with initialisations stored column-wise.
#'
#' The initial state estimates (for the different lambda-values) are returned as a matrix (array) (row = coordinate, column = lambda-value, (slice = initialisation)) and the parameter estimates are returned as a list of (lists with) sparse matrices (row = coordinate, column = lambda-value). The convergence codes, steps and losses are stored as matrices (row = lambda, column = initialisation).
#'
#' For details on the loss function, optimisation, etc. See documentation of \code{\link{opt}}.
#'
#' If explicitly setting \code{x0} to \code{NULL}, the system uses the first observation of each context (throws error if non finites).
#'
#' If explicitly setting \code{params} to \code{NULL}, 0 initialisations are used.
#'
#' @examples
#' set.seed(123)
#' # Example: Michaelis-Menten system with two 0-rate reactions
#' A <- matrix(c(1, 1, 0, 0,
#'               0, 0, 1, 0,
#'               0, 0, 1, 0,
#'               0, 1, 0, 0,
#'               0, 0, 0, 1), ncol = 4, byrow = TRUE)
#' B <- matrix(c(0, 0, 1, 0,
#'               1, 1, 0, 0,
#'               1, 0, 0, 1,
#'               0, 0, 0, 1,
#'               0, 0, 1, 0), ncol = 4, byrow = TRUE)
#' k <- c(1, 2, 0.5, 0, 0); x0 <- c(E = 2, S = 8, ES = 0.5, P = 0.5)
#' Time <- seq(0, 10, by = 1)
#'
#' # Simulate data, in second context the catalytic rate has been doubled
#' m <- mak(A, B)
#' contexts <- cbind(1, c(1, 1, 2, 1, 1))
#' y <- numsolve(m, c(Time, Time), cbind(x0, x0 + c(2, -2, 0, 0)), contexts * k)
#' y[, -1] <- y[, -1] + matrix(rnorm(prod(dim(y[, -1])), sd = .1), nrow = nrow(y))
#'
#' # Fit data using rodeo on mak-object
#' op <- opt(y)
#' fit <- rodeo(m, op, x0 = NULL, params = NULL)
#'
#' # Example: fit data knowing doubled catalytic rate
#' m_w_doubled <- mak(A, B, r = reg(contexts = contexts))
#' fit <- rodeo(m_w_doubled, op, x0 = NULL, params = NULL)
#'
#' @export
#' @useDynLib episode, .registration = TRUE
#' @importFrom Rcpp sourceCpp
rodeo.ode <- function(x, op, x0, params, trace = FALSE, ...) {
  .o_op_check(x, op)

  s <- op$s
  d <- x$d
  ps <- lapply(x$rs, getElement, name = "p")

  # params NULL handling
  if (is.null(params)) {
    params <- sapply(ps[-1], rep, x = 0, simplify = FALSE)
  }

  if (!methods::is(params, "list")) {
    stop("params is not list.")
  } else if (length(ps) - 1 != length(params)) {
    stop(paste0("params does not have length equal to that dictated by ode-object (", length(ps) - 1, ")."))
  }

  N <- 1
  are_matrices <- unlist(lapply(params, function(par) {is.matrix(par) || methods::is(par, "dgCMatrix")}))
  if (all(are_matrices)) {
    N <- unlist(lapply(params, ncol))
    if (range(N)[1] != range(N)[2]) {
      stop("Not all parameters have the same number of columns.")
    }
    N <- N[1]
    if (!all(unlist(Map("==", x = lapply(params, nrow), y = ps[-1])))) {
      stop(paste0("Number of rows in params (", paste(unlist(lapply(params, nrow)), collapse = ", "), ") are not equal to length of parameters according to ode-object (", paste(unlist(ps[-1]), collapse = ", "), ")"))
    }
    multi <- TRUE
  } else if (all(!are_matrices)) {
    if (!all(unlist(Map("==", x = lapply(params, length), y = ps[-1])))) {
      stop(paste0("Number of rows in params (", paste(unlist(lapply(params, length)), collapse = ", "), ") are not equal to length of parameters according to ode-object (", paste(unlist(ps[-1]), collapse = ", "), ")"))
    }
    multi <- FALSE
  } else {
    stop("params must only contain matrices or vectors and not both.")
  }

  # x0 NULL handling
  if (is.null(x0)) {
    x0 <- as.vector(t(op$y[c(0, which(diff(op$y[, 1]) < 0)) + 1, -1]))
    if (any(!is.finite(x0))) stop("x0 is NULL, but no default replacement exists, since some contexts have first observations which are non-finite.")
    x0 <- pmax(x0, x$rs[[1]]$ctrl$tol)
    if (multi) {
      x0 <- matrix(x0, nrow = length(x0), ncol = N)
    }
  }

  # Check params and x0 compatibility (either both matrix or vector)
  if (any(unlist(lapply(params, function(para) any(!is.finite(para)))))) stop("params contains non-finite values.")

  ## Check x0
  if (multi) {
    if (ncol(x0) != N) {
      stop(paste0("Number of columns in x0 (", ncol(x0), ") does not equal number of columns in parameters (", N, ")"))
    }
    if (nrow(x0) != (d * s)) stop(paste0("The number of rows in x0 (", nrow(x0), ") is not equal to d*s (", d * s, "), where d is dimension of state space, s = number of contexts."))
    x0 <- array(x0, dim = c(d, s, ncol(x0)))
  } else {
    if (length(x0) != (d * s)) stop(paste0("The length of x0 (", length(x0), ") is not equal to d*s (", d * s, "), where d is dimension of state space, s = number of contexts."))
    x0  <- array(x0, dim = c(d, s, 1))
  }

  op <- .set_defaults(x, op)#, x0, set_lambda = TRUE) # not before here, since x0 is first available here

  # make params into list of lists
  names_param <- names(params)
  if (multi) {
    params <- lapply(params, function(para) as.list(data.frame(para)))
    params <- .switch_hierarchy(params) # switch hierarchy, so it is: particle, parameter
  } else {
    params <- list(params)
  }

  # The user operates with params on scales provided in rs, but cpp code does not
  params <- lapply(params, function(para) {
    Map(function(p, s) {
      if (!is.null(s)) {
        p <- p / s
      }
      return(p)
      }, p = para, s = lapply(x$rs[-1], getElement, "scales"))
  })

  sc_ <- lapply(x$rs[-1], .get_sc, op = op)

  ## Actual rodeo part
  if (any(unlist(lapply(x$rs, function(re) length(re$exclude) > 0)))) {
    sc_ <- Map(function(sc, exc) {
      if (!is.null(sc) & !is.null(exc) & length(exc) > 0) sc <- sc[-exc, , drop = FALSE]
      return(sc)
      }, sc = sc_, exc = lapply(x$rs[-1], getElement, "exclude"))
    params <- lapply(params, function(particle) {
      Map(function(para, exc) {
        if (!is.null(exc) & length(exc) > 0) para <- para[-exc]
        return(para)
      }, para = particle, exc = lapply(x$rs[-1], getElement, "exclude"))
    })

    ode_struct <- x
    ode_struct <- .exclude(ode_struct)

    rodeo_ret <- bull(ode_struct = ode_struct, opt_struct = op, sc_ = sc_, params_ = params, x0s = x0, trace = trace)
  } else {
    rodeo_ret <- bull(ode_struct = x, opt_struct = op, sc_ = sc_, params_ = params, x0s = x0, trace = trace)
  }

  # make x0 estimates into array
  x0s <- array(do.call(c, lapply(rodeo_ret$params, function(x) as.vector(x[[1]]))),
    dim = c(nrow(rodeo_ret$params[[1]][[1]]), ncol(rodeo_ret$params[[1]][[1]]), length(rodeo_ret$params)))
  param_est <- lapply(rodeo_ret$params, function(x) x[-1])

  ## Prepare return
  if (multi && length(param_est) > 1) {
    # Means N > 1

    # Switcharoo (so now order is: parameter, then particle)
    param_est <- .switch_hierarchy(param_est)

    # add excluded (if null, nothing happens)
    param_est <- Map(function(para, re) {
      lapply(para, .add_exclude, re$exclude)
    }, para = param_est, re = x$rs[-1])

    # back on scale
    param_est <- Map(function(para, re) {
      if (!is.null(re$scales)) {
        lapply(para, function(pa) pa * re$scales)
      } else {
        para
      }
      }, para = param_est, re = x$rs[-1])

    # names
    names(param_est) <- names_param

    ret <- structure(list(o = x, op = op, params = param_est, x0s = x0s,
      dfs = rodeo_ret$dfs, codes = rodeo_ret$codes, steps = rodeo_ret$steps,
      losses = rodeo_ret$losses, penalties = rodeo_ret$penalties,
      jerr = rodeo_ret$jerr), class = "rodeo")
  } else {
    # N = 1
    param_est <- param_est[[1]]
    x0s <- x0s[,,1]

    # add excluded (if null, nothing happens)
    param_est <- Map(.add_exclude, x = param_est, exclude = lapply(x$rs[-1], getElement, "exclude"))

    # back on scale
    param_est <- Map(function(para, re) {
      if (!is.null(re$scales)) {
        para * re$scales
      } else {
        para
      }
      }, para = param_est, re = x$rs[-1])

    # names
    names(param_est) <- names_param

    ret <- structure(list(o = x, op = op, params = param_est, x0s = x0s,
      dfs = rodeo_ret$dfs[,,1], codes = rodeo_ret$codes[,,1], steps = rodeo_ret$steps[,,1],
      losses = rodeo_ret$losses[1,], penalties = rodeo_ret$penalties[,,1],
      jerr = rodeo_ret$jerr[,,1]), class = "rodeo")
  }

  # adjust lambda factors
  ret$op$lambda  <- as.vector(rodeo_ret$lambda)
  ret$op$nlambda <- length(ret$op$lambda)
  if (!is.null(rodeo_ret$lambda_fac)) {
    ret$o$rs <- Map(function(re, lf) {
      re$lambda_factor <- lf
      re
    }, ret$o$rs, lf = as.list(rodeo_ret$lambda_fac))
  }

  return(ret)
}


