###-----------------------------------------###
### aim (adaptive integral matching) object ###
###-----------------------------------------###

# switch hierarchy in list of lists, i.e.,
# list(l1 = (l1$a, l2$b, ...), l2 = (l2$a, l2$b, ...), ...) becomes
# list(a = (l1$a, l2$a, ...), b = (l2$b, l2$b, ...), ...)
# requires length of all sub lists to be equal
.switch_hierarchy <- function(ll) {
  if (!is.list(ll)) stop("ll is not list.")
  if (!all(unlist(lapply(ll, is.list)))) stop("Not all elements of ll are lists.")
  lens <- lapply(ll, length)
  if (range(lens)[1] != range(lens)[2]) stop("Not all lists in ll have same length.")

  return(do.call(function(...) Map(list, ...), ll))
}


# add excluded rows to sparse matrix x
.add_exclude <- function(x, exclude) {
  if (length(exclude) > 0) {
    p <- nrow(x) + length(exclude)
    if (!all(exclude %in% seq_len(p))) stop("Some 'exclude's are oob.")
    if (!identical(exclude, unique(exclude))) stop("Exclude does not contain unique elements.")
    x <- methods::as(x, "dgCMatrix")
    x@Dimnames <- list(NULL, NULL)
    x@Dim[1] <- x@Dim[1] + length(exclude)
    x@i <- setdiff(seq_len(p), exclude)[x@i + 1] - 1L
  }
  return(x)
}

# Sanity check of aim object
# @param x Object of class \code{\link{aim}}.
# @param ... Additional arguments passed to \code{.check.aim}.
# @inherit .check description return
# @keywords internal
.check.aim <- function(x, ...) {
  with(x, {
    .o_op_check(o, op)
    if (!is.null(op$lambda)) {
      .check_paramx0s(params, x0s, length(op$lambda))
    } else {
      .check_paramx0s(params, x0s, op$nlambda)
    }
  })
}

# Checks dimension on params and x0s
.check_paramx0s <- function(params, x0s, n_lambda) {
  if (n_lambda != ncol(x0s)) stop(paste0("Number of lambda values (nlambda = ", n_lambda , ") does not match number of initial state initialisations (ncol(x0s) = ", ncol(x0s), ")"))

  lapply(params, function(para) {
    if (n_lambda != ncol(para)) stop(paste0("Number of lambda values (", n_lambda , ") does not match number of parameter initialisations (ncol(params) = ", ncol(para), ")"))
    if (!methods::is(para, "dgCMatrix")) stop("'params' in 'aim'-object is not of class 'dgCMatrix'.")
  })
}

# x-smooth format check
.check_x <- function(x, op) {
  if (!is.matrix(x) || !is.numeric(x)) {
    stop("The smoothed curve, x, is not a matrix nor numeric.")
  }
  if (any(!is.finite(x[, 1]))) {
    stop("First column of smoothed curve, x, (time points) have non finites or missing values!")
  }
  if (nrow(x) <= 1) {
    stop("The smoothed curve, x, must contain at least two rows.")
  }
  if (ncol(x) != ncol(op$y)) {
    stop(paste0("Number of columns of the smoothed curve, x, (", ncol(x), ") does not match number of columns of y (", ncol(op$y), ")."))
  }
  s_smooth  <- sum(diff(x[, 1]) < 0) + 1
  if (s_smooth != op$s) {
    stop(paste0("Number of series in the smoothed curve, x, (", s_smooth, ") does not match number of series in y (", op$s, ")."))
  }
}

# x-smooth NA replacement
.na_x <- function(x, was_NULL) {
  times <- x[, 1]
  x_obs <- x[, -1, drop = FALSE]
  brks  <- c(1, which(diff(times) < 0) + 1, length(times) + 1)  # each series is brks[j] to brks[j+1]-1
  s <- length(brks) - 1
  n <- nrow(x_obs)

  na_coor <- apply(!is.na(x_obs), 2, sum)
  for (j in seq_len(s)) {
    ind <- seq(brks[j], brks[j + 1] - 1)

    x_obs[ind, ] <- apply(x_obs[ind, , drop = FALSE], 2, function(x) {
      isnt_na <- !is.na(x)
      if (sum(isnt_na) == 0) {
        stop(ifelse(was_NULL,
          "AIM not compatible with contexts of y with fully unobserved coordinates, when x is not supplied.",
          "AIM not compatible with contexts of x with fully unobserved coordinates."))
      }

      if (sum(isnt_na) == 1) {
        x <- rep(x[which(isnt_na)], length(ind))
      } else if (sum(isnt_na) < length(ind)) {
        x <- stats::approx(x = times[ind], y = x, xout = times[ind], rule = 2)$y
      }
      return(x)
    })
  }
  if (any(na_coor < n)) {
    # Warning if some where interpolated
    str_ <- "Some coordinate(s) of"
    str_ <- paste(str_, ifelse(was_NULL, "y (which will be used to create x)", "x"))
    str_ <- paste(str_, "had missing values. Linear interpolation will be used to replace these.")
    warning(str_)
  }

  return(cbind(times, x_obs))
}

# x-smooth NULL handling, checking and NA handling
.prepare_x <- function(x, op) {
  if (is.null(x)) {
    x  <- op$y
    return(.na_x(x, was_NULL = TRUE))
  } else {
    .check_x(x, op)
    return(.na_x(x, was_NULL = FALSE))
  }
}



#' Integral Matching Design
#'
#' @description Get the design matrix, formated data and weights used in integral matching.
#' @param o An object of class \code{\link{ode}} (created via \code{\link{mak}}, \code{\link{plk}}, etc.).
#' @param op An object of class \code{\link{opt}} (created via \code{\link{opt}}). Compatibility check with \code{o} is conducted.
#' @param x Matrix of dimension m-x-(d+1) containing custom time and smoothed values of the processes used for the integral matching, see details. If \code{NULL} (default) a linear-interpolation smoother is employed.
#' @param params List of vectors of initial parameter values. Default (\code{NULL}) translates to the zero-vector. If the ODE system is linear in the parameter vector and the boundary constraints on are not unusual, then \code{params} is ignored.
#' @seealso aim, rodeo.aim, numint
#' @details
#' The design matrix is as follows:
#' \deqn{X = (\int_{t_{i-1}}^{t_{i}}{\frac{df}{dparam}(x(s), context_l * param) * context_l ds})_{i=1}^n}
#' Here f is the vector field of the ODE-system, x is the smoothed process and params is (internally) scaled with \code{scales} in \code{\link{reg}} objects in \code{o}.
#'
#' Similiarly the observations and weights are concatinateed as follows:
#' \deqn{Y = (x_{t_i} - x_{t_{i-1}})_{i=1}^n}
#' \deqn{W = ((w_{t_i} + w_{t_{i-1}}) / 2)_{i=1}^n}
#'
#' The number of decreases in time in \code{x} must match the number of decreases in time in \code{y} in \code{op}. The process x is simply a linear interpolation of \code{x}, hence finer discretisations give more refined integral matching. Each context is handled seperately. Moreover, the compatibility checks in \code{rodeo} are also conducted here.
#'
#' @return A list with:
#' \item{X}{List of design matrices, one for each parameter argument.}
#' \item{Y}{A vector of observations.}
#' \item{W}{A vector of weights.}
#' \item{X0}{A matrix of initial states, one for each context. Interpolated from \code{x}.}
#' \item{xout}{A matrix containing the values of the process at the time points in \code{y}, interpolated from \code{x}. Use it to check that \code{x} works as intended.}
#'
#'
#' @examples
#' set.seed(123)
#' # Michaelis-Menten system with two 0-rate reactions
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
#' k <- c(2.1, 2.25, 1.5, 0, 0); x0 <- c(E = 8, S = 10, ES = 1.5, P = 1.5)
#' Time <- seq(0, 10, by = 1)
#'
#' # Simulate data, in second context the catalytic rate has been doubled
#' contexts <- cbind(1, c(1, 1, 2, 1, 1))
#' m <- mak(A, B, r = reg(contexts = contexts))
#' y <- numsolve(m, c(Time, Time), cbind(x0, x0 + c(2, -2, 0, 0)), contexts * k)
#' y[, -1] <- y[, -1] + matrix(rnorm(prod(dim(y[, -1])), sd = .5), nrow = nrow(y))
#'
#' # Get the design matrix used in integral matching
#' d <- imd(m, opt(y))
#' head(d$X[[1]])
#'
#' # Compare with glmnet
#' lasso <- glmnet::glmnet(x = d$X[[1]], y = d$Y, intercept = FALSE, lower.limits = 0)
#' a <- aim(m, opt(y, nlambda = 100), adapts = "scales")
#' all.equal(lasso$beta, a$params$rate)
#'
#' @export
#' @useDynLib episode, .registration = TRUE
#' @importFrom Rcpp sourceCpp
imd <- function(o, op, x = NULL, params = NULL) {
  ### Do compatibility check of dimension between ode and ro
  .o_op_check(o = o, op = op)

  ### x_smooth NULL handling or checking ###
  x <- .prepare_x(x, op)

  ### Set defaults, signal to .set_defaults that we need to use im-loss by x0 as matrix
  op <- .set_defaults(o = o, op = op)#, x0 = x, set_lambda = FALSE)

  # bring params on raw scale
  if (!is.null(params)) {
    params <- Map(function(para, re) {
      if (!is.null(re$scales)) {
        para <- para / re$scales
      }
      para
    }, para = params, re = o$rs[-1])
  }

  # actual design part
  sc  <- lapply(o$rs[-1], .get_sc, op = op)
  des <- aim_design(ode_struct = o, opt_struct = op, sc_ = sc, x = x, param_ = params)
  if (all(des$W == 1)) des$W <- NULL
  des$Y <- as.vector(des$Y)
  des$W <- as.vector(des$W)
  des$X0 <- matrix(des$X0, nrow = o$d)

  ## cpp operates on unscaled scale (param-cpp = param-here / scales),
  # bring it to scaled setup
  des$X <- Map(function(X, re) {
    if (!is.null(re$scales)) {
      X <- t(t(X) / re$scales)
    }
    X
  }, X = des$X, re = o$rs[-1])

  return(des)
}


## Helper function that does actual aim fitting ###
# aimfit must contain beta (sparse matrix, lambda in column), lambda, jerr,
# dfs (col = lambda, row = parameter, other than x0)
# lambda_factors (vector of new lambda_factors)
# fits (list of n-x-nlambda matrices)
# here params are on the scale as elsewhere in R part, it is scaled internally below
.get_aimfit <- function(o, op, x, alpha, standardize, params = NULL, ...) {
  ps  <- lapply(o$rs[-1], getElement, "p")
  sc  <- lapply(o$rs[-1], .get_sc, op = op)
  n   <- nrow(op$y)

  # bring params on raw scale
  if (!is.null(params)) {
    params <- Map(function(para, re) {
      if (!is.null(re$scales)) {
        para <- para / re$scales
      }
      para
    }, para = params, re = o$rs[-1])
  }

  # Get design
  des <- aim_design(ode_struct = o, opt_struct = op, sc_ = sc, x = x, param_ = params)
  if (all(des$W == 1)) des$W <- NULL
  des$Y <- as.vector(des$Y)
  des$W <- as.vector(des$W)

  # Determine if other packages handles the problem
  use_optim <- TRUE
  if (.get_ode_name(o) %in% c("mak", "plk")) {
    # param-linear setup
    if (o$rs[[2]]$reg_type %in% c("l1", "l2", "elnet")) {
      if (all(o$rs[[2]]$lower <= 0) & all(o$rs[[2]]$upper >= 0)) {
        if (!is.null(op$lambda)){
          op$lambda <- op$lambda * (n - op$s) / n # brings from aim-level to glmnet-level
          if (!is.null(o$rs[[2]]$lambda_factor)) op$lambda <- op$lambda * o$rs[[2]]$lambda_factor
        }

        # use glmnet
        if (ps[[1]] == 1) {
          v   <- ifelse(is.null(o$rs[[2]]$penalty_factor), 1, max(0.001, o$rs[[2]]$penalty_factor))
          yx  <- as.vector(t(des$Y) %*% des$X[[1]])
          xx  <- sum(des$X[[1]]^2)
          lambda_max  <- max(min(abs(yx) / (v * alpha[[2]]), 0.9 * xx/ (v * (1 - alpha[[2]]))), 1e-7) # Should be bulletproof (1st correspond to l1, 2nd to l2)
          lambda  <- exp(log(lambda_max) + seq(0, log(op$lambda_min_ratio), length.out = 5))
          aimfit  <- list(beta = Matrix::Matrix((yx - lambda * v * alpha[[2]] * sign(yx)) *
              (abs(yx) >= lambda * v * alpha[[2]]) / (xx - lambda * v * (1 - alpha[[2]])), nrow = 1), lambda = lambda)
        } else {
          if (is.null(o$rs[[2]]$penalty_factor)) {
            if (is.null(des$W)) {
              aimfit  <- glmnet::glmnet(x = des$X[[1]], y = des$Y, intercept = FALSE, lambda = op$lambda,
                nlambda = op$nlambda, lambda.min.ratio = op$lambda_min_ratio,
                lower.limits = o$rs[[2]]$lower, upper.limits = o$rs[[2]]$upper, alpha = alpha[[2]],
                exclude = o$rs[[2]]$exclude, standardize = standardize, ...)
            } else {
              aimfit  <- glmnet::glmnet(x = des$X[[1]], y = des$Y, intercept = FALSE, lambda = op$lambda,
                nlambda = op$nlambda, lambda.min.ratio = op$lambda_min_ratio,
                lower.limits = o$rs[[2]]$lower, upper.limits = o$rs[[2]]$upper, alpha = alpha[[2]],
                exclude = o$rs[[2]]$exclude, weights = des$W, standardize = standardize, ...)
            }
          } else {
            if (is.null(des$W)) {
              aimfit  <- glmnet::glmnet(x = des$X[[1]], y = des$Y, intercept = FALSE, lambda = op$lambda,
                nlambda = op$nlambda, lambda.min.ratio = op$lambda_min_ratio,
                lower.limits = o$rs[[2]]$lower, upper.limits = o$rs[[2]]$upper, alpha = alpha[[2]],
                exclude = o$rs[[2]]$exclude, penalty.factor = o$rs[[2]]$penalty_factor, standardize = standardize, ...)
            } else {
              aimfit  <- glmnet::glmnet(x = des$X[[1]], y = des$Y, intercept = FALSE, lambda = op$lambda,
                nlambda = op$nlambda, lambda.min.ratio = op$lambda_min_ratio,
                lower.limits = o$rs[[2]]$lower, upper.limits = o$rs[[2]]$upper, alpha = alpha[[2]],
                exclude = o$rs[[2]]$exclude, weights = des$W, penalty.factor = o$rs[[2]]$penalty_factor, standardize = standardize, ...)
            }
          }
        }

        aimfit$lambda <- aimfit$lambda * n / (n - op$s) # brings from glmnet-level to maker-level
        if (!is.null(o$rs[[2]]$lambda_factor)) op$lambda <- op$lambda / o$rs[[2]]$lambda_factor
        aimfit$dfs    <- matrix(aimfit$df, nrow = 1)
        use_optim <- FALSE
      }
    } else if (o$rs[[2]]$reg_type == "none") {
      if (all(o$rs[[2]]$lower == 0) & all(o$rs[[2]]$upper == Inf)) {
        # use nnls
        if (is.null(des$W)) {
          nnlsfit <- nnls::nnls(A = des$X[[1]], b = des$Y)
        } else {
          nnlsfit <- nnls::nnls(A = des$X[[1]] * sqrt(des$W), b = des$Y * sqrt(des$W))
        }
        aimfit <- list(beta = Matrix::Matrix(nnlsfit$x, ncol = 1), lambda = 1, jerr = nnlsfit$mode,
          dfs = matrix(sum(nnlsfit$x != 0), nrow = 1))
        use_optim <- FALSE
      } else if (all(o$rs[[2]]$lower == -Inf) & all(o$rs[[2]]$upper == Inf)) {
        # use lsfit
        lsfitfit <- stats::lsfit(x = des$X[[1]], y = des$Y, wt = des$W, intercept = FALSE)
        aimfit <- list(beta = Matrix::Matrix(lsfitfit$coef, ncol = 1), lambda = 1, jerr = NULL,
          dfs = matrix(sum(lsfitfit$coef != 0), nrow = 1))
        use_optim <- FALSE
      }
    }
  }

  if (use_optim) {
    # either it is param-linear and none of the pre-packaged solutions apply, or it is non param-linear
    # either way, we can send it to cattle (in the linear-param case, 'param' is ignored)

    if (!is.null(params)) {
      if (!all(unlist(Map("==", x = lapply(params, length), y = ps)))) {
        stop(paste0("Length of 'params' (", paste(unlist(lapply(params, length)), collapse = ", "), "), are not equal to ps from reg objects (", paste(unlist(ps), collapse = ", "),")."))
      }

      # Exclude handle
      params <- Map(function(para, exc) {
          if (!is.null(exc) & length(exc) > 0) para <- para[-exc]
          return(para)
        }, para = params, exc = lapply(o$rs[-1], getElement, "exclude"))
    }

    sc_ <- Map(function(sc, exc) {
      if (!is.null(sc) & !is.null(exc) & length(exc) > 0) sc <- sc[-exc, , drop = FALSE]
      return(sc)
    }, sc = sc, exc = lapply(o$rs[-1], getElement, "exclude"))

    ode_struct <- o
    ode_struct <- .exclude(ode_struct)

    cattlefit <- cattle(ode_struct = ode_struct, opt_struct = op, sc_ = sc_,
      param_ = params, x = x, trace = isTRUE(list(...)$trace))
    aimfit <- list(
      beta = Map(.add_exclude, x = cattlefit$params[-1], exclude = lapply(o$rs[-1], function(re) re$exclude)),
      lambda = as.vector(cattlefit$lambda), jerr = cattlefit$jerr,
      dfs = cattlefit$dfs, lambda_factors = cattlefit$lambda_fac)

    # Save the last des
    des <- aim_design(ode_struct = o, opt_struct = op, sc_ = sc, x = x,
      param_ = lapply(aimfit$beta, function(para) para[, ncol(para)]))
    if (all(des$W == 1)) des$W <- NULL
    des$Y <- as.vector(des$Y)
    des$W <- as.vector(des$W)

  } else {
    # if we used the ones above, we only had one parameter to return
    aimfit$beta <- list(aimfit$beta)

    # also no change in lambda_factors, so
    aimfit$lambda_factors <- NULL
  }

  # fitted values (matrices)
  aimfit$fits <- Map(function(X, beta) {
    X %*% beta
  }, X = des$X, beta = aimfit$beta)

  # aimfit back to scale
  aimfit$beta <- Map(function(para, re) {
    if (!is.null(re$scales)) {
      para <- para * re$scales
    }
    para
  }, para = aimfit$beta, re = o$rs[-1])


  return(list(aimfit = aimfit, des = des))
}



#' Adaptive Integral Matching (AIM)
#'
#' @description Gives approximate parameter estimates using integral matching and optionally adapts weights and scales to these. Feed this to \code{\link{rodeo}} for initialising exact estimation.
#' @param o An object of class \code{\link{ode}} (created via \code{\link{mak}}, \code{\link{plk}}, etc.).
#' @param op An object of class \code{\link{opt}} (created via \code{\link{opt}}). Compatibility check with \code{o} is conducted.
#' @param x Matrix of dimension mx(d+1) containing custom time and smoothed values of the processes used for the integral matching, see details. If \code{NULL} (default) a linear-interpolation smoother is employed.
#' @param adapts Character vector holding names of what quantities to adapt during algorithm. Possible quantities are: \code{"scales"}, \code{"weights"} and \code{"penalty_factors"}, see details. Default is \code{"scales"} and \code{"penalty_factor"}.
#' @param xout Logical indicating if a matrix containing the process values at the time points in \code{y} in \code{op}, linearly interpolated from \code{x} should be returned. Default is \code{FALSE}.
#' @param params List of vectors of initial parameter values. Default (\code{NULL}) translates to the zero-vector. If the ODE system is linear in the parameter vector and the boundary constraints on are not unusual, then \code{params} is ignored.
#' @param ... Additional arguments passed to \code{aim}.
#' @seealso rodeo, rodeo.ode, rodeo.aim, imd
#' @details
#' \subsection{Loss function}{
#' Integral matching requires a smoothed process in order to approximate the parameter estimates. More precisely, the loss function
#' \deqn{RSS / (2 * (n - s)) + lambda*penalty}
#' is minimised, where RSS is the sum of the squared 2-norms of
#' \deqn{x_i - x_{i-1} - \int_{t_{i-1}}^{t_{i}}{f(x(s), context_l * param) ds}}
#' Here f is the vector field of the ODE-system, x is the smoothed process and param is (internally) scaled with \code{scales} in \code{reg}.
#' }
#' \subsection{Custom smoother}{
#' The supplied \code{x} is a way of customising how x in the loss function is made. Firstly, \code{x} must have similiar layout as \code{y} in \code{op}, i.e., first column is time and the remaining columns contain smoothed values of the process at those time points, see \code{\link{opt}}-documentation for details.
#'
#' The number of decreases in time in \code{x} must match the number of decreases in time in \code{y} in \code{op}. The process x in the loss function is simply a linear interpolation of \code{x}, hence finer discretisations give more refined integral matching. Each context is handled seperately. Moreover, the compatibility checks in \code{rodeo} are also conducted here.
#' }
#' \subsection{Adaptations}{
#' The \code{adapts} are ways to adapt quantities in the \code{\link{opt}}-object and \code{\link{reg}}-objects in \code{o} to the data:
#' \describe{
#' \item{\code{"scales"}}{A standardisation of the columns in the linearisation takes place and carry over to \code{scales} in \code{\link{reg}}-objects in \code{o} (generally recommended and default).}
#' \item{\code{"weights"}}{The observational weights (\code{weights} in \code{op}) are adjusted coordinate-for-coordinate (column-wise) by reciprocal average residual sum of squares across penalty parameters.}
#' \item{\code{"penalty_factor"}}{The penalty factors in \code{\link{reg}} are adjusted by the reciprocal average magnitude of the estimates (parameters whose average magnitude is 0 join \code{exclude}). For extremely large systems, this option can heavily reduce further computations if the returned \code{\link{aim}}-object is passed to \code{\link{rodeo}}.}
#' }
#' If either \code{"penalty_factor"} or \code{"weights"} are in \code{adapts} a refitting takes place.
#' }
#' Finally, note that \code{lambda} and \code{nlambda} in returned \code{\link{opt}}-object may have changed.
#' @return An object with S3 class "aim":
#' \item{o}{Original \code{\link{ode}}-object with adapted quantities.}
#' \item{op}{Original \code{\link{opt}}-object with adapted quantities and default values for \code{lambda_min_ratio}, \code{lambda} and (if needed) \code{a} inserted, if these were originally \code{NULL}.}
#' \item{params}{A list of matrices (of type dgCMatrix) with the parameter estimates (scaled by \code{scales} in \code{reg}), one column for each lambda value. One element in the list for each parameter (other than initial state parameter).}
#' \item{x0s}{A matrix with the initial state estimates.}
#' \item{rss}{A vector of residual sum of squares.}
#' \item{x}{Original \code{x}, or \code{y} in \code{op} if \code{x} was \code{NULL}.}
#' \item{xout}{(If \code{xout} = \code{TRUE}) A matrix containing the values of the process at the time points in \code{y} in \code{op}, interpolated from \code{x}. Use it to check that \code{x} works as intended.}
#'
#'
#' @examples
#' set.seed(123)
#' # Michaelis-Menten system with two 0-rate reactions
#' A <- matrix(c(1, 1, 0, 0,
#'   0, 0, 1, 0,
#'   0, 0, 1, 0,
#'   0, 0, 1, 1,
#'   0, 0, 1, 1), ncol = 4, byrow = TRUE)
#' B <- matrix(c(0, 0, 1, 0,
#'   1, 1, 0, 0,
#'   1, 0, 0, 1,
#'   0, 1, 0, 0,
#'   1, 0, 0, 1), ncol = 4, byrow = TRUE)
#' k <- c(0.1, 1.25, 0.5, 0, 0); x0 <- c(E = 5, S = 5, ES = 1.5, P = 1.5)
#' Time <- seq(0, 10, by = 1)
#'
#' # Simulate data, in second context the catalytic rate has been inhibited
#' contexts <- cbind(1, c(1, 1, 0, 1, 1))
#' m <- mak(A, B, r = reg(contexts = contexts))
#' y <- numsolve(m, c(Time, Time), cbind(x0, x0 + c(2, -2, 0, 0)), contexts * k)
#' y[, -1] <- y[, -1] + matrix(rnorm(prod(dim(y[, -1])), sd = .25), nrow = nrow(y))
#'
#' # Create optimisation object via data
#' o <- opt(y, nlambda = 10)
#'
#' # Fit data using Adaptive Integral Matching on mak-object
#' a <- aim(m, o)
#' a$params$rate
#'
#' # Compare with true parameter on column vector form
#' matrix(k, ncol = 1)
#'
#'
#' # Example: Power Law Kinetics
#' A <- matrix(c(1, 0, 1,
#'   1, 1, 0), byrow = TRUE, nrow = 2)
#' p <- plk(A)
#' x0 <- c(10, 4, 1)
#' theta <- matrix(c(0, -0.25,
#'   0.75, 0,
#'   0, -0.1), byrow = TRUE, nrow = 3)
#' Time <- seq(0, 1, by = .025)
#'
#' # Simulate data
#' y <- numsolve(p, Time, x0, theta)
#' y[, -1] <- y[, -1] + matrix(rnorm(prod(dim(y[, -1])), sd = .25), nrow = nrow(y))
#'
#' # Estimation
#' a <- aim(p, opt(y, nlambda = 10))
#' a$params$theta
#'
#' # Compare with true parameter on column vector form
#' matrix(theta, ncol = 1)
#'
#'
#' # Example: use custom lowess smoother on data
#' # Smooth each coordinate of data to get curve
#' # on extended time grid
#' ext_time <- seq(0, 1, by = 0.001)
#' x_smooth <- apply(y[, -1], 2, function(z) {
#'   # Create loess object
#'   data <- data.frame(Time = y[, -1], obs = z)
#'   lo <- loess(obs ~ Time, data)
#'
#'   # Get smoothed curve on extended time vector
#'   predict(lo, newdata = data.frame(Time = ext_time))
#' })
#'
#' # Bind the extended time
#' x_smooth <- cbind(ext_time, x_smooth)
#'
#' # Run aim on the custom smoothed curve
#' a_custom <- aim(p, opt(y), x = x_smooth)
#'
#'
#' @export
#' @useDynLib episode, .registration = TRUE
#' @importFrom Rcpp sourceCpp
aim <- function(o, op, x = NULL, adapts = c("scales", "penalty_factor"), xout = FALSE, params = NULL, ...) {
  ### Adjustments
  if (!all(adapts %in% c("scales", "weights", "penalty_factor"))) warning("Not all adapts are among the options: 'scales', 'weights', 'penalty_factor'.")
  adapt_scales   <- isTRUE("scales" %in% adapts)
  adapt_weights  <- isTRUE("weights" %in% adapts)
  adapt_penalty  <- isTRUE("penalty_factor" %in% adapts)

  ### Do compatibility check of dimension between ode and ro
  .o_op_check(o = o, op = op)
  d <- o$d
  ps <- lapply(o$rs[-1], getElement, "p")
  n <- nrow(op$y)
  s <- op$s
  x0_fixed <- o$rs[[1]]$fixed
  o$rs[[1]]$fixed <- TRUE

  ### x_smooth NULL handling or checking ###
  x <- .prepare_x(x, op)

  ### Set defaults
  op <- .set_defaults(o = o, op = op)


  ## For standard cases of aimfit
  alpha <- lapply(o$rs, function(re) {
    if (re$reg_type == "l1") {
      1
    } else if (re$reg_type == "l2") {
      0
    } else if (re$reg_type == "elnet") {
      re$a
    } else {
      1
    }
  })


  ### Actual AIM part ###
  # note: it receives params on R scale, but .get_aimfit puts it to raw cpp scale
  aimfit <- .get_aimfit(o = o, op = op, x = x, alpha = alpha, standardize = adapt_scales, params = params, ...)
  des <- aimfit$des
  aimfit <- aimfit$aimfit
  # but des returned is X(eval in context * scales) * diag(context * scales)

  ## Scales
  if (adapt_scales) {
    o$rs[-1] <- Map(function(re, X) {
      colnorms <- sqrt(Matrix::colSums(X^2))
      if (mean(colnorms) == 0) colnorms <- rep(1, ncol(X))
      colnorms <- colnorms / mean(colnorms)
      colnorms[colnorms == 0] <- 1e-6

      # Since scales are allready in des$X, it carries over
      re$scales <- colnorms

      re
    }, re = o$rs[-1], X = des$X)
  }


  ## Weights and penalty_factor
  if (adapt_penalty | adapt_weights) {
    # Select model (gcv based ensemble)
    fitt <- Reduce("+", aimfit$fits)
    if (is.null(des$W)) {
      rss_fit <- apply((des$Y - fitt)^2, 2, function(x) apply(matrix(x, nrow = d), 1, sum))
    } else {
      rss_fit <- apply(des$W * (des$Y - fitt)^2, 2, function(x) apply(matrix(x, nrow = d), 1, sum))
    }
    dfs_fit <- colSums(aimfit$dfs)
    # gcv-weighted ensemble
    rec_gcv <- 1 / (apply(rss_fit, 2, sum) / (1 - pmin(1, dfs_fit / (d * n)))^2) # reciprok gcv

    # Weights
    if (adapt_weights) {
      w_const <- 1 / rss_fit %*% (rec_gcv / (pmax(1, d * (n - s) - dfs_fit)))
      if (is.null(op$weights)) {
        op$weights <- matrix(1, nrow = nrow(op$y), ncol = d)
      }
      row_means <- rowMeans(op$weights) # stored and used so that they will remain unchanged
      w_const <- as.vector(w_const) / mean(w_const)
      op$weights <- t(t(op$weights) * w_const)
      op$weights <- op$weights * row_means / rowMeans(op$weights)
    }

    # Penalty factor
    if (adapt_penalty) {
      v_const <- lapply(aimfit$beta, function(beta) abs(as.vector(beta %*% (rec_gcv / mean(rec_gcv)))))

      o$rs[-1] <- Map(function(re, v) {
        if (mean(v) > 0) {
          v <- v / mean(v)

          if (!is.null(re$penalty_factor)) {
            re$penalty_factor[v > 0] <- re$penalty_factor[v > 0] / v[v > 0]
            re$penalty_factor <- re$penalty_factor / mean(re$penalty_factor)
          } else {
            re$penalty_factor <- rep(1, length(v))
            re$penalty_factor[v > 0] <- 1 / v[v > 0]
          }

          # Non-positive penalties:
          if (methods::is(o, "mak")) {
            re$exclude <- union(re$exclude, which(v <= 0))
          } else {
            re$penalty_factor[v <= 0] <- 1e6
          }
        }
        re
      }, re = o$rs[-1], v = v_const)
    }

    # Refit
    aimfit <- .get_aimfit(o = o, op = op, x = x, alpha = alpha, standardize = adapt_scales, params = params, ...)
    # print(aimfit)
    des <- aimfit$des
    aimfit <- aimfit$aimfit
  }


  ### rss
  fitt <- Reduce("+", aimfit$fits)
  if (is.null(des$W)) {
    rss <- apply((des$Y - fitt)^2, 2, function(x) apply(matrix(x, nrow = d), 1, sum))
  } else {
    rss <- apply(des$W * (des$Y - fitt)^2, 2, function(x) apply(matrix(x, nrow = d), 1, sum))
  }


  ### Overwrite lambda
  op$lambda  <- aimfit$lambda
  op$nlambda <- length(op$lambda)
  if (!is.null(aimfit$lambda_factors)) {
    o$rs <- Map(function(re, lf) {
      re$lambda_factor <- lf
      re
    }, o$rs, lf = as.list(aimfit$lambda_factors))
  }

  o$rs[[1]]$fixed <- x0_fixed
  ret <- structure(list(o = o, op = op,
    params = lapply(aimfit$beta, methods::as, "dgCMatrix"),
    x0s = matrix(des$X0, nrow = length(des$X0), ncol = op$nlambda),
    #fitted = aimfit$fits,
    dfs = aimfit$dfs,
    rss = rss,
    x = x),
    class = "aim")
  if (xout) {
    ret$xout <- des$xout
  }
  if (!is.null(params) & !is.null(names(params))) {
    names(ret$params) <- names(params)
  } else {
    names(ret$params) <- names(o$rs[-1])
  }
  ret$jerr <- aimfit$jerr   # error flags, etc.
  if (all(unlist(lapply(o$rs[-1], getElement, "reg_type")) == "none")) {
    ret$params  <- lapply(ret$params, function(para) para[, which.min(ret$op$lambda), drop = FALSE])
    ret$x0s <- ret$x0s[, which.min(ret$op$lambda), drop = FALSE]
    ret$rss <- ret$rss[, which.min(ret$op$lambda), drop = FALSE]
    ret$op$lambda  <- 1
    ret$op$nlambda <- 1
  }
  return(ret)
}





#' Regularised Ordinary Differential Equation Optimisation (RODEO) initialised via Adaptive Integral Matching
#'
#' @description Fit the parameters for an ODE model with data sampled across different contexts.
#' @param x \code{aim}-object created via \code{\link{aim}}-function.
#' @param adjusts Character vector holding names of what quantities to adjust during algorithm. Possible quantities are: \code{"lambda"}, \code{"scales"} and \code{"weights"}.
#' @param trace Logical indicating if status messages should be printed during \code{rodeo}.
#' @param ... Additional arguments passed to \code{rodeo}.
#' @seealso rodeo, aim, rodeo.ode
#' @inherit rodeo return
#'
#' @details
#' The adapted quantities (\code{scales}, \code{weights}, \code{penalty_factor}) of \code{x} (returned by \code{\link{aim}}) are fed to the exact estimator \code{rodeo}. This estimator then traverses the \code{lambda} sequence in reverse order initialised in the last estimates from \code{aim}.
#'
#' If desired, the quantities \code{lambda}, \code{scales} and \code{weights} are adjusted as in \code{aim}.
#'
#' @examples
#' set.seed(123)
#' # Example: Power Law Kinetics
#' A <- matrix(c(1, 0, 1,
#'   1, 1, 0), byrow = TRUE, nrow = 2)
#' p <- plk(A)
#' x0 <- c(10, 4, 1)
#' theta <- matrix(c(0, -0.25,
#'   0.75, 0,
#'   0, -0.1), byrow = TRUE, nrow = 3)
#' Time <- seq(0, 1, by = .025)
#'
#' # Simulate data
#' y <- numsolve(p, Time, x0, theta)
#' y[, -1] <- y[, -1] + matrix(rnorm(prod(dim(y[, -1])), sd = .25), nrow = nrow(y))
#'
#' # Estimation via aim
#' a <- aim(p, opt(y, nlambda = 10))
#' a$params$theta
#'
#' # Supply to rodeo
#' rod <- rodeo(a)
#' rod$params$theta
#'
#' # Compare with true parameter on column vector form
#' matrix(theta, ncol = 1)
#'
#'
#' # Example: include data from an intervened system
#' # where the first complex in A is inhibited
#' contexts <- cbind(1, c(0, 0, 0, 1, 1, 1))
#' y2 <- numsolve(p, Time, x0 + 1, theta * contexts[, 2])
#' y2[, -1] <- y2[, -1] + matrix(rnorm(prod(dim(y2[, -1])), sd = .25), nrow = nrow(y2))
#'
#' # Estimation via aim
#' a <- aim(plk(A, r = reg(contexts = contexts)), opt(rbind(y, y2), nlambda = 10))
#' a$params$theta
#'
#' # Supply to rodeo
#' rod <- rodeo(a)
#' rod$params$theta
#'
#'
#' @export
#' @useDynLib episode, .registration = TRUE
#' @importFrom Rcpp sourceCpp
rodeo.aim <- function(x, adjusts = c("lambda"), trace = FALSE, ...) {
  stopifnot(methods::is(x, "aim"))
  .check(x)

  ### Adjustments
  if (!all(adjusts %in% c("lambda", "scales", "weights"))) warning("Not all adjusts are among the options: 'lambda', 'scales', 'weights'.")
  adjust_scales   <- isTRUE("scales" %in% adjusts)
  adjust_weights  <- isTRUE("weights" %in% adjusts)
  adjust_lambda   <- isTRUE("lambda" %in% adjusts)


  ## Prepare for rodeo ##
  s <- x$op$s
  d <- x$o$d
  ps <- lapply(x$o$rs[-1], getElement, "p")
  nlams <- unlist(lapply(x$params, ncol))
  if (diff(range(nlams)) != 0) stop("Number of columns in parameter arguments are not equal.")
  indices <- nlams[1]
  indices <- unique(sort(indices))
  if (length(indices) == 0) stop("Indices has length 0.")
  if (any(!(indices %in% seq_len(nlams[1])))) stop(paste0("Some indices are not among 1:", nlams[1], " (seq_along(lambda))."))

  # Remove scales (cpp part operates on this) and subset the indices
  pnames <- names(x$params)
  x$params  <- Map(function(para, re) {
    para <- para[, indices, drop = FALSE]
    if (!is.null(re$scales)) {
      para <- para / re$scales
    }
    as.list(data.frame(as.matrix(para))) # turn into list of vectors
  }, para = x$params, re = x$o$rs[-1])

  # Switch hierarchy, so that it is: particles, parameter argument
  x$params <- .switch_hierarchy(x$params)

  # Arrayify x0
  x0 <- x$x0s[, indices, drop = FALSE]
  x0 <- array(x0, dim = c(d, s, ncol(x0)))

  # scales and contexts
  sc <- lapply(x$o$rs[-1], .get_sc, op = x$op)

  # Exclude handle
  x$params <- lapply(x$params, function(particle) {
    Map(function(para, exc) {
      if (!is.null(exc) & length(exc) > 0) para <- para[-exc]
      return(para)
    }, para = particle, exc = lapply(x$o$rs[-1], getElement, "exclude"))
  })
  sc_ <- Map(function(sc, exc) {
    if (!is.null(sc) & !is.null(exc) & length(exc) > 0) sc <- sc[-exc, , drop = FALSE]
    return(sc)
  }, sc = sc, exc = lapply(x$o$rs[-1], getElement, "exclude"))
  ode_struct <- x$o
  ode_struct <- .exclude(ode_struct)


  ## Actual fitting
  broncfit <- bronc(ode_struct = ode_struct, opt_struct = x$op, sc_ = sc_,
                    params_ = x$params, x0s = x0,
                    indices = indices,
                    adjust_lambda = adjust_lambda,
                    adjust_scales = adjust_scales,
                    adjust_weights = adjust_weights,
                    trace = isTRUE(trace))


  ## Put adjustments back into the return object
  if (adjust_scales) {
    x$o$rs[-1] <- Map(function(re, exc, new_sc) {
      new_sc <- as.vector(new_sc)
      if (is.null(re$scales)) re$scales <- rep(1, re$p)
      if (is.null(exc) | length(exc) == 0) {
        re$scales <- re$scales * new_sc * mean(re$scales) / mean(re$scales * new_sc)
      } else {
        re$scales[-exc] <- re$scales[-exc] * new_sc * mean(re$scales[-exc]) / mean(re$scales[-exc] * new_sc)
      }
      return(re)
    }, re = x$o$rs[-1], exc = lapply(x$o$rs[-1], getElement, "exclude"), new_sc = broncfit$new_scales)
  }
  if (adjust_weights) {
    x$op$weights <- broncfit$new_weights
  }
  if (adjust_lambda) {
    x$op$lambda <- as.vector(broncfit$lambda)
  }

  ## lambda_factor needs to be replaced in either case,
  x$o$rs <- Map(function(re, lf) {
    re$lambda_factor <- lf
    re
  }, re = x$o$rs, lf = as.list(broncfit$lambda_fac))


  ## Hierarchy becomes: parameter, particle
  broncfit$params <- .switch_hierarchy(broncfit$params)

  ## correct scales and add excluded variables
  x0s <- array(do.call(c, lapply(broncfit$params[[1]], function(x) as.vector(x))),
    dim = c(nrow(broncfit$params[[1]][[1]]), ncol(broncfit$params[[1]][[1]]), length(broncfit$params[[1]])))
  broncfit$params <- Map(function(para, re) {
    para <- lapply(para, .add_exclude, exclude = re$exclude)
    if (!is.null(re$scales)) {

      para <- lapply(para, function(x, y) {y * x}, y = re$scales)
    }
    para
  }, para = broncfit$params[-1], re = x$o$rs[-1])
  broncfit$params <- stats::setNames(broncfit$params, pnames)

  ## return
  multi <- dim(x0s)[3] > 1
  if (multi) {
    ret <- structure(list(o = x$o, op = x$op, params = broncfit$params, x0s = x0s,
      dfs = broncfit$dfs, codes = broncfit$codes, steps = broncfit$steps,
      losses = broncfit$losses, penalties = broncfit$penalties,
      jerr = broncfit$jerr), class = "rodeo")
  } else {
    ret <- structure(list(o = x$o, op = x$op,
      params = lapply(broncfit$params, function(para) para[[1]]), x0s = x0s[,,1],
      dfs = broncfit$dfs[,,1], codes = broncfit$codes[,,1], steps = broncfit$steps[,,1],
      losses = broncfit$losses[1,], penalties = broncfit$penalties[,,1],
      jerr = broncfit$jerr[,,1]), class = "rodeo")
  }
  return(ret)
}




#' Estimation with Penalisation In Systems of Ordinary Differential Equations.
#'
#' This package provide the tools for approximate and exact parameter estimation in ODE models with regularisation via penalisation.
#'
#' @section Specify your ODE:
#' The ode system is specified via the \code{\link{ode}}-subclasses: \code{\link{mak}}, \code{\link{plk}}, \code{\link{rlk}} and \code{\link{ratmak}}. In creating these you can also specify numerical solver type via \code{\link{solver}} and regularisation type of the parameter arguments via \code{\link{reg}}. To numerically solve the ODE use \code{\link{numsolve}} and to evaluate the ODE field use \code{\link{field}}. The differentials of both quantities can also be evaluated.
#'
#' @section Specify loss function:
#' To specify the loss function use the \code{\link{reg}} in the ODE object to control regularisation and \code{\link{opt}} to control the observations, their weights and the tuning parameter of the regularisation.
#'
#' @section Optimise the loss function:
#' Having an \code{\link{ode}} object and an \code{\link{opt}} object, there are two methods for estimating the parameters: approximate estimation via inverse collocation methods, \code{\link{aim}}, or exact estimation via interior point methods, \code{\link{rodeo}}. If desired, call \code{\link{rodeo}} on the results from \code{\link{aim}} to use the approximate estimates for initialising the exact estimation.
#'
#' @docType package
#' @name episode
NULL




# #' Regularised Optimisation of Linear Least Squares
# #'
# #' @description Solves penalised least squares problems, using the interior point method employed in this package. This is not a computationally efficient algorithm, hence if other packages can solve your problem (e.g., 'glmnet' or 'ncvreg') consider using those instead. The main reason for using this solver, is that it can handle some cases that classic solvers do not cover, e.g., box-constrained SCAD and MCP.
# #' \subsection{Loss function:}{
# #' The loss function
# #' \deqn{RSS / (2 * (n - 1)) + lambda*penalty}
# #' is minimised, where RSS is the sum of the squared 2-norms of
# #' \deqn{y - X * param}
# #' where y is a n-vector, X is a nxp matrix and param is p-vector.
# #'
# #' Note: no centering or standardisation takes place!
# #' }
# #'\subsection{Regularisation:}{
# #'The type of regularisation is chosen via reg_type:
# #'\describe{
# #'\item{\code{l1}:}{Least Absolute Shrinkage and Selection Operator (Lasso) penalty. The penalty is the absolute value: \eqn{pen(param_j)=|param_j|}}
# #'\item{\code{l2}:}{Ridge penalty. The penalty is the squaring: \eqn{pen(param_j)=param_j^2/2}}
# #'\item{\code{elnet}:}{Elastic Net. A convex combination of \code{l1} and \code{l2} penalties: \eqn{pen(param_j)=(1-a)param_j^2/2+a|param_j|}, 0<=a<=1. Note if a = 0 or a = 1, then the penalty is automatically reduced to "l2" and "l1" respectively.}
# #'\item{\code{scad}:}{Smoothly Clipped Absolute Deviation penalty: \deqn{pen(param_j)=\int^{param_j}_{0}{min( max((a\lambda-|param|)/(\lambda(a-1)), 0), 1) dparam}, a > 2.}}
# #'\item{\code{mcp}:}{Maximum Concave Penalty penalty: \deqn{pen(param_j)=\int^{param_j}_{0}{max(1 - |param|/(\lambda a), 0) dparam}, a > 1.}}
# #'\item{\code{none}:}{No penalty. Not recommended for large systems. This overwrites both user-supplied and automatically generated lambda-sequences.}
# #'}
# #'}
# #' @param X Design matrix of dimension nxp.
# #' @param y Observation vector of length n.
# #' @param w (Optional) observation weights vector of length n.
# #' @param reg_type Character string determining the regularisation. Must be one of: "l1" (default), "l2", "elnet", "scad", "mcp" or "none". See details.
# #' @param a Numeric value of tuning parameter in \code{elnet} (must be between 0 and 1, default is 0.5), \code{scad} (must be larger than 2, default = 3.7) or \code{mcp} (must be larger than 1, default = 3).
# #' @param nlambda Number of lambda values.
# #' @param lambda_min_ratio Ratio between smallest and largest value of lambda (latter derived using data). If n  > p, the default is 0.0001, else 0.01.
# #' @param lambda A custom lambda sequence. If not \code{NULL}, \code{nlambda} and \code{lambda_min_ratio} are ignored. The optimisation reuses last optimal parameter value for next step, so \code{lambda} should preferably be monotonic.
# #' @param lower_param Either numeric vector of length 1 or p, of lower limit(s) for parameters. Defaults to -Inf.
# #' @param upper_param Either numeric vector of length 1 or p, of upper limit(s) for parameters. Defaults to Inf.
# #' @param exclude A vector indices to exclude from model (this is how to specify an infinite penalty factor). Default is none.
# #' @param penalty_factor A non-negative vector (p) of individual penalty weights. Defaults to 1 for each parameter. Will always be rescaled so its mean is 1.
# #' @param screen Logical indicating if a faster optimisation relying on screening should be adapted. If \code{NULL}, it i set to \code{TRUE} iff p > 20.
# #' @param step_max Positive integer giving the maximal number of steps in optimisation procedure, per lambda.
# #' @param step_screen Positive integer giving the number of steps between screenings.
# #' @param backtrack_max Positive integer giving the maximal number of backtracking steps taken in each optimisation step.
# #' @param tol_l Positive numeric tolerance level used for stopping criterion (max-norm of gradients).
# #' @param tol_param Positive numeric tolerance level used for parameter space.
# #' @param tau_init Positive initial step length for backtracking.
# #' @param tau_min Positive numeric giving minimal value of step length in backtracking.
# #' @param tau_scale Scaling parameter of step length, must be in (0,1).
# #' @param trace Logical indicating if status messages should be printed during \code{rolls}.
# #' @param ... Additional arguments passed to \code{rolls}.
# #'
# #' @return A list holding:
# #' \item{params}{A matrix (of type dgCMatrix) with the parameter estimates, one column for each lambda value.}
# #' \item{dfs}{A vector of degrees of freedom estimates.}
# #' \item{lambda}{A vector of tuning parameter values.}
# #' \item{codes}{A vector of convergence codes for each lambda.
# #' \describe{
# #' \item{0:}{The convergence criteria is met (see details in \code{\link{ro}}). Current estimate is probably a local minimum.}
# #' \item{1:}{Backtracking in the last iteration yields no numerical improvement, but no unsual behavior observed. Current estimate is probably a local minimum. However, if \code{exact_gradient} = 4, changing this may improve the code. Alternatively one can adjust backtracking via \code{backtrack_max} and \code{tau_min}.}
# #' \item{2:}{The optimisation procedure exceeded maximal number of steps.}
# #' \item{3:}{The last gradient was unsually large. The tolerances may be off.}
# #'}}
# #' \item{steps}{A vector holding number of steps used in optimisation procedure for each lambda.}
# #' \item{losses}{A vector of penalised losses at optimum for each lambda value.}
# #' \item{rss}{A vector of residual sum of squares at optimum for each lambda value.}
# #' \item{jerr}{Summary codes (for internal debugging).}
# #'
# #' @examples
# #' X <- matrix(rnorm(10 * 12), nrow = 10)
# #' y <- as.vector(X %*% c(5, 2, 1, .5, rep(0, 8))) + rnorm(10, sd = .5)
# #'
# #' # Lasso
# #' rolls(X, y)
# #'
# #' # Ridge regression
# #' rolls(X, y, reg_type = "l2")
# #'
# #' # SCAD with positive coefficients
# #' rolls(X, y, reg_type = "scad", lower_param = 0)
# #'
# #' @export
# rolls <- function(X, y, w = NULL, reg_type = "l1", a = NULL, nlambda = 100,
#   lambda_min_ratio = ifelse(n < p, 0.01, 0.0001), lambda = NULL,
#   lower_param = -Inf, upper_param = Inf,
#   exclude = NULL, penalty_factor = NULL, screen = (p > 20),
#   step_max = 10000, step_screen = 100, backtrack_max = 50,
#   tol_l = 1e-7, tol_param = 1e-7, tau_init = 0.1, tau_min = 1e-12, tau_scale = .5,
#   trace = FALSE, ...) {
#
#   # Penalty factor rescale
#   if (!is.null(penalty_factor)) {
#     m <- mean(penalty_factor)
#     if (m <= 0 || !is.finite(m)) {
#       stop("'penalty_factor' is not null, but its mean is non-positive or infinite.")
#     }
#     penalty_factor <- penalty_factor / m
#   }
#
#   # First n checks
#   n <- length(y)
#   if (!is.null(w)) if (length(w) != n) stop("Length of y and w do not match.")
#   if (nrow(X) != n) stop("Number of rows in X does not match length of y.")
#
#   # p
#   p <- ncol(X)
#   exclude <- unique(exclude)
#
#   # Call barrel_rolls
#   o <- plk(A = matrix(0, nrow = p, ncol = 1), lower_param = lower_param, upper_param = upper_param)
#   r <- ro(cbind(1, y), reg_type = reg_type, a = a, nlambda = nlambda, lambda_min_ratio = lambda_min_ratio,
#     lambda = lambda, exclude = exclude, penalty_factor = penalty_factor, screen = screen,
#     step_max = step_max, step_screen = step_screen, backtrack_max = backtrack_max,
#     tol_l = tol_l, tol_param = tol_param, tol_x = tol_param,
#     tau_init = tau_init, tau_min = tau_min, tau_scale = tau_scale)
#   r <- .set_defaults(o, r)#, matrix(0, 1, 1), set_lambda = FALSE)
#   if (is.null(r$lambda)) {
#     # Determine a
#     if (r$reg_type == "l2") {
#       alpha <- 0.001
#     } else if (r$reg_type == "elnet") {
#       alpha <- max(0.001, r$a)
#     } else {
#       alpha <- 1
#     }
#
#     # Get SS diff at 0
#     if (is.null(w)) {
#       diffat0 <- - t(y) %*% X / (n - 1)
#     } else {
#       diffat0 <- - t(y * w) %*% X / (n - 1)
#     }
#
#     # Remove not useful ones
#     useful <- is.finite(diffat0) & ((diffat0 > 0 & o$lower_param < 0) | (diffat0 < 0 & o$upper_param > 0))
#     if (is.null(r$penalty_factor)) {
#       lam_max <- max(abs(diffat0[useful] / alpha))
#     } else {
#       useful <- useful & (r$penalty_factor > 0)
#       lam_max <- max(abs(diffat0[useful] / (alpha * r$penalty_factor[useful])))
#     }
#
#     r$lambda  <- max(lam_max, 1e-7) * exp(seq(0, log(r$lambda_min_ratio), length.out = r$nlambda))
#       exp(seq(log(lam_max), log(lam_max) + log(r$lambda_min_ratio), length.out = r$nlambda))
#   }
#
#   if (length(exclude) > 0) {
#     o <- exclude(o, r$exclude)
#     r <- exclude(r, r$exclude)
#     X <- X[, -r$exclude, drop = FALSE]
#   }
#   res_rolls <- barrel_rolls(ode_struct = o, ro_struct = r, X = X, Y = y, W = w, trace = trace)
#   if (length(exclude) > 0) {
#     res_rolls$params <- .add_exclude(res_rolls$params, exclude)
#   }
#
#   return(res_rolls)
# }
#

