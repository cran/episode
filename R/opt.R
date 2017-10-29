###-----------------------###
###  Optimisation object  ###
###-----------------------###


#' Create 'opt' (optimisation) object
#'
#' This function creates an object of class \code{opt}, which holds data, weights, tuning parameters and control list for optimisation. This is basically the control panel for the loss function optimised in \code{\link{rodeo}} and \code{\link{aim}}.
#'
#'
#' @param y Numeric matrix (nx(d+1)). First column is time and remaining columns hold observations from different species, see details. Missing values (except in first column) are allowed and are marked with \code{NA} or \code{NaN}.
#' @param weights Numeric matrix (nxd) of observation weights (optional).
#' @param nlambda Number of lambda values.
#' @param lambda_min_ratio Ratio between smallest and largest value of lambda (latter derived using data). If (n - s) * d > p, the default is 0.0001, else 0.01.
#' @param lambda_decrease Logical indicating if automatically generated lambda sequence should be decreasing (default is \code{TRUE}). Consider switching to \code{FALSE} if nonzero initialisations are given to \code{rodeo.ode}.
#' @param lambda A custom lambda sequence. If not \code{NULL}, \code{nlambda}, \code{lambda_min_ratio} and \code{decrease.lambda} are ignored. The optimisation reuses last optimal parameter value for next step, so \code{lambda} should preferably be monotonic.
#' @param tol_l Positive numeric tolerance level used for stopping criterion (max-norm of gradients).
#' @seealso aim, rodeo, rodeo.aim, rodeo.ode
#' @inherit reg details
#' @details
#' \subsection{Data format}{
#' Whenever time (first column of \code{y}) decreases, the system restarts. Hence the data is assumed generated from s different contexts, where s - 1 is the number of decreases in the time vector.
#'
#' Each context has its own initial condition and parameter vector specified through \code{contexts} (see \code{\link{reg}} for details).
#' }
#' \subsection{Loss function}{
#' The loss function optimised in \code{\link{rodeo}} is:
#' \deqn{RSS/(2 * (n - s)) + \lambda*\sum_{parameter \ argument}\lambda_{factor}*\sum_{j=1}^p v_jpen(param_j)}
#' where \eqn{\lambda} belongs to the \code{lambda}-sequence and v is \code{penalty_factor}. Moreover, the residual sum of squares, RSS, is given as:
#' \deqn{RSS = \sum^{n}_{i=1}||(y(t_i) - x(t_i, {x_0}_l, context_l * param))*\sqrt{w(t_i)}||_2^2}
#' where \code{param} has been (internally) scaled with \code{scales}, and \eqn{w(t_i)} and \eqn{y(t_i)} refers to the i'th row of \code{weights} and \code{y} (with first column removed), respectively. The solution to the ODE system is the x()-function. The subscript 'l' refers to the context, i.e., the columns of \code{contexts} in \code{reg}-object and \code{x0} in \code{rodeo}-functions (\code{x0} is the initial state of the system at the first time point after a decrease in the time-vector). All products are Hadamard products.
#' }
#' \subsection{Tuning parameter}{
#' The lambda sequence can either be fully customised through \code{lambda} or automatically generated. In the former case, a monotonic \code{lambda}-sequence is highly recommended. Throughout all optimisations, each optimisation step re-uses the old optimum, when sweeping through the \code{lambda}-sequence.
#'
#' If \code{lambda} is \code{NULL}, an automatically generated sequence is used. A maximal value of lambda (the smallest at which 0 is a optimum in the rate parameter) is calculated and log-equi-distant sequence of length \code{nlambda} is generated, with the ratio between the smallest and largest lambda value being \code{lambda_min_ratio}. Note: when the \code{opt}-object is passed to \code{\link{rodeo.ode}}, one may want to initialise the optimisation at a non-zero parameter and run the optimisation on the reversed lambda sequence. This is indicated by setting \code{decrease.lambda} = FALSE. If, however, the \code{opt}-object is passed to \code{\link{aim}}, \code{\link{glmnet}} ultimately decides on the lambda sequence, and may cut it short.
#' }
#' \subsection{Optimisation}{
#'  A proximal-gradient type of optimisation is employed. The step length is denoted \eqn{\tau}. The convergence criteria is \eqn{\Delta \eta < 10^3 * \max(|(\Delta param \neq 0)| + |(\Delta x0 \neq 0)|, 1)} AND \eqn{\Delta loss / \Delta \eta < tol_l}, where
#'  \deqn{\Delta \eta = ||\Delta param||/tol_{param} + ||\Delta x0|| / tol_{x_0}}
#'
#' }
#'
#' @return An object with S3 class "opt".
#'
#' @examples
#'
#' # Generate some data
#' set.seed(123)
#' A <- matrix(c(1, 0, 1,
#'               1, 1, 0), byrow = TRUE, nrow = 2)
#' p <- plk(A)
#' x0 <- c(10, 4, 1)
#' theta <- matrix(c(0, -0.25,
#'                   0.75, 0,
#'                   0, -0.1), byrow = TRUE, nrow = 3)
#' Time <- seq(0, 1, by = .025)
#' y <- numsolve(p, Time, x0, theta)
#' y[, -1] <- y[, -1] + matrix(rnorm(prod(dim(y[, -1])), sd = .1), nrow = nrow(y))
#'
#' # Minimally, you need to supply data:
#' op <- opt(y)
#'
#' # More weight on early observations
#' w <- outer(1 / seq_len(nrow(y)), rep(1, length(x0)))
#' op <- opt(y, weights = w)
#'
#' # Less weight on first coordinate
#' w <- outer(rep(1, nrow(y)), c(1, 2, 2))
#' op <- opt(y, weights = w)
#'
#' # Adjust tuning parameter sequence
#' op <- opt(y, nlambda = 10, lambda_min_ratio = 0.05)
#'
#' @export
opt <- function(y, weights = NULL, nlambda = 25, lambda_min_ratio = NULL, lambda_decrease = TRUE, lambda = NULL, tol_l = 1e-5) {
  ret <- structure(list(y = y, weights = weights, nlambda = nlambda,
    lambda_min_ratio = lambda_min_ratio, lambda_decrease = isTRUE(lambda_decrease), lambda = lambda,
    s = sum(diff(y[, 1]) < 0) + 1,
    ctrl = list(tol_l = tol_l)),
    class = "opt")
  .check(ret)

  return(ret)
}

# Sanity check of opt object
# @param x Object of class \code{\link{opt}}.
# @param ... Additional arguments passed to \code{.check.opt}.
# @inherit .check description return
# @keywords internal
.check.opt <- function(x, ...) {
  with(x, {
    stopifnot(s == sum(diff(y[, 1]) < 0) + 1)
    stopifnot(is.numeric(y),
      is.matrix(y))
    if (is.logical(lambda_decrease)) {
      if (is.na(lambda_decrease)) stop("lambda_decrease in 'opt'-object must be TRUE or FALSE.")
    } else {
      stop("lambda_decrease in 'opt'-object must be TRUE or FALSE.")
    }
    stopifnot(length(nlambda) == 1,
      is.numeric(nlambda),
      is.finite(nlambda))
    if (!is.null(weights)) {
      stopifnot(is.matrix(weights),
        all(is.finite(weights)))
    }
    if (!is.null(lambda_min_ratio)) {
      stopifnot(length(lambda_min_ratio) == 1,
        is.numeric(lambda_min_ratio),
        all(lambda_min_ratio > 0),
        all(is.finite(lambda_min_ratio)))
    }
    with(ctrl, {
      stopifnot(is.numeric(tol_l),
        length(tol_l) == 1,
        tol_l > 0)
    })
    if (ncol(y) <= 1) {
      stop("y in 'opt'-object must have at least two columns (the first being time).")
    }
    if (any(!is.finite(y[, 1]))) {
      stop("First column of y in 'opt'-object (time points) have non finites or missing values.")
    }
    if (!is.null(weights)) {
      if (any(dim(weights) != (dim(y) - c(0, 1)))) {
        stop("Dimension of weights and y in 'opt'-object do not match (y should have one more column that weights).")
      }
      if (any(weights < 0 | !is.finite(weights))) {
        stop("Some weights in 'opt'-object are either negative or non-finite.")
      }
    }
    if (is.null(lambda)) {
      if (nlambda < 1) stop("nlambda is smaller than one.")
      if (!is.null(lambda_min_ratio)) {
        if (lambda_min_ratio >= 1 || lambda_min_ratio <= 0) stop("lambda_min_ratio in 'opt'-object is not in (0,1).")
      }
    } else {
      if (any(!is.numeric(lambda))) stop("Some lambda values in 'opt'-object are not numeric.")
      if (any(!is.finite(lambda))) stop("Some lambda values in 'opt'-object are non-finite.")
      if (any(lambda <= 0)) stop("Some lambda values in 'opt'-object are negative.")
    }

    # Check for number of obs
    if (any(diff(which(diff(y[, 1]) < 0)) == 1) || (diff(y[, 1]) < 0)[length(y[, 1]) - 1] || (diff(y[, 1]) < 0)[1]) {
      stop("Some context(s) of data in 'opt'-object only has one observation. A sharp decrement in the first column of y indicates a new context.")
    }
  })
}


### Compatibility check between ro and ode object
.o_op_check <- function(o, op) {
  if (methods::is(o, "ode")) {
    .check(o)
  } else {
    stop("'o' is not of class ode.")
  }
  if (methods::is(op, "opt")) {
    .check(op)
  } else {
    stop("'op' is not of class opt.")
  }
  d <- o$d
  ps <- lapply(o$rs, function(re) re$p)
  n <- nrow(op$y) - op$s
  s <- op$s

  lapply(o$rs, function(re) {
    if (!is.null(re$contexts)) {
      if (ncol(re$contexts) != s) stop("Number of contexts according to 'y' does not match that specified in 'reg' object.")

      stopifnot(is.matrix(re$contexts))
      if (ncol(re$contexts) != s) {
        stop(paste0("'contexts' in reg-object is non-null, but the number of columns (", ncol(re$contexts), ") does not match s (", s, "), where s is number of contexts according to 'y'."))
      }
    }
  })

  if (n <= 0) {
    stop("Time (first column of y) is decreasing, thus no degrees of freedom left for estimation (see details on er-function).")
  }
  if (ncol(op$y) - 1 != d) {
    stop(paste0("Number of observed coordinates (", ncol(op$y) - 1, ") (ncol(y) - 1 in opt-object) is not equal to state dimension, d = ", d, ", in ode-object."))
  }
}


# Here x0 is cube at lambda is based on max over them all
.set_defaults <- function(o, op) {
  .o_op_check(o, op)
  d <- o$d
  ps <- lapply(o$rs, function(re) re$p)
  n <- nrow(op$y) - op$s

  if (is.null(op$lambda_min_ratio)) {
    if (n * d < Reduce("+", ps)) {
      op$lambda_min_ratio <- 0.01
    } else {
      op$lambda_min_ratio <- 0.0001
    }
  }

  if (all(unlist(lapply(o$rs, function(re) re$reg_type == "none")))) {
    op$lambda  <- 1
    op$nlambda <- 1
  }

  return(op)
}

