###-------------------------###
###  Regularisation object  ###
###-------------------------###


#' Create 'reg' (regularisation) object
#'
#' This function creates an object of class \code{reg}, which holds regularisation type, tuning parameter scales, penalty factors and control list for optimisation. This is basically the control panel for the optimisation in \code{\link{rodeo}} and \code{\link{aim}}.
#'
#'
#' @param reg_type Character string determining the regularisation. Must be one of: "l1" (default), "l2", "elnet", "scad", "mcp" or "none". See details.
#' @param a Numeric value of tuning parameter in \code{elnet} (must be between 0 and 1, default is 0.5), \code{scad} (must be larger than 2, default = 3.7) or \code{mcp} (must be larger than 1, default = 3).
#' @param lower Either numeric vector of length 1 or p, of lower limit(s) for parameters, must be smaller or equal to \code{upper}.
#' @param upper Either numeric vector of length 1 or p, of upper limit(s) for parameters, must be larger or equal to \code{lower}.
#' @param lambda_factor Non-negative numeric value which will be multiplied with the tuning parameter in \code{\link{opt}}.
#' @param exclude A vector indices to exclude from model (this is how to specify an infinite penalty factor). Default is none.
#' @param penalty_factor A non-negative vector (p) of individual penalty weights. Defaults to 1 for each parameter. Will always be rescaled so its mean is 1.
#' @param contexts A non-negative matrix (pxs) specifying design on context-level, see details. Defaults to a matrix of ones.
#' @param scales A positive vector (p) of scales, see details. Defaults to a vector of ones. Will be rescaled to mean to 1.
#' @param fixed Logical indicating if this parameter is fixed or should be optimised.
#' @param screen Logical indicating if a faster optimisation relying on screening should be adapted. If \code{NULL}, it is set to \code{TRUE} iff p > 20.
#' @param exact_gradient Logical indicating if exact gradient should be used. If \code{NULL}, it is set to \code{TRUE} iff p > 20.
#' @param step_max Positive integer giving the maximal number of steps in optimisation procedure, per lambda.
#' @param step_screen Positive integer giving the number of steps between screenings (defaults to \code{ceiling(step_max / 50)}).
#' @param step_cycle Positive integer giving the number of steps per optimisation cycle, defaults to \code{step_screen}.
#' @param backtrack_max Positive integer giving the maximal number of backtracking steps taken in each optimisation step.
#' @param tol Positive numeric tolerance level used for parameter space.
#' @param tau_init Positive initial step length for backtracking.
#' @param tau_min Positive numeric giving minimal value of step length in backtracking.
#' @param tau_scale Scaling parameter of step length, must be in (0,1).
#' @seealso aim, rodeo, rodeo.aim, rodeo.ode
#' @details
#' \subsection{Data format}{
#' The data (\code{y} in \code{\link{opt}}) is assumed generated from s different contexts. A new context begins whenever the time (the first column of \code{y} in \code{\link{opt}}) decreases. Hence s - 1 is the number of decreases in the time points.
#'
#' Each context has its own initial condition and parameter vector specified through \code{contexts} in \code{\link{reg}}. More precisely, the effective parameter in context l is the element-wise product of a baseline parameter (scaled by \code{scales} in \code{\link{reg}}) and the l'th column of \code{contexts}. This enables the user to pre-specify different scales for each parameter, as well as different scales for the same parameter across contexts.
#'
#' The default setup sets \code{contexts} to a matrix of ones. The following are examples of case-specific modifications for the \code{\link{mak}} ODE class:
#' If reaction \code{j} does not take place in context \code{l} then set \eqn{contexts_{j,l} = 0}.
#' If reaction \code{j} has a \eqn{50\%} increase in rate in context \code{l} then set \eqn{contexts_{j,l} = 1.5}.
#' If reaction \code{j} has independent rates in contexts \code{l} and \code{l}', then write that reaction twice in \code{\link{mak}}-object (with \code{j}' denoting its second appearance) and set \eqn{contexts_{j,l'} = 0} and \eqn{contexts_{j',l} = 0}.
#'}
#'\subsection{Loss function}{
#'The loss function optimised in \code{\link{rodeo}} is:
#'\deqn{RSS/(2 * (n - s)) + \lambda*\sum_{parameter \ argument}\lambda_{factor}*\sum_{j=1}^p v_jpen(param_j)}
#' where \eqn{\lambda} belongs to the \code{lambda}-sequence in \code{\link{opt}}-object and v is \code{penalty_factor} in \code{\link{reg}}. Moreover, the residual sum of squares, RSS, is given as:
#'\deqn{RSS = \sum^{n}_{i=1}||(y(t_i) - x(t_i, {x_0}_l, context_l * param))*\sqrt{w(t_i)}||_2^2}
#' where \code{param} has been (internally) scaled with \code{scales}, and \eqn{w(t_i)} and \eqn{y(t_i)} refers to the i'th row of \code{weights} and \code{y} in \code{\link{opt}} (\code{y} with first column removed), respectively. The solution to the ODE system is the x()-function. The subscript 'l' refers to the context, i.e., the columns of \code{contexts}  and \code{x0} in \code{rodeo}-functions (\code{x0} is the initial state of the system at the first time point after a decrease in the time-vector). All products are Hadamard products.
#'}
#'\subsection{Regularisation}{
#'The type of regularisation is chosen via \code{reg_type} in \code{\link{reg}}:
#'\describe{
#'\item{\code{l1}:}{Least Absolute Shrinkage and Selection Operator (Lasso) penalty. The penalty is the absolute value: \eqn{pen(param_j)=|param_j|}}
#'\item{\code{l2}:}{Ridge penalty. The penalty is the squaring: \eqn{pen(param_j)=param_j^2/2}}
#'\item{\code{elnet}:}{Elastic Net. A convex combination of \code{l1} and \code{l2} penalties: \eqn{pen(param_j)=(1-a)param_j^2/2+a|param_j|}, 0<=a<=1. Note if a = 0 or a = 1, then the penalty is automatically reduced to "l2" and "l1" respectively.}
#'\item{\code{scad}:}{Smoothly Clipped Absolute Deviation penalty: \deqn{pen(param_j)=\int^{param_j}_{0}{min( max((a\lambda-|param|)/(\lambda(a-1)), 0), 1) dparam}, a > 2.}}
#'\item{\code{mcp}:}{Maximum Concave Penalty penalty: \deqn{pen(param_j)=\int^{param_j}_{0}{max(1 - |param|/(\lambda a), 0) dparam}, a > 1.}}
#'\item{\code{none}:}{No penalty. Not recommended for large systems. This overwrites both user-supplied and automatically generated lambda-sequences.}
#'}
#'}
#'\subsection{Optimisation}{
#'  A proximal-gradient type of optimisation is employed. The step length is denoted \eqn{\tau}.
#'
#'  The flag \code{screen} (in \code{\link{reg}}) is passed onto the optimisation, it simply tells the optimisation to focus entirely on optimising a subset of the parameters, which where selected due to large gradients. At most \code{step_screen} (in \code{\link{reg}}) steps are taken before a re-evaluation of the screened subset takes place.
#'
#' The convergence criteria is \eqn{\Delta \eta < 10^3 * \max(|(\Delta param \neq 0)| + |(\Delta x0 \neq 0)|, 1)} AND \eqn{\Delta loss / \Delta \eta < tol_l}, where
#'  \deqn{\Delta \eta = ||\Delta param||/tol_{param} + ||\Delta x0|| / tol_{x_0}}
#'
#'}
#'\subsection{Tuning parameter}{
#' The lambda sequence can either be fully customised through \code{lambda} in \code{\link{opt}} or automatically generated. In the former case, a monotonic \code{lambda}-sequence is highly recommended. Throughout all optimisations, each optimisation step re-uses the old optimum, when sweeping through the \code{lambda}-sequence.
#'
#' If \code{lambda} is \code{NULL}, an automatically generated sequence is used. A maximal value of lambda (the smallest at which 0 is a optimum in the rate parameter) is calculated and log-equidistant sequence of length \code{nlambda} is generated, with the ratio between the smallest and largest lambda value being \code{lambda.min.ratio}. Note: when the \code{opt}-object is passed to \code{\link{rodeo.ode}}, one may want to initialise the optimisation at a non-zero parameter and run the optimisation on the reversed lambda sequence. This is indicated by setting \code{decrease.lambda} = FALSE. If, however, the \code{opt}-object is passed to \code{\link{aim}}, \code{\link{glmnet}} ultimately decides on the lambda sequence, and may cut it short.
#'}
#'\subsection{Gradient evaluations}{
#' Two gradient evaluations are implemented: exact and approximate. The computational time of exact gradient is generally longer than that of approximate gradient, but its relative magnitude depends on the \code{solver} in the \code{ode}-object passed to \code{\link{rodeo}}.
#' }
#' @return An object with S3 class "reg".
#'
#' @examples
#' # Use 'reg' object to specify regularisation when creating 'ode' objects
#'
#' # Example: power law kinetics with SCAD penalty on rate parameter
#' A <- matrix(
#' c(1, 1, 0, 0,
#'   0, 0, 1, 0,
#'   0, 0, 1, 0), ncol = 4, byrow = TRUE)
#' p <- plk(A, r = reg("scad"))
#'
#' # Example: power law kinetics as above, with lower bound of -1
#' p <- plk(A, r = reg("scad", lower = -1))
#'
#' # Example: rational mass action kinetics with
#' # MCP penalty on first parameter argument and l2 on second
#' B <- matrix(
#' c(0, 0, 1, 0,
#'   1, 1, 0, 0,
#'   1, 0, 0, 1), ncol = 4, byrow = TRUE)
#' rmak <- ratmak(A, B, r1 = reg("mcp"), r2 = reg("l2"))
#'
#' # Example: mass action kinetics with
#' # no penalty on rate parameter and l2 on initial state
#' m <- mak(A, B, r = reg("none"), rx0 = reg("l2"))
#'
#' @export
reg <- function(reg_type = "l1", a = NULL, lower = -Inf, upper = Inf, lambda_factor = 1, exclude = NULL, penalty_factor = NULL, contexts = NULL, scales = NULL, fixed = FALSE, screen = NULL, exact_gradient = NULL, step_max = 200, step_screen = ceiling(step_max / 10), step_cycle = step_screen, backtrack_max = 50, tol = 1e-5, tau_init = .1, tau_min = 1e-10, tau_scale = .5) {
  ret <- structure(list(reg_type = reg_type, p = NULL, a = a, lower = lower, upper = upper,
    lambda_factor = lambda_factor, exclude = unique(exclude),
    penalty_factor = penalty_factor, contexts = contexts, scales = scales, fixed = isTRUE(fixed),
    screen = screen, exact_gradient = exact_gradient,
    ctrl = list(step_max = step_max, step_cycle = step_cycle, step_screen = step_screen,
      backtrack_max = backtrack_max,
      tol = tol,
      tau_init = tau_init, tau_min = tau_min, tau_scale = tau_scale)),
    class = "reg")
  .check(ret)

  # penalty_factor
  if (!is.null(penalty_factor)) {
    m <- mean(penalty_factor)
    if (m <= 0 || !is.finite(m)) {
      stop("'penalty_factor' is not null, but its mean is non-positive or infinite.")
    }
    ret$penalty_factor <- penalty_factor / m
  }

  # scales
  if (!is.null(scales)) {
    m <- mean(scales)
    if (m <= 0 || !is.finite(m)) {
      stop("'scales' is not null, but its mean is non-positive or infinite.")
    }
    ret$scales <- scales / m
  }

  # a
  if (reg_type == "l1") {
    ret$a <- 1
  } else if (reg_type == "l2") {
    ret$a <- 0
  } else if (reg_type == "elnet") {
    if (is.null(a)) {
      ret$a <- 0.5
    } else {
      if (length(a) != 1) stop("a is not null and does not have length 1.")
      if (a > 1 || a < 0) {
        stop("Tuning parameter, a, in elnet must be in [0, 1].")
      } else if (a == 1) {
        warning("Supplied 'a' was 1, reg_type changed from elnet to l1.")
        ret$reg_type <- "l1"
        ret['a'] <- list(NULL)
      } else if (a == 0) {
        warning("Supplied 'a' was 0, reg_type changed from elnet to l2.")
        ret$reg_type <- "l2"
        ret['a'] <- list(NULL)
      }
    }
  } else if (reg_type == "scad") {
    if (is.null(a)) {
      ret$a <- 3.7
    } else {
      if (length(a) != 1) stop("a is not null and does not have length 1.")
      if (a <= 2) {
        stop("Tuning parameter, a, in SCAD must be larger than 2.")
      }
    }
  } else if (reg_type == "mcp") {
    if (is.null(a)) {
      ret$a <- 3
    } else {
      if (length(a) != 1) stop("a is not null and does not have length 1.")
      if (a <= 1) {
        stop("Tuning parameter, a, in MCP must be larger than 1.")
      }
    }
  } else if (reg_type == "none") {
    ret$a <- 0
  } else {
    stop("error in reg -> reg_type not recognised.")
  }

  return(ret)
}

# Default reg for x0
.reg_x0 <- function(reg_type = "none", a = NULL, lower = -Inf, upper = Inf, lambda_factor = 1, exclude = NULL, penalty_factor = NULL, contexts = NULL, scales = NULL, fixed = TRUE, exact_gradient = NULL, step_max = 200, step_screen = ceiling(step_max / 10), step_cycle = step_screen, backtrack_max = 50, tol = 1e-5, tau_init = .1, tau_min = 1e-10, tau_scale = .5) {
  return(reg(reg_type = reg_type, a = a, lower = lower, upper = upper, lambda_factor = lambda_factor, exclude = exclude, penalty_factor = penalty_factor, contexts = contexts, scales = scales, fixed = fixed, screen = FALSE, exact_gradient = exact_gradient, step_max = step_max, step_screen = step_screen, step_cycle = step_cycle, backtrack_max = backtrack_max, tol = tol, tau_init = tau_init, tau_min = tau_min, tau_scale = tau_scale))
}


# Sanity check of reg object
# @param x Object of class \code{\link{reg}}.
# @param ... Additional arguments passed to \code{.check.reg}.
# @inherit .check description return
# @keywords internal
.check.reg <- function(x, ...) {
  with(x, {
    if (!(reg_type %in% .reg_names())) {
      stop(paste("reg_type not recognised. Only accepts:", paste0(.solver_names(), collapse = ", ")))
    }
    if (!is.null(screen)) {
      if (is.na(screen) | !is.logical(screen)) stop("screen must be TRUE or FALSE.")
    }
    if (!is.null(fixed)) {
      if (is.na(fixed) | !is.logical(fixed)) stop("screen must be TRUE or FALSE.")
    }
    if (!is.null(exact_gradient)) {
      if (is.na(exact_gradient) | !is.logical(exact_gradient)) stop("exact_gradient must be TRUE or FALSE.")
    }
    if (!is.null(a)) {
      stopifnot(is.numeric(a))
    }

    # Bounds param
    if (!is.null(lower)) {
      if (!all(is.finite(lower) | is.infinite(lower))) stop("Elements of 'lower' are neither finite nor infinite.")
    }
    if (!is.null(upper)) {
      if (!all(is.finite(upper) | is.infinite(upper))) stop("Elements of 'upper' are neither finite nor infinite.")
    }
    if (!is.null(upper) & !is.null(lower)) {
      if (any(lower > upper)) stop("Some entries of 'lower' exceeds those of 'upper'.")
    }
    if (!is.null(penalty_factor)) {
      stopifnot(is.numeric(penalty_factor),
        all(penalty_factor >= 0),
        all(is.finite(penalty_factor)))
    }
    stopifnot(length(lambda_factor) == 1, lambda_factor > 0)
    if (!is.null(contexts)) {
      stopifnot(is.numeric(contexts),
        is.matrix(contexts),
        all(contexts >= 0),
        all(is.finite(contexts)))
    }
    if (!is.null(scales)) {
      stopifnot(is.numeric(scales),
        all(scales > 0),
        all(is.finite(scales)))
    }
    with(ctrl, {
      stopifnot(is.numeric(step_max),
        length(step_max) == 1,
        step_max >= 1)
      stopifnot(is.numeric(step_cycle),
        length(step_cycle) == 1,
        step_cycle >= 1)
      stopifnot(is.numeric(step_screen),
        length(step_screen) == 1,
        step_screen >= 1)
      stopifnot(is.numeric(backtrack_max),
        length(backtrack_max) == 1,
        backtrack_max >= 1)
      stopifnot(is.numeric(tol),
        length(tol) == 1,
        tol > 0)
      stopifnot(is.numeric(tau_init),
        length(tau_init) == 1,
        tau_init > 0)
      stopifnot(is.numeric(tau_min),
        length(tau_min) == 1,
        tau_min > 0)
      if (tau_init <= tau_min) stop("tau_init must be larger than tau_min.")
    })

    # scale check
    if (ctrl$tau_scale >= 1 || ctrl$tau_scale <= 0) stop("tau_scale not in (0,1).")

    # p check
    if (!is.null(p)) {
      if (!(length(lower) == 0 | length(lower) == 1 | length(lower) == p)) stop("Length of lower in 'reg' is neither 0, 1 or p.")
      if (!(length(upper) == 0 | length(upper) == 1 | length(upper) == p)) stop("Length of upper in 'reg' is neither 0, 1 or p.")
      if (!is.null(exclude)) {
        if (!all(exclude %in% seq_len(p))) stop("Not all exclude in 'reg' are in 1:p.")
      }
      if (!is.null(contexts)) {
        if (nrow(contexts) != p) stop("The number of rows in context in 'reg' is not p.")
      }
      if (!is.null(scales)) {
        if (length(scales) != p) stop("The length of scales in 'reg' is not p.")
      }
    }

    # contexts and scales check
    if (!is.null(contexts)) {
      if (any(contexts < 0)) stop("Some entries in 'contexts' in reg-object are negative.")
    }
    if (!is.null(scales)) {
      if (any(scales <= 0)) stop("Some entries in 'scales' in reg-object are non-positive.")
    }
  })
}


# Exclude parameters in reg object
#
# @description Adapts the \code{\link{reg}} object to the 'exclude' argument in \code{\link{reg}}-objects. Probably irrelevant for most users.
#
# @param x Object of class \code{\link{reg}} to adapt.
# @param ... Additional arguments passed to \code{.exclude.reg}.
# @return An object of class \code{reg}, but with entries in \code{penalty_factor} indexed by \code{.exclude} removed.
# @keywords internal
.exclude_reg <- function(x) {
  exc <- x$exclude
  if (!is.null(exc) & length(exc) > 0) {
    if (!is.null(x$penalty_factor)) x$penalty_factor <- x$penalty_factor[-exc]
    if (length(x$lower) > 1) x$lower <- x$lower[-exc]
    if (length(x$upper) > 1) x$upper <- x$upper[-exc]
    if (!is.null(x$p)) x$p <- x$p - length(exc)
  }
  return(x)
}


.get_sc <- function(re, op) {
  with(re, {
    if (is.null(scales)) {
      if (is.null(contexts)) {
        NULL
      } else {
        contexts
      }
    } else {
      if (is.null(contexts)) {
        matrix(scales, nrow = length(scales), ncol = op$s)
      } else {
        scales * contexts
      }
    }
  })
}





#' Print 'reg' object
#'
#' @description Prints objects of class \code{reg}.
#'
#' @param x Object of class \code{reg}.
#' @param ... Additional arguments passed to \code{\link{print}}.
#' @seealso reg
#'
#' @examples
#' # Example: default regularisation
#' reg()
#'
#' @export
print.reg <- function(x, ...) {
  catname <- .reg_names()
  catname <- names(catname[catname == x$reg_type])
  if (x$reg_type %in% c("scad", "mcp")) catname <- paste0(catname, " (a = ", x$a, ")")

  catctrl <- c(
    paste0("Maximal no. of steps      ", x$ctrl$step_max, "\n "),
    paste0("Steps per cycle           ", x$ctrl$step_cycle, "\n "),
    ifelse(x$screen, paste0("Steps between screenings  ", x$ctrl$step_screen, "\n"), ""),
    paste0("Maximal backtracking      ", x$ctrl$backtrack_max, "\n "),
    paste0("Tolerance level           ", x$ctrl$tol, "\n "),
    paste0("Initial step length       ", x$ctrl$tau_init, "\n "),
    paste0("Minimal step length       ", x$ctrl$tau_min, "\n "),
    paste0("Step length factor        ", x$ctrl$tau_scale, "\n "))
  catctrl <- paste0("\n Control parameters\n ",
    paste0(rep("-", max(nchar(catctrl))), collapse = ""), "\n ",
    paste0(catctrl, collapse = ""), collapse = "")

  cat("Penalty:", catname, ifelse(x$fixed, "", catctrl))
}

.reg_names <- function() {
  c('l1 norm' = "l1",
    'l2 norm' = "l2",
    'Elastic Net' = "elnet",
    'Smoothly Clipped Absolute Deviation' = "scad",
    'Maximum Concave Penalty' = "mcp",
    'None' = "none")
}


