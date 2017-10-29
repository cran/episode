###------------------###
###  Solver objects  ###
###------------------###

#' Create 'solver' object
#'
#' @description This function creates an object of class \code{solver}, which holds the basic information of numeric solver applied to the \code{ode}-systems.
#'
#' @param name Character string naming the ODE-solver. Must be one of: "rk23" (Runge-Kutta order 2/3), "bs23" (Bogacki-Shampine order 2/3), "dp45" (Dormand-Prince order 4/5) or "rkf45" (Runge-Kutta-Fehlberg order 4/5, default).
#' @param step_max Positive integer giving the maximal number of steps the solver may take between two consecutive time points.
#' @param tol Positive numeric tolerance level used for embedded pair solver.
#' @param h_init Positive numeric giving initial discretisation of time interval for solver.
#' @param ... Additional arguments passed to \code{\link{solver}}.
#' @seealso ode
#'
#' @return An object with S3 class "solver".
#'
#' @examples
#' # Use 'solver' object to specify numerical solver when creating 'ode' objects
#'
#' # Example: power law kinetics with Dormand-Prince order 4/5 solver
#' A <- matrix(
#' c(1, 1, 0, 0,
#'   0, 0, 1, 0,
#'   0, 0, 1, 0), ncol = 4, byrow = TRUE)
#' p <- plk(A, s = solver("dp45"))
#'
#' # Example: ... and with more steps
#' p <- plk(A, s = solver("dp45", step_max = 1e3))
#'
#' # Example: rational mass action kinetics with Runge-Kutta order 2/3 solver
#' B <- matrix(
#' c(0, 0, 1, 0,
#'   1, 1, 0, 0,
#'   1, 0, 0, 1), ncol = 4, byrow = TRUE)
#' rmak <- ratmak(A, B, s = solver("rk23"))
#'
#' @export
solver <- function(name = "rkf45", step_max = 100, tol = 1e-6, h_init = 1e-4, ...) {
  ret <- structure(list(name = name, ctrl = list(step_max = step_max, tol = tol, h_init = h_init)),
    class = "solver")
  .check(ret)
  return(ret)
}

# Sanity check of solver object
# @param x Object of class \code{\link{solver}}.
# @param ... Additional arguments passed to \code{.check.solver}.
# @inherit .check description return
# @keywords internal
.check.solver <- function(x, ...) {
  with(x, {
    if (!(name %in% .solver_names())) {
      stop(paste("'name' of solver not recognised. Must be one of:", paste0(.solver_names(), collapse = ", ")))
    }
    with(ctrl, {
      stopifnot(is.numeric(step_max),
        length(step_max) == 1,
        step_max >= 1)
      stopifnot(is.numeric(tol),
        length(tol) == 1,
        tol > 0)
      stopifnot(is.numeric(h_init),
        length(h_init) == 1,
        h_init > 0)
    })
  })
}
.solver_names <- function() {
  c('Runge-Kutta of order 2/3' = "rk23",
    'Bogacki-Shampine of order 2/3' = "bs23",
    'Dormand-Prince of order 4/5' = "dp45",
    'Runge-Kutta-Fehlberg of order 4/5' = "rkf45")
}



#' Print 'solver' object
#'
#' @description Prints objects of class \code{solver}.
#'
#' @param x Object of class \code{solver}.
#' @param ... Additional arguments passed to \code{\link{print}}.
#' @seealso solver
#'
#' @examples
#' # Dormand-Prince 4/5 scheme
#' solver("dp45")
#'
#' # Runge-Kutta-Fehlberg 4/5 scheme (default)
#' solver()
#'
#' @export
print.solver <- function(x, ...) {
  catname <- .solver_names()
  catname <- names(catname[catname == x$name])

  catctrl <- c(
    paste0("Maximum no. steps           ", x$ctrl$step_max, "\n"),
    paste0("Tolerance level             ", x$ctrl$tol, "\n"),
    paste0("Initial step discretisation ", x$ctrl$h_init))

  cat("Solver type:", catname, "\n",
    "Control parameters\n",
    paste0(rep("-", max(nchar(catctrl))), collapse = ""), "\n",
    catctrl)
}

