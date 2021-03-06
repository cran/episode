% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/solver.R
\name{solver}
\alias{solver}
\title{Create 'solver' object}
\usage{
solver(name = "rkf45", step_max = 100, tol = 1e-06, h_init = 1e-04, ...)
}
\arguments{
\item{name}{Character string naming the ODE-solver. Must be one of: "rk23" (Runge-Kutta order 2/3), "bs23" (Bogacki-Shampine order 2/3), "dp45" (Dormand-Prince order 4/5) or "rkf45" (Runge-Kutta-Fehlberg order 4/5, default).}

\item{step_max}{Positive integer giving the maximal number of steps the solver may take between two consecutive time points.}

\item{tol}{Positive numeric tolerance level used for embedded pair solver.}

\item{h_init}{Positive numeric giving initial discretisation of time interval for solver.}

\item{...}{Additional arguments passed to \code{\link{solver}}.}
}
\value{
An object with S3 class "solver".
}
\description{
This function creates an object of class \code{solver}, which holds the basic information of numeric solver applied to the \code{ode}-systems.
}
\examples{
# Use 'solver' object to specify numerical solver when creating 'ode' objects

# Example: power law kinetics with Dormand-Prince order 4/5 solver
A <- matrix(
c(1, 1, 0, 0,
  0, 0, 1, 0,
  0, 0, 1, 0), ncol = 4, byrow = TRUE)
p <- plk(A, s = solver("dp45"))

# Example: ... and with more steps
p <- plk(A, s = solver("dp45", step_max = 1e3))

# Example: rational mass action kinetics with Runge-Kutta order 2/3 solver
B <- matrix(
c(0, 0, 1, 0,
  1, 1, 0, 0,
  1, 0, 0, 1), ncol = 4, byrow = TRUE)
rmak <- ratmak(A, B, s = solver("rk23"))

}
\seealso{
ode
}
