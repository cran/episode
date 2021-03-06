% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rodeo.R
\name{rodeo}
\alias{rodeo}
\title{Regularised Ordinary Differential Equation Optimisation (RODEO) generic}
\usage{
rodeo(x, ...)
}
\arguments{
\item{x}{Object that specifies ode system.}

\item{...}{Additional arguments passed to \code{rodeo}.}
}
\value{
An object with S3 class "rodeo":
\item{o}{Original \code{ode}-object.}
\item{op}{Original \code{opt}-object with default values for \code{lambda_min_ratio}, \code{lambda} and (if needed) \code{a} inserted, if these were originally \code{NULL}.}
\item{params}{Parameter estimates, stored as list of sparse column format matrices, "dgCMatrix" (or a list of those if multiple initialisations). Rows represent coordinates and columns represent the \code{lambda} value.}
\item{x0s}{Initial state estimates stored in a matrix (or array). Rows represent coordinates, columns represent the \code{lambda} value and (if multiple initialisations) slices represent initialisations.}
\item{dfs}{A matrix (or array, if multiple initialisations) of degrees of freedom. Row represents a parameter (the first is always the initial state parameter), columns represent lambda, slices represent initialisation, if multiple are provided.}
\item{codes}{A matrix (or array) of convergence codes organised as \code{dfs}.
\describe{
\item{0:}{The convergence criteria is met (see details in \code{\link{opt}}). Current estimate is probably a local minimum.}
\item{1:}{Backtracking in the last iteration yields no numerical improvement, but no unsual behavior observed. Current estimate is probably a local minimum. However, if \code{exact_gradient} = FALSE in the \code{reg}-object in the \code{ode}-object, changing this may improve the code. Alternatively one can adjust backtracking via \code{backtrack_max} and \code{tau_min} in \code{reg} objects in \code{ode} object.}
\item{2:}{The optimisation procedure exceeded maximal number of steps (\code{step_max} in \code{reg} objects).}
\item{3:}{The last gradient was unsually large. Either the tolerances in \code{reg} objects are off or the ODE systems is very sensitive and runs over long time spans. In the latter case, initialisation(s) may have inappropriate zeros (change initialisation and/or make sure they start at smaller lambda value).}
\item{4:}{The numeric ODE solver exceeded maximal number of steps. Check if supplied initial states were out of bounds, if not increase \code{step_max} (or \code{tol}) in \code{reg}-objects in \code{ode}-object.}
}
}
\item{steps}{A matrix (or array) holding number of steps used in optimisation procedure. Organised as \code{dfs}.}
\item{losses}{A vector (or matrix) of unpenalised losses at optimum for each lambda value (stored row-wise if  multiple are provided).}
\item{penalties}{A matrix (or array) of penalties for each parameter, organised as \code{dfs}.}
\item{jerr}{A matrix (or array) of summary codes (for internal debugging), organised as \code{dfs}.}
}
\description{
Fit the parameters (and optionally initial states) for an Ordinary Differential Equation model with data sampled across contexts. Two implementations exist:
\describe{
\item{\code{\link{rodeo.ode}}}{The raw form of rodeo, based on \code{\link{ode}} (created via \code{\link{mak}}, \code{\link{plk}}, etc.) and \code{\link{opt}}-objects. If needed, everything can be user-specified through these two objects. However, parameter initialisations are always required and no action towards adjusting weights or scales of the parameters are taken.}
\item{\code{\link{rodeo.aim}}}{The more automatic form of rodeo, where automatic parameter initialisations and (optional) adjustments of scales and weights are available through \code{\link{aim}}.}
}
}
\details{
For details on the loss function, optimisation, etc. See documentation of \code{\link{opt}}.
}
\examples{
set.seed(123)
# Michaelis-Menten system with two 0-rate reactions
A <- matrix(c(1, 1, 0, 0,
 0, 0, 1, 0,
 0, 0, 1, 0,
 0, 0, 0, 1,
 1, 0, 0, 0), ncol = 4, byrow = TRUE)
B <- matrix(c(0, 0, 1, 0,
 1, 1, 0, 0,
 1, 0, 0, 1,
 0, 0, 1, 0,
 1, 0, 1, 0), ncol = 4, byrow = TRUE)
k <- c(1.1, 2.5, 0.25, 0, 0); x0 <- c(E = 2, S = 2, ES = 8.5, P = 2.5)
Time <- seq(0, 10, by = 1)

# Simulate data, in second context the catalytic rate has been doubled
contexts <- cbind(1, c(1, 1, 2, 1, 1))
m <- mak(A, B, r = reg(contexts = contexts))
y <- numsolve(m, c(Time, Time), cbind(x0, x0 + c(2, -2, 0, 0)), contexts * k)
y[, -1] <- y[, -1] + matrix(rnorm(prod(dim(y[, -1])), sd = .5), nrow = nrow(y))

# Example: fit data using rodeo on mak-object
op <- opt(y, nlambda = 10)
fit <- rodeo(m, op, x0 = NULL, params = NULL)
fit$params$rate

# Example: fit dat using rodeo on aim-object
a <- aim(m, op)
a$params$rate
fit <- rodeo(a)
fit$params$rate

}
\seealso{
rodeo.aim, rodeo.ode
}
