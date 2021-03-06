% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ode.R
\name{print.ratmak}
\alias{print.ratmak}
\title{Print 'ratmak' object}
\usage{
\method{print}{ratmak}(x, ...)
}
\arguments{
\item{x}{Object of class \code{ratmak}.}

\item{...}{Additional arguments passed to \code{\link{print}}.}
}
\description{
Prints objects of class \code{ratmak}, and presents the powers.
}
\examples{
# Rational mass action kinetics
A <- matrix(
c(1, 0, 0, 0,
  0, 1, 2, 0,
  1, 0, 0, 1), ncol = 4, byrow = TRUE)
x0 <- c(X = 1, Y = 4, Z = 0.1, W = 0.1)
time <- seq(0, 1, by = .1)

rmak <- ratmak(A, diag(4))
rmak
}
\seealso{
ratmak
}
