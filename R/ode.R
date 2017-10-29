###---------------###
###  ODE objects  ###
###---------------###


#' Abstract 'ode' object
#'
#' @description The abstract class \code{ode}, which determines the ordinary differential equation system. Some systems have multiple parameter arguments (so far only \code{\link{ratmak}}). The bare \code{ode} class is abstract, hence in practice you need to use one of
#'
#'\describe{
#' \item{\code{\link{mak}}}{Mass Action Kinetics reaction network system.}
#' \item{\code{\link{plk}}}{Power Law Kinetics system.}
#' \item{\code{\link{rlk}}}{Rational Law Kinetics system.}
#' \item{\code{\link{ratmak}}}{Rational Mass Action Kinetics.}
#' }
#'
#' to create an actual \code{ode} object, which you can use.
#'
#' @param ... Arguments passed to \code{\link{ode}}.
#' @seealso mak, plk, rlk, ratmak
#'
#' @return An object with S3 class "ode".
#'
#' @export
ode <- function(...) {
  return(structure(list(...), class = "ode"))
}
# things that must be implemented in new ode's: d, s (solver), rs (list of reg)


#' Print 'ode' object
#'
#' @description Prints objects of class \code{ode}. Since this class is abstract, it moves on the print method for the derived class.
#'
#' @param x Object of class \code{ode}.
#' @param ... Additional arguments passed to \code{\link{print}}.
#' @seealso ode
#'
#' @export
print.ode <- function(x, ...) {
  NextMethod()
}

# Exclude parameters in ode object
#
# @description Adapts the ode object to the 'exclude' argument in \code{\link{reg}}-objects. Probably irrelevant for most users.
#
# @param x Object of class \code{\link{ode}} to adapt.
# @param ... Additional arguments passed to \code{.exclude.ode}.
# @keywords internal
.exclude <- function (x, ...) {
  UseMethod(".exclude", x)
}
.exclude.ode <- function(x, ...) {
  NextMethod()
}


# Sanity check of ode object
# @param x Object of class \code{\link{ode}}.
# @param ... Additional arguments passed to \code{.check.ode}.
# @inherit .check description return
# @keywords internal
.check.ode <- function(x, ...) {
  with(x, {
    # Solver
    if (!methods::is(s, "solver")) {
      stop("'s' is not a 'solver'.")
    } else {
      .check(s)
    }

    # Reg
    lapply(rs, function(re) {
      if (!methods::is(re, "reg")) {
        stop("Elements in 'rs' are not 'reg'.")
      } else {
        .check(re)
      }
    })
  })

  NextMethod()

  # If one of limits is not singleton or empty, both must match if they are not null
  x$rs <- lapply(x$rs, function(re) {
    if (length(re$lower) > 1 | length(re$upper) > 1) {
      if (length(re$lower) == 1) re$lower <- rep(re$lower, length(re$upper))
      if (length(re$upper) == 1) re$upper <- rep(re$upper, length(re$lower))
    }

    # screen and exact_gradient
    if (is.null(re$screen) & !is.null(re$p)) {
      re$screen <- (re$p > 20)
    }
    if (is.null(re$exact_gradient) & !is.null(re$p)) {
      re$exact_gradient <- (re$p > 20)
    }

    re
  })

  return(x)
}



.get_ode_name <- function(o) {
  ode_names <- c("mak", "plk", "ratmak", "rlk")
  matches   <- sapply(ode_names, methods::is, object = o)
  if (all(!matches)) stop(paste0("ODE not recognised, must be one of: ", paste(ode_names, collapse = ", ")))
  if (sum(matches) > 1) stop("Multiple ODE classes recognised.")

  return(ode_names[matches][1])
}


#' @title Numerical solver for Ordinary Differential Equation (ODE) systems.
#'
#' @description Calculates numerical solution to the ODE system specified in the \code{ode} object and returns the solution path (and optionally its derivatives).
#'
#' @param o An object of class \code{\link{ode}} (created via \code{\link{mak}} or \code{\link{plk}}).
#' @param time A numeric vector holding desired time-points at which to evaluate the solution path. By convention of this package, whenever \code{time} decreases the parameters and initial conditions are reset, see details.
#' @param x0 A non-negative numeric vector or matrix holding the initial conditons for the ODE system.
#' @param param A non-negative numeric vector or matrix holding the parameter values for the ODE system (or a list of these if the ode system has multiple arguments, like \code{\link{ratmak}}).
#' @param approx_sensitivity Logical or NULL. If NULL (default) the sensitivity equations are not solved and no derivatives are returned. Otherwise if logical, it will solve the sensitivity equations. If TRUE an approximate solution is returned, else a more exact solution (based on the numerical solver) is returned.
#' @param ... Additional arguments passed to numeric solver.
#' @seealso ode, mak, plk, rlk, ratmak
#' @details As mentioned above, whenever \code{time} decreases the system restarts at a new initial condition and with a new parameter vector, called a context/environment. The number of contexts is s, i.e., the number of decreases + 1. To support these s contexts, the new initial conditions and parameter vectors can be supplied through \code{param} and \code{x0}. For both of these, either a matrix form holding the new initial conditions or parameters in the columns or a concatinated vector will work. In either case \code{param} and \code{x0} are coarced to matrices with p (= number of parameters) and d (= dimension of state space) rows respectively (with p and d extracted from \code{o}). Hence a check of whether \code{param} and \code{x0} have length divisible by p and d respectively is conducted. Both parameter vectors and initial conditions will be recycled if s exceeds the number of these divisors. Therefore, if the same parameter vector and/or initial condition is desired across resets, supply only that vector. A warning will be thrown if p or d is neither a multiple nor a sub-multiple of s.
#' @examples
#' # Example: Michaelis-Menten system
#' A <- matrix(
#' c(1, 1, 0, 0,
#'   0, 0, 1, 0,
#'   0, 0, 1, 0), ncol = 4, byrow = TRUE)
#' B <- matrix(
#' c(0, 0, 1, 0,
#'   1, 1, 0, 0,
#'   1, 0, 0, 1), ncol = 4, byrow = TRUE)
#' k <- c(1, 2, 0.5)
#' x0 <- c(E = 1, S = 4, ES = 0, P = 0)
#' Time <- seq(0, .5, by = .1)
#'
#' m <- mak(A, B)
#'
#' # Solution for one context
#' numsolve(m, Time, x0, k)
#'
#' # Solution for two contexts (the latter with faster rate)
#' numsolve(m, c(Time, Time), x0, cbind(k, k * 1.5))
#'
#' # Solution for two contexts (the latter with different initial condition)
#' numsolve(m, c(Time, Time), cbind(x0, x0 + 1.5), k)
#'
#' # As above, but with sensitivity equations are solved (using approximate solution)
#' numsolve(m, c(Time, Time), cbind(x0, x0 + 1.5), k, TRUE)
#'
#'
#' # Example: Power law kinetics
#' A <- matrix(c(1, 0, 1,
#'               1, 1, 0), byrow = TRUE, nrow = 2)
#' p <- plk(A)
#' x0 <- c(10, 4, 1)
#' theta <- matrix(c(0, -0.25,
#'                   0.75, 0,
#'                   0, -0.1), byrow = TRUE, nrow = 3)
#' numsolve(p, Time, x0, theta)
#'
#' @return If \code{approx_sensitivity} is \code{NULL} a matrix with the first column holding the \code{time} vector and all consecutive columns holding the states. A 0/1 convergence code of the solver is given as an attribute.
#'
#' If \code{approx_sensitivity} is logical a list is returned with
#' \item{trajectory}{A matrix containing the trajectory of the system (as with \code{approx_sensitivity = NULL}).}
#' \item{sensitivity_param}{A list with arrays of sensitivity solutions of parameters (row = time, colum = state, slice = parameter coordinate).}
#' \item{sensitivity_x0}{An array with sensitivity solutions of initial state (row = time, colum = state, slice = initial state coordinate).}
#'
#' For a convergence code of 1, try decreasing \code{tol} or increasing \code{step_max} in \code{\link{solver}} object in \code{ode} object. Note that the sensitivity equations may have non-finites, even if the convergence code is 0. This is typically due to zero initial states at which the state derivative of the field may be undefined.
#'
#' @export
numsolve <- function(o, time, x0, param, approx_sensitivity = NULL, ...) {
  if (!methods::is(o, "ode")) stop("'o' is not an 'ode'-object.")

  ode_name <- .get_ode_name(o)

  if (!methods::is(param, "list")) {
    param <- list(param)
  }
  if (!all(unlist(lapply(param, function(para) is.vector(para) || is.matrix(para))))) {
    stop("Not all parameters are vectors or matrices.")
  }
  if (!all(unlist(lapply(param, function(para) is.finite(para) || is.numeric(para))))) {
    stop("Not all parameters are finite or numeric.")
  }
  param <- lapply(param, as.vector)


  stopifnot(is.numeric(time),
    all(is.finite(time)),
    is.numeric(x0),
    all(is.finite(x0)),
    is.matrix(x0) || is.vector(x0))
  sensitivity <- !is.null(approx_sensitivity)
  approx <- sensitivity & isTRUE(approx_sensitivity)

  s     <- sum(diff(time) < 0) + 1
  o$rs[[1]]$p <- o$d * s

  if (length(param) != length(o$rs) - 1) {
    stop("Number of parameter arguments in 'param' is not the same as number of parameter arguments in 'o' (length of o$rs).")
  } else {
    Map(function(p1, p2) {
      if (p2 %% p1 != 0) stop(paste0("Length of param (", p1, ") is multiple of number of parameters (", p2, ") according to reg-objects.")) }, lapply(o$rs[-1], function(re) re$p), lapply(param, length))
  }

  if (length(x0) %% o$d != 0) stop(paste0("Length of x0 (", length(x0), ") is not a multiple of number of states (", o$d, ")."))
  param <- Map(function(para, re) matrix(para, nrow = re$p, ncol = s),
    para = param, re = o$rs[-1])
  param <- lapply(param, function(para) as.list(data.frame(para)))
  x0_   <- x0
  x0    <- matrix(x0, nrow = o$d, ncol = s)

  if (any(x0 < o$rs[[1]]$lower)) stop("Some coordinates of 'x0' are smaller than 'lower' limit of state.")
  if (any(x0 > o$rs[[1]]$upper)) stop("Some coordinates of 'x0' are larger than 'upper' limit of state.")
  Map(function(part, re) {
    lapply(part, function(para) {
      if (any(para < re$lower)) stop("Some coordinates of 'param' are smaller than 'lower' in 'o'.")
      if (any(para > re$upper)) stop("Some coordinates of 'param' are larger than 'upper' in 'o'.")
    })
  }, part = param, re = o$rs[-1])

  param <- .switch_hierarchy(param)

  ret <- ode_solve(ode_struct = o, x0 = x0, param_ = param, time = time,
    sensitivity = sensitivity, approx = approx)

  conv_code   <- ret$conv_code
  trajectory  <- cbind(Time = time, ret$trajectory)
  if (is.matrix(x0_)) {
    if (!is.null(rownames(x0_))){
      colnames(trajectory) <- c("Time", rownames(x0_)[seq_len(o$d)])
    } else {
      colnames(trajectory) <- c("Time", paste0("x", seq_len(o$d)))
    }
  } else {
    if (!is.null(names(x0_))){
      colnames(trajectory) <- c("Time", names(x0_)[seq_len(o$d)])
    } else {
      colnames(trajectory) <- c("Time", paste0("x", seq_len(o$d)))
    }
  }
  attr(trajectory, "conv_code") <- as.numeric(any(conv_code != 0))

  if (sensitivity) {
    return(list(trajectory = trajectory, sensitivity_param = ret$sens_dparam, sensitivity_x0 = ret$sens_dx0))
  } else {
    return(trajectory)
  }
}



#' @title Field of Ordinary Differential Equation (ODE) systems.
#'
#' @description Evaluates the vector field of the ODE system specified in the \code{ode} object (and optionally its differentials with respect to the state and parameter vectors).
#'
#' @param o An object of class \code{\link{ode}} (created via \code{\link{mak}} or \code{\link{plk}}).
#' @param x A non-negative numeric vector holding the state of the ODE system.
#' @param param A non-negative numeric vector of parameter values for the ODE system (or list of these if multiple arguments are required).
#' @param differentials Logical indicating if the differentials with respect to the state and parameter vector ought to be returned as well.
#' @param ... Additional arguments passed to field.
#' @seealso ode, mak, plk, rlk, ratmak
#' @examples
#' # Example: Michaelis-Menten system
#' A <- matrix(
#' c(1, 1, 0, 0,
#'   0, 0, 1, 0,
#'   0, 0, 1, 0), ncol = 4, byrow = TRUE)
#' B <- matrix(
#' c(0, 0, 1, 0,
#'   1, 1, 0, 0,
#'   1, 0, 0, 1), ncol = 4, byrow = TRUE)
#' k <- c(1, 2, 0.5)
#' x <- c(E = 1, S = 4, ES = 0, P = 0)
#'
#' m <- mak(A, B)
#'
#' # Vector field
#' field(m, x, k)
#'
#' # ... with differentials
#' field(m, x, k, TRUE)
#'
#'
#' # Example: Power law kinetics
#' A <- matrix(c(1, 0, 1,
#'               1, 1, 0), byrow = TRUE, nrow = 2)
#' p <- plk(A)
#' x <- c(10, 4, 1)
#' theta <- matrix(c(0, -0.25,
#'                   0.75, 0,
#'                   0, -0.1), byrow = TRUE, nrow = 3)
#'
#' # Vector field
#' field(p, x, theta)
#'
#' # ... with differentials
#' field(p, x, theta, TRUE)
#'
#' @return If \code{differentials} is \code{FALSE} a vector (of length \code{d}) holding the evaluated value of the field.
#'
#' If \code{differentials} is \code{TRUE} a list is returned with
#' \item{f}{A vector (length \code{d}) holding the evaluated value of the field.}
#' \item{f_dx}{A sparse matrix (dxd) holding the state differential of the field.}
#' \item{f_dparam}{A (list of) sparse matrix (dxp) holding the parameter differential of the field.}
#'
#' @export
field <- function(o, x, param, differentials = FALSE, ...) {
  if (!methods::is(o, "ode")) stop("'o' is not an 'ode'-object.")

  ode_name <- .get_ode_name(o)

  if (!methods::is(param, "list")) {
    param <- list(param)
  }
  param <- lapply(param, as.vector)

  lapply(param, function(para) stopifnot(is.numeric(para), all(is.finite(para))))

  stopifnot(is.numeric(x),
    all(is.finite(x)),
    is.vector(x))

  if (length(param) != length(o$rs) - 1) {
    stop("Length of 'param' is not the same as number of parameter arguments in 'o' (length of o$rs).")
  } else {
    Map(function(p1, p2) {
      if (p1 %% p2 != 0) stop(paste0("Length of param (", p1, ") is not a multiple of number of parameters (", p2, ").")) }, lapply(o$rs[-1], function(re) re$p), lapply(param, length))
  }
  if (length(x) != o$d) stop(paste0("Length of x (", length(x), ") is not equal to number of states (", o$d, ")."))

  if (any(x < o$rs[[1]]$lower)) stop("Some coordinates of 'x' are smaller than 'lower' limit in state.")
  if (any(x > o$rs[[1]]$upper)) stop("Some coordinates of 'x' are larger than 'upper' limit in state.")
  Map(function(para, re) {
    if (any(para < re$lower)) stop("Some coordinates of 'param' are smaller than 'lower' in 'o'.")
    if (any(para > re$upper)) stop("Some coordinates of 'param' are larger than 'upper' in 'o'.")
  }, para = param, re = o$rs[-1])

  Matrix::Matrix()
  ret <- ode_field(ode_struct = o, x = x, param_ = param, differentials = differentials)
  ret$f <- as.vector(ret$f)

  if (!is.null(names(x))) {
    names(ret$f) <- names(x)
    if (differentials) {
      colnames(ret$f_dx) <- rownames(ret$f_dx) <- names(x)
      ret$f_dparam <- lapply(ret$f_dparam, function(r) {
        rownames(r) <- names(x)
        r
      })
    }
  }

  if (differentials) {
    return(ret)
  } else {
    return(ret$f)
  }
}



#' Create 'mak' (Mass Action Kinetics) object
#'
#' @description This function creates an object of class \code{mak} (subclass of \code{ode}), which holds the basic information of the Mass Action Kinetics system in question.
#'
#' @param A The reactant stoichiometric matrix (pxd) containing non-negative values. Here p is the number of parameters and d the number of species.
#' @param B The product stoichiometric matrix (pxd) containing non-negative values. Here p is the number of parameters and d the number of species.
#' @param s \code{\link{solver}} object.
#' @param r An object of class \code{\link{reg}} giving info about how to regularise and bound the rate parameters. If not provided, the default one is used.
#' @param rx0 An object of class \code{\link{reg}} giving info about how to regularise and bound the initial state parameter. If not provided, the default one is used. This default \code{reg} sets \code{fixed = TRUE}, which is generally recommended.
#'
#' @details Mass Action Kinetics is a class of ODE systems, having the following vector field:
#' \deqn{\frac{dx}{dt} = (B - A)^T diag(x^A) k}
#' with \eqn{x^A = (\prod_{i=1}^dx_i^{A_{ji}})_{j=1}^p} and \eqn{k} estimatable and non-negative.
#'
#' @seealso ode, numsolve, field
#' @examples
#' # Michaelis-Menten system
#' A <- matrix(
#' c(1, 1, 0, 0,
#'   0, 0, 1, 0,
#'   0, 0, 1, 0), ncol = 4, byrow = TRUE)
#' B <- matrix(
#' c(0, 0, 1, 0,
#'   1, 1, 0, 0,
#'   1, 0, 0, 1), ncol = 4, byrow = TRUE)
#' k <- c(1, 2, 0.5)
#' x0 <- c(E = 1, S = 4, ES = 0, P = 0)
#' Time <- seq(0, 1, by = .1)
#'
#' m <- mak(A, B)
#'
#' # Solve system
#' numsolve(m, Time, x0, k)
#'
#' # Evaluate field
#' field(m, x0, k)
#'
#' @return An object with S3 class "mak" and "ode".
#' @export
mak <- function(A, B, s = solver(), r = NULL, rx0 = reg("none", lower = 0, upper = Inf, fixed = TRUE)) {
  if (!methods::is(A, "matrix")) stop("A is not a matrix")
  if (!methods::is(B, "matrix")) stop("B is not a matrix")
  rx0$lower <- pmax(0, rx0$lower)

  if (is.null(r)) {
    r <- reg(lower = 0, upper = Inf)
    r$p <- nrow(A)
  } else {

    if (is.null(r$p)) {
      r$p <- nrow(A)
    } else if (r$p != nrow(A)) {
      stop("Number of rows in A does not match 'p' in r.")
    }

    # Check limits
    if (is.null(r$lower)) {
      r$lower <- 0
    } else {
      r$lower[r$lower < 0] <- 0
      if (length(r$lower) != r$p & length(r$lower) != 1) stop("Number of lower limits is not 1 or number of rows in A.")
    }
    if (is.null(r$upper)) {
      upper <- Inf
    } else {
      if (length(r$upper) != r$p & length(r$upper) != 1) stop("Number of upper limits is not 1 or number of rows in A.")
    }
  }
  r$nrow <- 0

  .check(r)
  .check(rx0)

  r <- list(x0 = rx0, rate = r)
  r$x0$nrow <- ncol(A)
  r$x0$p    <- ncol(A)

  ret <- structure(list(A = A, B = B, d = ncol(A), s = s,
    rs = r,
    sparse = methods::is(A, "dgCMatrix")),
    class = c("ode", "mak"))

  ret <- .check(ret)

  return(ret)
}

# Sanity check of mak object
# @param x Object of class \code{\link{mak}}.
# @param ... Additional arguments passed to \code{.check.mak}.
# @inherit .check description return
# @keywords internal
.check.mak <- function(x, ...) {
  Matrix::Matrix() # now "maker" works without prior call of Matrix::Matrix

  with(x, {
    # MAK specifics
    stopifnot(is.matrix(A) || methods::is(A, "dgCMatrix"),
      is.numeric(A) || methods::is(A, "dgCMatrix"),
      all(is.finite(A)),
      all(dim(A) > 0))
    stopifnot(is.matrix(B) || methods::is(B, "dgCMatrix"),
      is.numeric(B) || methods::is(B, "dgCMatrix"),
      all(is.finite(B)),
      all(dim(B) > 0))
    if (!identical(dim(A), dim(B))) {
      stop("Dimensions of 'A' and 'B' must match.")
    }
    if (!identical(class(A), class(B))) {
      stop("Matrix 'A' and 'B' must be of same class (either 'matrix' or 'dgCMatrix').")
    }
    if (any(A < 0)) stop("'A' has negative entries.")

    # MAK specific bounds on x
    if (any(rs[[1]]$lower < 0)) stop("Negatives in 'lower' in reg-object for x0.")
    if (any(rs[[1]]$upper < 0)) stop("Negatives in 'upper' in reg-object for x0.")
    if (any(rs[[2]]$lower < 0)) stop("Negatives in 'lower' in reg-object for rate.")
    if (any(rs[[2]]$upper < 0)) stop("Negatives in 'upper' in reg-object for rate.")

    # nrow
    if (rs[[2]]$nrow != 0) stop("nrow in second 'rs' in mak-object is not 0.")
  })
}


# Exclude parameters in mak object
#
# @description Adapts the \code{\link{mak}} object to the 'exclude' argument in \code{\link{reg}}-objects. Probably irrelevant for most users.
#
# @param x Object of class \code{\link{mak}} to adapt.
# @param ... Additional arguments passed to \code{.exclude.mak}.
# @return An object of class \code{mak}, but with rows of stoichiometric matrices indexed by \code{exclude} removed, \code{p} reduced and \code{lower} and \code{upper} shortened.
# @keywords internal
.exclude.mak <- function(x, ...) {
  exc <- x$rs[[2]]$exclude
  if (!is.null(exc) & length(exc) > 0) {
    x$A <- x$A[-exc, , drop = FALSE]
    x$B <- x$B[-exc, , drop = FALSE]
  }
  x$rs <- lapply(x$rs, .exclude_reg)
  return(x)
}



#' Print 'mak' object
#'
#' @description Prints objects of class \code{mak}, including visualising the reactions.
#'
#' @param x Object of class \code{mak}.
#' @param ... Additional arguments passed to \code{\link{print}}.
#' @seealso mak
#' @examples
#' # Michaelis-Menten system
#' A <- matrix(
#' c(1, 1, 0, 0,
#'   0, 0, 1, 0,
#'   0, 0, 1, 0), ncol = 4, byrow = TRUE)
#' B <- matrix(
#' c(0, 0, 1, 0,
#'   1, 1, 0, 0,
#'   1, 0, 0, 1), ncol = 4, byrow = TRUE)
#' m <- mak(A, B)
#' m
#' @export
print.mak <- function(x, ...) {
  # Print reactions
  A <- .formatcn(x$A); B <- .formatcn(x$B)
  reactions <- paste0(.formatAB(A), " -> ", .formatAB(B))
  cat("Reactions:\n",
    paste0(rep("-", max(nchar(reactions)) + 1), collapse = ""), "\n",
    paste0(reactions, "\n", collapse = " "))
  .print_rest(x)
  cat("\n\nRate Parameter "); print(x$rs[[2]]);
}
.formatcn <- function(A) {
  if (is.null(colnames(A))) {
    colnames(A) <- paste0("X", seq_len(ncol(A)))
  }
  A
}
.formatAB <- function(A, sign = "+") {
  if (sign == "+") {
    nA <- apply(A, 1, function(a) {
      n <- colnames(A)
      n[a == 0] <- ""
      a[a == 1 | a == 0] <- ""
      res <- paste0(a, n)

      # Add plusses
      addplus <- nchar(res) > 0 & cumsum(nchar(res) > 0) > 1
      res[-1] <- paste0((paste0(" ", c(" ", sign), " "))[addplus[-1] + 1], res[-1])
      res
    })
  } else if (sign == "*") {
    nA <- apply(A, 1, function(a) {
      n <- colnames(A)
      a_ <- paste0("^", a)
      n[a == 0] <- ""
      a_[a == 1 | a == 0] <- ""
      res <- paste0(n, a_)

      # Add plusses
      addplus <- nchar(res) > 0 & cumsum(nchar(res) > 0) > 1
      res[-1] <- paste0((paste0(" ", c(" ", sign), " "))[addplus[-1] + 1], res[-1])
      res
    })
  } else {
    return()
  }

  if (ncol(A) == 1) {
    nA <- matrix(nA, ncol = 1)
  } else {
    nA <- t(nA)
  }

  if (nrow(nA) > 1) nA <- apply(nA, 2, format)
  apply(nA, 1, paste0, collapse = "")
}



#' Create 'plk' (Power Law Kinetics) object
#'
#' @description This function creates an object of class \code{plk} (subclass of \code{ode}), which holds the basic information of the Power Law Kinetics system in question.
#'
#' @param A The matrix of powers (pxd). Here d the number of species.
#' @param s \code{\link{solver}} object.
#' @param r An object of class \code{\link{reg}} giving info about how to regularise and bound the rate parameters. If not provided, the default one is used.
#' @param rx0 An object of class \code{\link{reg}} giving info about how to regularise and bound the initial state parameter. If not provided, the default one is used. This default \code{reg} sets \code{fixed = TRUE}, which is generally recommended.
#'
#' @details Power Law Kinetics is a class of ODE systems, having the following vector field:
#' \deqn{\frac{dx}{dt} = \theta x^A}
#' with \eqn{x^A = (\prod_{i=1}^dx_i^{A_{ji}})_{j=1}^p} and \eqn{\theta} an estimatable parameter matrix of dimension dxp. By convention theta will only be reported as a vector (concatinated column-wise).
#'
#' @seealso ode, numsolve, field
#' @examples
#' # Power law system
#' A <- matrix(
#' c(1, 0, 0, 0,
#'   0, 1, 2, 0,
#'   0, 0, 0, 1), ncol = 4, byrow = TRUE)
#' theta <- matrix(
#' c(0, 2, -0.5, 0,
#'   1, 0, 0, 1,
#'   -1, -1, -1, -1), ncol = 3, byrow = TRUE)
#' x0 <- c(X = 1, Y = 4, Z = 0.1, W = 0.1)
#' Time <- seq(0, 1, by = .1)
#'
#' p <- plk(A)
#'
#' # Solve system
#' numsolve(p, Time, x0, theta)
#'
#' # Evaluate field
#' field(p, x0, theta)
#'
#' @return An object with S3 class "plk" and "ode".
#' @export
plk <- function(A, s = solver(), r = NULL, rx0 = reg("none", lower = 0, upper = Inf, fixed = TRUE)) {
  if (!methods::is(A, "matrix")) stop("A is not a matrix")
  if (is.null(r)) {
    r <- reg(lower = -Inf, upper = Inf)
    r$p <- prod(dim(A))
  } else {

    if (is.null(r$p)) {
      r$p <- prod(dim(A))
    } else if (r$p != prod(dim(A))) {
      stop("Number of elements in A does not match 'p' in r.")
    }

    # Check limits
    if (is.null(r$lower)) {
      r$lower <- -Inf
    } else {
      if (length(r$lower) != r$p & length(r$lower) != 1) stop("Number of lower limits is not 1 or number of elments in A.")
    }
    if (is.null(r$upper)) {
      upper <- Inf
    } else {
      if (length(r$upper) != r$p & length(r$upper) != 1) stop("Number of upper limits is not 1 or number of elements in A.")
    }
  }
  r$nrow <- ncol(A)

  .check(r)
  .check(rx0)

  r <- list(x0 = rx0, theta = r)
  r$x0$nrow <- ncol(A)
  r$x0$p    <- ncol(A)

  ret <- structure(list(A = A, d = ncol(A), s = s,
    rs = r,
    sparse = methods::is(A, "dgCMatrix")),
    class = c("ode", "plk"))

  ret <- .check(ret)

  return(ret)
}

# Sanity check of plk object
# @param x Object of class \code{\link{plk}}.
# @param ... Additional arguments passed to \code{.check.plk}.
# @inherit .check description return
# @keywords internal
.check.plk <- function(x, ...) {
  Matrix::Matrix() # now "maker" works without prior call of Matrix::Matrix

  with(x, {
    # PLK specifics
    stopifnot(is.matrix(A) || methods::is(A, "dgCMatrix"),
      is.numeric(A) || methods::is(A, "dgCMatrix"),
      sparse == methods::is(A, "dgCMatrix"),
      all(is.finite(A)),
      all(dim(A) > 0))

    # PLK specific bounds on x

    # Exclude not implemented for plk
    if (!is.null(rs[[2]]$exclude) | length(rs[[2]]$exclude) > 0) stop("Exclude not implemented for 'plk'.")
  })
}


# Exclude parameters in plk object
#
# @description Adapts the \code{\link{plk}} object to the 'exclude' argument in \code{\link{reg}}-objects. Probably irrelevant for most users.
#
# @param x Object of class \code{\link{plk}} to adapt.
# @param ... Additional arguments passed to \code{.exclude.plk}.
# @return An object of class \code{plk}, but with rows of stoichiometric matrix indexed by \code{exclude} removed, \code{p} reduced and \code{lower} and \code{upper} shortened.
# @keywords internal
.exclude.plk <- function(x, ...) {
  # exc <- x$rs[[2]]$exclude
  # if (!is.null(exc)) {
  #   x$A <- x$A[-exc, , drop = FALSE]
  # }
  # x$rs <- lapply(x$rs, .exclude_reg)

  return(x)
}

#' Print 'plk' object
#'
#' @description Prints objects of class \code{plk}, including visualising the powers.
#'
#' @param x Object of class \code{plk}.
#' @param ... Additional arguments passed to \code{\link{print}}.
#' @seealso plk
#' @examples
#' # Power law system
#' A <- matrix(
#' c(1, 0, 0, 0,
#'   0, 1, 2, 0,
#'   0, 0, 0, 0), ncol = 4, byrow = TRUE)
#' p <- plk(A)
#' p
#' @export
print.plk <- function(x, ...) {
  # Print reactions
  A <- .formatcn(x$A);
  reactions <- .formatAB(A, sign = "*")
  cat("Powers:\n",
    paste0(rep("-", max(nchar(reactions)) + 1), collapse = ""), "\n",
    paste0(reactions, "\n", collapse = " "))
  .print_rest(x)
  cat("\n\nTheta "); print(x$rs[[2]]);
}


#' Create 'ratmak' (Rational Mass Action Kinetics) object
#'
#' @description This function creates an object of class \code{ratmak} (subclass of \code{ode}), which holds the basic information of the Rational Mass Action Kinetics system in question.
#'
#' @param A The matrix of powers (bxd). Here d the number of species.
#' @param C The stoichiometric coefficient matrix (rxd). Here r is the number of reactions.
#' @param s \code{\link{solver}} object.
#' @param r1 An object of class \code{\link{reg}} giving info about how to regularise and bound the first parameter. If not provided, the default one is used.
#' @param r2 An object of class \code{\link{reg}} giving info about how to regularise and bound the second  parameters. If not provided, the default one is used.
#' @param rx0 An object of class \code{\link{reg}} giving info about how to regularise and bound the initial state parameter. If not provided, the default one is used. This default \code{reg} sets \code{fixed = TRUE}, which is generally recommended.
#'
#' @details Rational Mass Action Kinetics is a class of ODE systems, having the following vector field:
#' \deqn{\frac{dx}{dt} = C^T * (\theta_1 * x^A) / (1 + \theta_2 * x^A)}
#' with \eqn{x^A = (\prod_{i=1}^dx_i^{A_{ji}})_{j=1}^b}, \eqn{\theta_1} and \eqn{\theta_2} estimatable parameter matrices of dimension rxb. \eqn{\theta_2} is restricted to non-negative values. The fraction is entry-wise. By convention the theta's will only be reported as vectors (concatinated column-wise).
#'
#' @seealso ode, numsolve, field
#'
#' @return An object with S3 class "ratmak" and "ode".
#'
#' @examples
#' # Rational mass action kinetics
#' A <- matrix(
#' c(1, 0, 0, 0,
#'   0, 1, 2, 0,
#'   1, 0, 0, 1), ncol = 4, byrow = TRUE)
#' x0 <- c(X = 1, Y = 4, Z = 0.1, W = 0.1)
#' time <- seq(0, 1, by = .1)
#'
#' rmak <- ratmak(A, diag(4))
#'
#' # Solve system
#' numsolve(o = rmak, time = time, x0 = x0,
#'   param = list(theta1 = t(A * 1:3),
#'       theta2 = t((A + 1) * 3:1)))
#'
#' # Evaluate field
#' field(rmak, x0,
#'   param = list(theta1 = t(A * 1:3),
#'       theta2 = t((A + 1) * 3:1)))
#'
#' @export
ratmak <- function(A, C, s = solver(), r1 = NULL, r2 = NULL, rx0 = reg("none", lower = 0, upper = Inf, fixed = TRUE)) {
  if (!methods::is(A, "matrix")) stop("A is not a matrix")
  if (!methods::is(C, "matrix")) stop("C is not a matrix")
  if (ncol(A) != ncol(C)) {
    stop("'A' and 'C' must have same number of columns.")
  }

  r1 <- .prepare_ratmak(re = r1, b = nrow(A), r = nrow(C))
  r2 <- .prepare_ratmak(re = r2, b = nrow(A), r = nrow(C))

  # Check r1 limits
  if (is.null(r1$lower)) {
    r1$lower <- -Inf
  } else {
    if (length(r1$lower) != r1$p & length(r1$lower) != 1) stop("Number of lower limits is not 1 or r*b.")
  }
  if (is.null(r1$upper)) {
    r1$upper <- Inf
  } else {
    if (length(r1$upper) != r1$p & length(r1$upper) != 1) stop("Number of upper limits is not 1 or r*b.")
  }

  # Check r2 limits
  if (is.null(r2$lower)) {
    r2$lower <- 0
  } else {
    if (length(r2$lower) != r2$p & length(r2$lower) != 1) stop("Number of lower limits is not 1 or r*b.")
    r2$lower[r2$lower < 0] <- 0
  }
  if (is.null(r2$upper)) {
    r2$upper <- Inf
  } else {
    r2$upper[r2$upper < 0] <- 0
    if (length(r2$upper) != r2$p & length(r2$upper) != 1) stop("Number of upper limits is not 1 or r*b.")
  }

  .check(r1)
  .check(r2)
  .check(rx0)

  r <- list(x0 = rx0, theta1 = r1, theta2 = r2)
  r$x0$nrow <- ncol(A)
  r$x0$p    <- ncol(A)

  ret <- structure(list(A = A, C = C, d = ncol(A), s = s,
    rs = r,
    sparseA = methods::is(A, "dgCMatrix"),
    sparseC = methods::is(C, "dgCMatrix")),
    class = c("ode", "ratmak"))

  ret <- .check(ret)

  return(ret)
}

# Prepares (and updates) the reg objects for ratmak
.prepare_ratmak <- function(re, b, r) {
  if (is.null(re)) {
    re <- reg(lower = -Inf, upper = Inf)
    re$p <- r * b
  } else if (methods::is(re, "reg")) {
    if (is.null(re$p)) {
      re$p <- r * b
    } else if (re$p != r * b) {
      stop("r*b does not match 'p' in reg-object.")
    }
  } else {
    stop("Supplied r is not reg-object.")
  }

  re$nrow <- r
  return(re)
}



# Sanity check of plk object
# @param x Object of class \code{\link{plk}}.
# @param ... Additional arguments passed to \code{.check.plk}.
# @inherit .check description return
# @keywords internal
.check.ratmak <- function(x, ...) {
  Matrix::Matrix() # now "maker" works without prior call of Matrix::Matrix

  with(x, {
    # Ratmak specifics
    stopifnot(is.matrix(A) || methods::is(A, "dgCMatrix"),
      is.numeric(A) || methods::is(A, "dgCMatrix"),
      (methods::is(A, "dgCMatrix") && sparseA) || (is.matrix(A) && !sparseA),
      all(is.finite(A)),
      all(dim(A) > 0))
    stopifnot(is.matrix(C) || methods::is(C, "dgCMatrix"),
      is.numeric(C) || methods::is(C, "dgCMatrix"),
      (methods::is(C, "dgCMatrix") && sparseC) || (is.matrix(C) && !sparseC),
      all(is.finite(C)),
      all(dim(C) > 0))
    if (ncol(A) != ncol(C)) {
      stop("Number of columns of 'A' and 'C' must match.")
    }
    if (any(A < 0)) stop("'A' has negative entries.")

    # ratmak specific bounds on x
    if (any(rs[[3]]$lower < 0)) stop("Negatives in 'lower' in reg-object for theta2.")
    if (any(rs[[3]]$upper < 0)) stop("Negatives in 'upper' in reg-object for theta2.")

    # nrow
    if (rs[[2]]$nrow != nrow(C)) stop(paste0("nrow in second 'rs' (", rs[[2]]$nrow , ") in mak-object is not r (", nrow(C) ,")."))
    if (rs[[3]]$nrow != nrow(C))  stop(paste0("nrow in second 'rs' (", rs[[3]]$nrow , ") in mak-object is not r (", nrow(C) ,")."))
  })
}


# Exclude parameters in ratmak object
#
# @description Adapts the \code{\link{ratmak}} object to the 'exclude' argument in \code{\link{reg}}-objects. Probably irrelevant for most users.
#
# @param x Object of class \code{\link{ratmak}} to adapt.
# @param ... Additional arguments passed to \code{.exclude.ratmak}.
# @return An object of class \code{ratmak}, but with rows of stoichiometric matrix indexed by \code{exclude} removed, \code{p} reduced and \code{lower} and \code{upper} shortened.
# @keywords internal
.exclude.ratmak <- function(x, ...) {
  # exc <- x$rs[[2]]$exclude
  # if (!is.null(exc)) {
  #   x$A <- x$A[-exc, , drop = FALSE]
  # }
  # x$rs <- lapply(x$rs, .exclude_reg)

  return(x)
}


#' Print 'ratmak' object
#'
#' @description Prints objects of class \code{ratmak}, and presents the powers.
#'
#' @param x Object of class \code{ratmak}.
#' @param ... Additional arguments passed to \code{\link{print}}.
#' @seealso ratmak
#'
#' @examples
#' # Rational mass action kinetics
#' A <- matrix(
#' c(1, 0, 0, 0,
#'   0, 1, 2, 0,
#'   1, 0, 0, 1), ncol = 4, byrow = TRUE)
#' x0 <- c(X = 1, Y = 4, Z = 0.1, W = 0.1)
#' time <- seq(0, 1, by = .1)
#'
#' rmak <- ratmak(A, diag(4))
#' rmak
#' @export
print.ratmak <- function(x, ...) {
  # Print reactions
  A <- .formatcn(x$A); C <- .formatcn(x$C)
  reactions <- .formatAB(A, sign = "*")
  cat("Powers:\n",
    paste0(rep("-", max(nchar(reactions)) + 1), collapse = ""), "\n",
    paste0(reactions, "\n", collapse = " "))
  .print_rest(x)
  cat("\n\nTheta1 "); print(x$rs[[2]]);
  cat("\n\nTheta2 "); print(x$rs[[3]]);
}




#' Create 'rlk' (Rational Law Kinetics) object
#'
#' @description This function creates an object of class \code{rlk} (subclass of \code{ode}), which holds the basic information of the Rational Law Kinetics system in question.
#'
#' @param A The matrix of powers (pxd). Here d the number of species.
#' @param B The matrix of powers (pxd). Here d the number of species.
#' @param s \code{\link{solver}} object.
#' @param r An object of class \code{\link{reg}} giving info about how to regularise and bound the rate parameters. If not provided, the default one is used.
#' @param rx0 An object of class \code{\link{reg}} giving info about how to regularise and bound the initial state parameter. If not provided, the default one is used. This default \code{reg} sets \code{fixed = TRUE}, which is generally recommended.
#'
#' @details Rational Law Kinetics is a class of ODE systems, having the following vector field:
#' \deqn{\frac{dx}{dt} = \theta (x^A / (1 + x^B))}
#' with \eqn{x^A = (\prod_{i=1}^dx_i^{A_{ji}})_{j=1}^p} and \eqn{\theta} an estimatable parameter matrix of dimension dxp. By convention theta will only be reported as a vector (concatinated column-wise).
#'
#' @seealso ode, numsolve, field
#'
#' @return An object with S3 class "rlk" and "ode".
#'
#' @examples
#' # Rational law kinetics
#' A <- matrix(
#' c(1, 0, 0, 0,
#'   0, 1, 2, 0,
#'   0, 0, 0, 1), ncol = 4, byrow = TRUE)
#' theta <- matrix(
#' c(0, 2, -0.5, 0,
#'   1, 0, 0, 1,
#'   -1, -1, -1, -1), ncol = 3, byrow = TRUE)
#' x0 <- c(X = 1, Y = 4, Z = 0.1, W = 0.1)
#' time <- seq(0, 1, by = .1)
#' r <- rlk(A, A[c(2, 1, 3), ])
#'
#' # Solve system
#' numsolve(o = r, time = time, x0 = x0, param = theta)
#'
#' # Evaluate field
#' field(r, x0, theta)
#'
#' @export
rlk <- function(A, B, s = solver(), r = NULL, rx0 = reg("none", lower = 0, upper = Inf, fixed = TRUE)) {
  if (!methods::is(A, "matrix")) stop("A is not a matrix")
  if (!methods::is(B, "matrix")) stop("B is not a matrix")
  if (!all(dim(A) == dim(B))) {
    stop("Dimension of A and B do not match.")
  }
  if (is.null(r)) {
    r <- reg(lower = -Inf, upper = Inf)
    r$p <- prod(dim(A))
  } else {

    if (is.null(r$p)) {
      r$p <- prod(dim(A))
    } else if (r$p != prod(dim(A))) {
      stop("Number of elements in A does not match 'p' in r.")
    }

    # Check limits
    if (is.null(r$lower)) {
      r$lower <- -Inf
    } else {
      if (length(r$lower) != r$p & length(r$lower) != 1) stop("Number of lower limits is not 1 or number of elments in A.")
    }
    if (is.null(r$upper)) {
      upper <- Inf
    } else {
      if (length(r$upper) != r$p & length(r$upper) != 1) stop("Number of upper limits is not 1 or number of elements in A.")
    }
  }
  r$nrow <- ncol(A)

  .check(r)
  .check(rx0)

  r <- list(x0 = rx0, theta = r)
  r$x0$nrow <- ncol(A)
  r$x0$p    <- ncol(A)

  ret <- structure(list(A = A, B = B, d = ncol(A), s = s,
                        rs = r,
                        sparseA = methods::is(A, "dgCMatrix"),
                        sparseB = methods::is(B, "dgCMatrix")),
                   class = c("ode", "rlk"))

  ret <- .check(ret)

  return(ret)
}

# Sanity check of rlk object
# @param x Object of class \code{\link{rlk}}.
# @param ... Additional arguments passed to \code{.check.rlk}.
# @inherit .check description return
# @keywords internal
.check.rlk <- function(x, ...) {
  Matrix::Matrix()

  with(x, {
    # RLK specifics
    stopifnot(is.matrix(A) || methods::is(A, "dgCMatrix"),
              is.numeric(A) || methods::is(A, "dgCMatrix"),
              sparseA == methods::is(A, "dgCMatrix"),
              all(is.finite(A)),
              all(dim(A) > 0))
    stopifnot(is.matrix(B) || methods::is(B, "dgCMatrix"),
              is.numeric(B) || methods::is(B, "dgCMatrix"),
              sparseB == methods::is(B, "dgCMatrix"),
              all(is.finite(B)),
              all(dim(B) > 0))

    # PLK specific bounds on x

    # Exclude not implemented for plk
    if (!is.null(rs[[2]]$exclude) | length(rs[[2]]$exclude) > 0) stop("Exclude not implemented for 'rlk'.")
  })
}


# Exclude parameters in rlk object
#
# @description Adapts the \code{\link{rlk}} object to the 'exclude' argument in \code{\link{reg}}-objects. Probably irrelevant for most users.
#
# @param x Object of class \code{\link{rlk}} to adapt.
# @param ... Additional arguments passed to \code{.exclude.rlk}.
# @return An object of class \code{rlk}, but with rows of stoichiometric matrix indexed by \code{exclude} removed, \code{p} reduced and \code{lower} and \code{upper} shortened.
# @keywords internal
.exclude.rlk <- function(x, ...) {
  # exc <- x$rs[[2]]$exclude
  # if (!is.null(exc)) {
  #   x$A <- x$A[-exc, , drop = FALSE]
  # }
  # x$rs <- lapply(x$rs, .exclude_reg)

  return(x)
}


#' Print 'rlk' object
#'
#' @description Prints objects of class \code{rlk} and presents the fractions.
#'
#' @param x Object of class \code{rlk}.
#' @param ... Additional arguments passed to \code{\link{print}}.
#' @seealso rlk
#'
#' @examples
#' # Rational law kinetics
#' A <- matrix(
#' c(1, 0, 0, 0,
#'   0, 1, 2, 0,
#'   0, 0, 0, 1), ncol = 4, byrow = TRUE)
#' r <- rlk(A, A[c(2, 1, 3), ])
#' r
#' @export
print.rlk <- function(x, ...) {
  # Print reactions
  A <- .formatcn(x$A); B <- .formatcn(x$B)
  Aside <- .formatAB(A, sign = "*")
  Aside[!apply(x$A != 0, 1, any)] <- "1"
  Aside <- format(Aside, justify = "right")
  Bside <- paste0(" / 1 + ", .formatAB(B, sign = "*"))
  Bside[!apply(x$B != 0, 1, any)] <- ""
  reactions <- paste0(Aside, Bside)
  cat("Fractions:\n",
    paste0(rep("-", max(nchar(reactions)) + 1), collapse = ""), "\n",
    paste0(reactions, "\n ", collapse = ""))
  .print_rest(x)
  cat("\n\nTheta "); print(x$rs[[2]]);
}

.print_rest <- function(x) {
  # Print solver
  cat("\n"); print(x$s)

  # Print reg
  cat("\n\nInitial State "); print(x$rs[[1]]);
}

