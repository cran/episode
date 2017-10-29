###------------###
###  Test ODE  ###
###------------###

library(episode)
context("ode object")


A <- matrix(c(1, 1, 0, 0,
  0, 0, 1, 0,
  0, 0, 1, 0), ncol = 4, byrow = TRUE)
B <- matrix(c(0, 0, 1, 0,
  1, 1, 0, 0,
  1, 0, 0, 1), ncol = 4, byrow = TRUE)


test_that("Constructor", {
  m <- mak(A, B)

  expect_identical(class(m), c("ode", "mak"))

  expect_equal(m$s, solver())
  expect_equal(m$d, ncol(A))
  expect_equal(m$rs[[2]]$p, nrow(A))

  expect_error(mak(A, B[, -1]), "Dimensions of 'A' and 'B' must match.")
  expect_error(mak(A, B[-1, ]), "Dimensions of 'A' and 'B' must match.")
  expect_error(mak(A - 1, B), "'A' has negative entries.")


  m <- mak(A, B, r = reg(lower = seq_len(nrow(A)), upper = Inf))
  expect_equal(m$rs[[2]]$upper, rep(Inf, nrow(A)))
  expect_equal(m$d, ncol(A))
  expect_equal(m$rs[[2]]$p, nrow(A))


  p <- plk(A)
  expect_identical(class(p), c("ode", "plk"))

  expect_equal(p$s, solver())
  expect_equal(p$d, ncol(A))
  expect_equal(p$rs[[2]]$p, prod(dim(A)))

  p <- plk(A, r = reg(lower = -Inf, upper = seq_along(A)))
  expect_equal(p$rs[[2]]$lower, rep(-Inf, p$rs[[2]]$p))
  expect_equal(p$d, ncol(A))
  expect_equal(p$rs[[2]]$p, prod(dim(A)))
})



test_that("numsolve", {
  A <- matrix(c(1, 1, 0, 0,
    0, 0, 1, 0,
    0, 0, 1, 0), ncol = 4, byrow = TRUE)
  B <- matrix(c(0, 0, 1, 0,
    1, 1, 0, 0,
    1, 0, 0, 1), ncol = 4, byrow = TRUE)

  # Should all give same system by tweeking parameters
  m <- mak(A, B)
  p <- plk(A)
  r <- rlk(A, matrix(0, nrow = 3, ncol = 4))
  rmak <- ratmak(A, diag(4))

  time <- seq(0, 1, by = 0.01)
  x0 <- c(X = 1, Y = 2, Z = 0, W = 1)

  expect_equal(numsolve(o = m, time = time, x0 = x0, param = 1:3),
    numsolve(o = p, time = time, x0 = x0, param = t((B - A) * 1:3)))
  expect_equal(numsolve(o = m, time = time, x0 = x0, param = 1:3),
    numsolve(o = r, time = time, x0 = x0, param = 2 * t((B - A) * 1:3)))
  expect_equal(numsolve(o = m, time = time, x0 = x0, param = 1:3),
    numsolve(o = rmak, time = time, x0 = x0,
      param = list(theta1 = t((B - A) * 1:3), theta2 = matrix(0, ncol = 3, nrow = 4))))

  expect_equal(numsolve(o = m, time = time, x0 = x0, param = 2:0),
    numsolve(o = p, time = time, x0 = x0, param = t((B - A) * 2:0)))

  # Extra check that rlk and ratmak works for nontrivial systems
  r <- rlk(A, B)
  numsolve(o = r, time = time, x0 = x0, param = 2 * t((B - A) * 1:3), approx_sensitivity = FALSE)

  rmak <- ratmak(A, t(B) %*% B)
  numsolve(o = rmak, time = time, x0 = x0,
    param = list(theta1 = t((B - A) * 1:3), theta2 = t(B * 2:0)), approx_sensitivity = TRUE)
})



test_that("field", {
  A <- matrix(c(1, 1, 0, 0,
    0, 0, 1, 0,
    0, 0, 1, 0), ncol = 4, byrow = TRUE)
  B <- matrix(c(0, 0, 1, 0,
    1, 1, 0, 0,
    1, 0, 0, 1), ncol = 4, byrow = TRUE)

  m <- mak(A, B)
  p <- plk(A)

  x <- c(X = 1, Y = 2, Z = 0, W = 1)
  k <- c(k1 = 1, k2 = 2, k3 = .5)
  fm <- field(o = m, x = x, param = k, differentials = TRUE)
  fp <- field(o = p, x = x, param = t((B - A) * k), differentials = TRUE)
  expect_equal(fm$f, fp$f)
  expect_equal(fm$f_dx, fp$f_dx)

  x <- c(X = 1, Y = 2, Z = 0.25, W = 1)
  fm <- field(o = m, x = x, param = k, differentials = TRUE)
  fp <- field(o = p, x = x, param = t((B - A) * k), differentials = TRUE)
  expect_equal(fm$f, fp$f)
  expect_equal(sum((fm$f_dparam[[1]] - t((B - A) * apply(p$A, 1, function(a) prod(x^a))))^2), 0)
  expect_equal(sum((fp$f_dparam[[1]] -
    do.call(cbind, sapply(apply(p$A, 1, function(a) prod(x^a)), diag, nrow = p$d, ncol = p$d, simplify = FALSE))))^2, 0)
  expect_equal(fm$f_dx, fp$f_dx)


  A <- matrix(c(1, 1, 0, 0,
    0, 1, 1, 0,
    0, 0, 1, 1,
    2, 1, 0, 0,
    2, 0, 0, 1,
    0, 2, 2, 1,
    0, 1, 1, 2), ncol = 4, byrow = TRUE)
  B <- A[c(4, 3, 2, 1, 5, 7, 6), c(1, 3, 2, 4)]

  m <- mak(A, B)
  p <- plk(A)

  k <- c(k1 = 1, k2 = 2, k3 = .5, k4 = 0, k5 = 1, k6 = 0, k7 = .1)
  x <- c(X = 1, Y = 2, Z = .1, W = 2)
  f <- field(o = m, x = x, param = k, differentials = TRUE)


  expect_equal(f$f,
    setNames(as.vector(t(B - A) %*% (apply(A, 1, function(a) prod(x^a)) * k)), names(x)))
  f_dx <- t((B - A) * k) %*% t(t(apply(A, 1, function(a) prod(x^a)) * A) / x)
  rownames(f_dx) <- colnames(f_dx) <- names(x)
  expect_equal(as.matrix(f$f_dx), f_dx)
  f_dk <- t((B - A) * apply(A, 1, function(a) prod(x^a)))
  rownames(f_dk) <- names(x)
  # colnames(f_dk) <- names(k)
  expect_equal(as.matrix(f$f_dparam[[1]]), f_dk)


  x <- c(X = 1, Y = 2, Z = 0, W = 2)
  f <- field(o = m, x = x, param = k, differentials = TRUE)

  expect_equal(f$f,
    setNames(as.vector(t(B - A) %*% (apply(A, 1, function(a) prod(x^a)) * k)), names(x)))
  f_dk <- t((B - A) * apply(A, 1, function(a) prod(x^a)))
  rownames(f_dk) <- names(x)
  # colnames(f_dk) <- names(k)
  expect_equal(as.matrix(f$f_dparam[[1]]), f_dk)

})

