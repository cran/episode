###-------------###
###  Test aim   ###
###-------------###

library(episode)
context("aim object")

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
x0 <- c(X = 1, Y = 2, Z = .1, W = 2)

set.seed(1234)
times_1 <- c(seq(0, 1, by = .1))
curve_1 <- numsolve(m, times_1, x0, k)
y_safe1 <- curve_1
y_safe1[, -1] <- y_safe1[, -1] + rnorm(prod(dim(y_safe1[, -1])), sd = .1)

times_2 <- c(seq(0, .1, by = .01))
context <- cbind(1, c(1, 2, 1, 2, 1, 2, 1))
curve_2 <- rbind(numsolve(m, times_1, x0, k * context[, 1]),
  numsolve(m, times_2, x0, k * context[, 2]))
y_safe2 <- curve_2
y_safe2[, -1] <- y_safe2[, -1] + rnorm(prod(dim(y_safe2[, -1])), sd = .1)

curve_x <- rbind(numsolve(m, sort(runif(42)), x0, k * context[, 1]),
  numsolve(m, sort(runif(42, min = -.01, max = 0.11)), x0, k * context[, 2]))

test_that("imd", {
  w <- matrix(runif(prod(dim(y_safe2) - c(0, 1)), min = .9, max = 1.1), nrow = nrow(y_safe2))
  v <- runif(length(k), min = .9, max = 1.1)

  m <- mak(A, B,
    r = reg(reg_type = "l1", step_max = 200, a = 3.25,
      upper = seq_along(k) / .1, backtrack_max = 25,
      lambda_factor = .2, exclude = c(6),
      screen = TRUE, penalty_factor = v,
      contexts = cbind(1, rep(c(1, 2), length.out = 7)),
      scales = c(10.1, 1.2, 1.2, 1.2, rep(1.2, 2), 1.2),
      exact_gradient = TRUE),
    rx0 = reg(reg_type = "none", fixed = TRUE))
  ro <- opt(y_safe2, weights = w, nlambda = 44, lambda_min_ratio = 0.0012, tol = 1e-6)
  foo <- imd(m, ro)

  colnames(foo$xout) <- colnames(y_safe2)
  expect_equal(foo$xout, y_safe2)
  expect_equal(as.vector(t(w[c(2:11, 13:22),] + w[c(1:10, 12:21), ]) / 2), foo$W)

  rownames(foo$X0) <- colnames(foo$xout)[-1]
  expect_equal(t(foo$xout[c(TRUE, diff(foo$xout[, 1]) < 0), -1]), foo$X0)

  expect_equal(as.vector(t(foo$xout[c(2:11, 13:22),-1] - foo$xout[c(1:10, 12:21),-1])), foo$Y)


  # since linear it should not use params
  foo2 <- imd(m, ro, params = list(k = seq_along(k)))
  colnames(foo2$xout) <- colnames(y_safe2)
  rownames(foo2$X0) <- colnames(foo2$xout)[-1]
  expect_identical(foo, foo2)

  # custom x
  foo2 <- imd(m, ro, x = curve_x)
  expect_equal(t(foo2$xout[c(TRUE, diff(foo2$xout[, 1]) < 0), -1]), foo2$X0)
  expect_equal(as.vector(t(foo2$xout[c(2:11, 13:22),-1] - foo2$xout[c(1:10, 12:21),-1])), foo2$Y)
})


test_that("Adapts", {
  #  Minimal
  o <- mak(A, B)
  ro <- opt(y_safe2)
  a1 <- aim(o = o, op = ro, adapts = c("scales"))
  expect_true(is.null(a1$o$rs[[2]]$penalty_factor))
  a2 <- aim(o = o, op = ro, adapts = c("scales", "penalty_factor", "weights"))
  expect_identical(a1$o$rs[[2]]$scales, a2$o$rs[[2]]$scales)

  # maximal
  w <- matrix(runif(prod(dim(y_safe2) - c(0, 1)), min = .9, max = 1.1), nrow = nrow(y_safe2))
  v <- runif(length(k), min = .9, max = 1.1)
  m <- mak(A, B,
    r = reg(reg_type = "scad", step_max = 200, a = 3.25,
      upper = seq_along(k) / .1, backtrack_max = 25,
      lambda_factor = 1, exclude = c(6),
      screen = TRUE, penalty_factor = v,
      contexts = cbind(1, rep(c(1, 2), length.out = 7)),
      scales = c(10.1, 1.2, 1.2, 1.2, rep(1.2, 2), 1.2),
      exact_gradient = TRUE),
    rx0 = reg(reg_type = "none", fixed = TRUE))
  ro <- opt(y_safe2, weights = w, nlambda = 44, lambda_min_ratio = 0.0012, tol = 1e-6)
  a2 <- aim(o = m, op = ro, adapts = NULL)


  cc <- cattle(ode_struct = m, opt_struct = ro, sc_ = list(NULL),
    param_ = NULL, x = curve_x)
  cc$params
  cc$lambda


  colnames(a2$x) <- colnames(y_safe2)
  expect_identical(a2$x, y_safe2)

  # custom x
  w <- matrix(runif(prod(dim(y_safe2) - c(0, 1)), min = .9, max = 1.1), nrow = nrow(y_safe2))
  v <- runif(length(k), min = .9, max = 1.1)
  m <- mak(A, B,
    r = reg(reg_type = "scad", step_max = 200, a = 3.25,
      upper = seq_along(k) / .1, backtrack_max = 25,
      lambda_factor = .2, exclude = c(2),
      screen = TRUE, penalty_factor = v,
      contexts = cbind(1, rep(c(1, 2), length.out = 7)),
      scales = c(10.1, 1.2, 1.2, 1.2, rep(1.2, 2), 1.2),
      exact_gradient = TRUE),
    rx0 = reg(reg_type = "none", fixed = TRUE))
  ro <- opt(y_safe2, weights = w, nlambda = 44, lambda_min_ratio = 0.0012, tol = 1e-6)
  a2 <- aim(o = m, op = ro, adapts = c("scales", "penalty_factor"), x = curve_x)

  expect_true(all(a2$params[[1]][2, ] == 0))

  p <- plk(B,
    r = reg(reg_type = "none", step_max = 200, a = 3.25),
    rx0 = reg(reg_type = "none", fixed = TRUE))
  ro <- opt(y_safe2, weights = w, nlambda = 44, lambda_min_ratio = 0.0012, tol = 1e-6)
  a3 <- aim(o = p, op = ro, adapts = c("scales", "penalty_factor", "weights"))
})


test_that("Reg types", {
  o <- mak(A, B, r = reg(reg_type = "l2"))
  a <- aim(o = o, op = opt(y_safe1))

  o <- mak(A, B, r = reg(reg_type = "elnet"))
  a2 <- aim(o = o, op = opt(y_safe1))

  o <- mak(A, B, r = reg(reg_type = "mcp"))
  a2 <- aim(o = o, op = opt(y_safe1))

  o <- mak(A, B, r = reg(reg_type = "scad"))
  a2 <- aim(o = o, op = opt(y_safe1))

})




test_that("Integral Matching Design", {
  o <- mak(A, B, r = reg(reg_type = "l2"))
  im <- imd(o = o, op = opt(y_safe2))
  expect_equal(sum((im$xout - y_safe2)^2), 0)
  expect_identical(
    t(im$xout[c(TRUE, (diff(im$xout[, 1]) < 0)), -1]),
    im$X0)

  im <- imd(o = o, op = opt(y_safe2), x = curve_x)
  expect_identical(
    t(im$xout[c(TRUE, (diff(im$xout[, 1]) < 0)), -1]),
    im$X0)
})





test_that("rodeo.aim", {
  # Estimate x0
  o <- mak(A, B, rx0 = reg(reg_type = "l2", fixed = FALSE))
  ro <- opt(y_safe2, nlambda = 15)
  a1 <- aim(o = o, op = ro, adapts = c("scales"))
  r1 <- rodeo(a1, adjusts = c("scales", "weights", "lambda"))

  expect_true(all(is.finite(r1$dfs)))
  expect_true(all(is.finite(r1$x0s)))
  expect_true(all(is.finite(r1$params[[1]])))


  a2 <- aim(o = o, op = ro, adapts = c("scales", "penalty_factor", "weights"))
  a2$params
  a2$x0s
  r2 <- rodeo(a2)
  r2$params
  r2$dfs
  r2$x0s

  expect_true(all(is.finite(r2$dfs)))
  expect_true(all(is.finite(r2$x0s)))
  expect_true(all(is.finite(r2$params[[1]])))


  # custom x
  w <- matrix(runif(prod(dim(y_safe2) - c(0, 1)), min = .9, max = 1.1), nrow = nrow(y_safe2))
  v <- runif(length(k), min = .9, max = 1.1)
  m <- mak(A, B,
    r = reg(reg_type = "scad", step_max = 200, a = 3.25,
      upper = seq_along(k) / .1, backtrack_max = 25,
      lambda_factor = .21, exclude = c(2),
      screen = TRUE, penalty_factor = v,
      contexts = cbind(1, rep(c(1, 2), length.out = 7)),
      scales = c(10.1, 1.2, 1.2, 1.2, rep(1.2, 2), 1.2),
      exact_gradient = TRUE),
    rx0 = reg(reg_type = "none", fixed = TRUE))
  ro <- opt(y_safe2, weights = w, lambda_min_ratio = 0.000012, tol = 1e-6, nlambda = 15)
  a2 <- aim(o = m, op = ro, adapts = c("scales", "penalty_factor"), x = curve_x)
  a2$params
  r2 <- rodeo(a2, adjusts = c("scales", "weights", "lambda"))
  r2$params

  expect_true(all(is.finite(r2$dfs)))
  expect_true(all(is.finite(r2$x0s)))
  expect_true(all(is.finite(r2$params[[1]])))



  # plk
  w <- matrix(runif(prod(dim(y_safe2) - c(0, 1)), min = .9, max = 1.1), nrow = nrow(y_safe2))
  v <- runif(length(A), min = .9, max = 1.1)
  p <- plk(A,
    r = reg(reg_type = "scad", step_max = 200, a = 3.25,
      upper = seq_along(A) / .1, backtrack_max = 25,
      lambda_factor = .21,
      screen = TRUE, penalty_factor = v,
      contexts = cbind(1, rep(c(1, 2), length.out = length(A))),
      scales = c(10.1, rep(1.2, length(A) - 1)),
      exact_gradient = TRUE),
    rx0 = reg(reg_type = "none", fixed = TRUE))
  ro <- opt(y_safe2, weights = w, lambda_min_ratio = 0.0012, tol = 1e-6, nlambda = 15)
  a2 <- aim(o = p, op = ro, adapts = c("scales", "penalty_factor"), x = curve_x)
  a2$params
  r2 <- rodeo(a2, adjusts = c("scales", "weights", "lambda"))
  r2$params

  expect_true(all(is.finite(r2$dfs)))
  expect_true(all(is.finite(r2$x0s)))
  expect_true(all(is.finite(r2$params[[1]])))


  # ratmak
  A <- matrix(c(1, 1, 0, 0,
    0, 0, 1, 0,
    0, 0, 1, 0), ncol = 4, byrow = TRUE)
  CC <- matrix(c(0, 0, 1, 0,
    1, 1, 0, 0,
    1, 0, 0, 1), ncol = 4, byrow = TRUE) - A

  w <- matrix(runif(prod(dim(y_safe2) - c(0, 1)), min = .9, max = 1.1), nrow = nrow(y_safe2))
  v <- runif(nrow(A) * nrow(CC), min = .9, max = 1.1)

  rat <- ratmak(A, CC,
    r1 = reg(reg_type = "scad", step_max = 200, a = 3.25,
      upper = seq_along(v) / 1.1, backtrack_max = 25,
      step_screen = 3,
      lambda_factor = .8,
      screen = TRUE, penalty_factor = v,
      contexts = cbind(1, rep(c(1, 2), length.out = length(v))),
      scales = c(1.2, rep(1.2, length(v) - 1)),
      exact_gradient = TRUE),
    r2 = reg(reg_type = "mcp", step_max = 200, a = 3.25,
      upper = seq_along(v) / 1.1, backtrack_max = 25,
      step_screen = 5,
      lambda_factor = .8,
      screen = TRUE, penalty_factor = v,
      contexts = cbind(1, rep(c(1, 2), length.out = length(v))),
      scales = c(10.1, rep(1.2, length(v) - 1)),
      exact_gradient = TRUE),
    rx0 = reg(reg_type = "none", fixed = TRUE))
  ro <- opt(y_safe2, weights = w, nlambda = 19, lambda_min_ratio = 0.000012, tol = 1e-6)
  a2 <- aim(o = rat, op = ro, adapts = c("scales", "penalty_factor"), x = curve_x)
  a2$params
  r2 <- rodeo(a2, adjusts = c("scales", "weights", "lambda"))
  r2$params

  expect_true(all(is.finite(r2$dfs)))
  expect_true(all(is.finite(r2$x0s)))
  expect_true(all(is.finite(r2$params[[1]])))


})


