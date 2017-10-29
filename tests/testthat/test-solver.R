###---------------###
###  Test solver  ###
###---------------###

library(episode)
context("solver object")


test_that("Constructor", {
  expect_error(solver("fretz"), "'name' of solver not recognised. Must be one of: rk23, bs23, dp45, rkf45")
  expect_error(solver(step_max = 0), "step_max >= 1 is not TRUE")
  expect_error(solver(tol = 0), "tol > 0 is not TRUE")
})


test_that("Loss", {
  A <- matrix(c(1, 1, 0, 0,
    0, 0, 1, 0,
    0, 0, 1, 0), ncol = 4, byrow = TRUE)
  B <- matrix(c(0, 0, 1, 0,
    1, 1, 0, 0,
    1, 0, 0, 1), ncol = 4, byrow = TRUE)


  m <- mak(A, B)
  p <- plk(A)

  time <- seq(0, 1, by = 0.01)
  x0 <- c(X = 1, Y = 2, Z = 0, W = 1)

  y <- numsolve(o = m, time = time, x0 = x0, param = 1:3)
  x <- numsolve(o = m, time = time, x0 = x0 + .01, param = 1:3)


  w <- abs(matrix(rnorm(prod(dim(y) - c(0, 1))), nrow = nrow(y)))

  m$rs[[1]]$fixed <- FALSE
  # # Used test_loss to asses a few things
  # foo <- test_loss(ode_struct = m, opt_struct = opt(y, weights = w),
  #   x0 = matrix(x0), param_ = list(k = 1:3 - .5), sc = list(k = NULL), x = x)
  # head(foo, 14)
  # foo$l2_dx0
  #
  # f <- numsolve(o = m, time = time, x0 = x0, param = 1:3 - .5, approx_sensitivity = TRUE)
  # expect_equal(foo$l, sum(w * ((y - f$trajectory)[, -1])^2) / (2 * (101 - 1)))
  #
  # expect_equal(apply(f$sensitivity_param[[1]], 3, function(sl) - sum(sl * w * (y - f$trajectory)[, -1]) / (101-1) ),
  #   as.vector(foo$l3_dparam))
  #
  # expect_equal(apply(f$sensitivity_x0, 3, function(sl) - sum(sl * w * (y - f$trajectory)[, -1]) / (101-1) ),
  #   as.vector(foo$l4_dx0))
  #
  #
  # # integral matching
  # expect_equal(foo$i_dx0, matrix(0, nrow = 4))
  # expect_equal(foo$i2_dx0, matrix(0, nrow = 4))
  # expect_equal(foo$i3_dx0, matrix(0, nrow = 4))
  # expect_equal(foo$i4_dx0, matrix(0, nrow = 4))
  # expect_equal(foo$i5_dx0, matrix(0, nrow = 4))
  # foo[-length(foo)]
  #
  # expect_equal(as.vector(t(x[-1, -1] - x[-nrow(y), -1])),
  #   as.vector(foo$XYW$Y))
  # expect_equal(as.vector(t(w[-1,] + w[-nrow(w),]) / 2),
  #   as.vector(foo$XYW$W))
  #
  # expect_equal(-t(foo$XYW$X[[1]]) %*% ((foo$XYW$Y - foo$XYW$X[[1]] %*% (1:3 - .5)) * foo$XYW$W) / (101 - 1),
  #   foo$i3_dparam)
  #
  # expect_equal(-t(foo$XYW$X[[1]]) %*% ((foo$XYW$Y - foo$XYW$X[[1]] %*% (1:3)) * foo$XYW$W) / (101 - 1),
  #   foo$i5_dparam)
  #
  #
  #
  # # Exact gradients
  # m$rs[[1]]$exact_gradient <- TRUE
  # foo <- test_loss(ode_struct = m, opt_struct = opt(y),
  #   x0 = matrix(x0), param_ = list(k = 1:3 - .5), sc = list(k = NULL), x = x)
  #
  # f <- numsolve(o = m, time = time, x0 = x0, param = 1:3 - .5, approx_sensitivity = FALSE)
  # expect_equal(foo$l, sum(((y - f$trajectory)[, -1])^2) / (2 * (101 - 1)))
  #
  # expect_equal(apply(f$sensitivity_x0, 3, function(sl) - sum(sl * (y - f$trajectory)[, -1]) / (101-1) ),
  #   as.vector(foo$l4_dx0))
  #
  #
  # m$rs[[1]]$exact_gradient <- FALSE
  # m$rs[[2]]$exact_gradient <- TRUE
  # foo <- test_loss(ode_struct = m, opt_struct = opt(y),
  #   x0 = matrix(x0), param_ = list(k = 1:3 - .5), sc = list(k = NULL), x = x)
  #
  # f <- numsolve(o = m, time = time, x0 = x0, param = 1:3 - .5, approx_sensitivity = FALSE)
  # expect_equal(foo$l, sum(((y - f$trajectory)[, -1])^2) / (2 * (101 - 1)))
  #
  # expect_equal(apply(f$sensitivity_param[[1]], 3, function(sl) - sum(sl * (y - f$trajectory)[, -1]) / (101-1) ),
  #   as.vector(foo$l3_dparam))


  # PLK
  time <- seq(0, 1, by = 0.01)
  x0 <- c(X = 1, Y = 2, Z = 0, W = 1)

  expect_equal(numsolve(o = p, time = time, x0 = x0, param = t((B - A) * 1:3)),
    numsolve(o = m, time = time, x0 = x0, param = 1:3))
  y <- numsolve(o = p, time = time, x0 = x0, param = t((B - A) * 1:3))

  # foo <- test_loss(ode_struct = p, opt_struct = opt(y),
  #   x0 = matrix(x0), param_ = list(k = t((B - A) * 3:1) - .5), sc = list(k = NULL), x = x)
  #
  # f <- numsolve(o = p, time = time, x0 = x0, param = t((B - A) * 3:1) - .5, approx_sensitivity = TRUE)
  # expect_equal(foo$l, sum(((y - f$trajectory)[, -1])^2) / (2 * (101 - 1)))
  #
  # expect_equal(apply(f$sensitivity_param[[1]], 3, function(sl) - sum(sl * (y - f$trajectory)[, -1]) / (101-1) ),
  #   as.vector(foo$l3_dparam))
  #
  # expect_equal(apply(f$sensitivity_x0, 3, function(sl) - sum(sl * (y - f$trajectory)[, -1]) / (101-1) ),
  #   as.vector(foo$l4_dx0))

})


test_that("ratmak", {
  A <- rbind(diag(4), c(1, 1, 0, 1))
  CC <- matrix(c(0, 0, 1, 0,
    1, 1, 0, 0,
    1, 0, 0, 1), ncol = 4, byrow = TRUE) -
    matrix(c(1, 1, 0, 0,
      0, 0, 1, 0,
      0, 0, 1, 0), ncol = 4, byrow = TRUE)

  rat <- ratmak(A, CC)
  expect_equal(rat$rs[[2]]$nrow, nrow(CC))
  expect_equal(rat$rs[[2]]$p, nrow(A) * nrow(CC))
  expect_equal(rat$rs[[3]]$nrow, nrow(CC))
  expect_equal(rat$rs[[3]]$p, nrow(A) * nrow(CC))

  rat <- ratmak(A, CC, r1 = reg("scad", lower = -1))
  expect_equal(rat$rs[[2]]$nrow, nrow(CC))
  expect_equal(rat$rs[[2]]$p, nrow(A) * nrow(CC))
  expect_equal(rat$rs[[3]]$nrow, nrow(CC))
  expect_equal(rat$rs[[3]]$p, nrow(A) * nrow(CC))


  time <- seq(0, 1, by = 0.01)
  x0 <- c(X = 1, Y = 2, Z = 0, W = 1)

  z <- numsolve(o = rat, time = time, x0 = x0,
    param = list(theta1 = rep(1, rat$rs[[2]]$p), theta2 = rep(1, rat$rs[[3]]$p)))
  expect_equal(attr(z, "conv_code"), 0)

  z <- numsolve(o = rat, time = c(time, time), x0 = cbind(x0, x0 + .1),
    param = list(theta1 = cbind(rep(1, rat$rs[[2]]$p), .5), theta2 = cbind(rep(1, rat$rs[[3]]$p), .5)),
    approx_sensitivity = TRUE)
  z$sensitivity_param
  z$sensitivity_x0


  # compare to mak
  A <- matrix(c(1, 1, 0, 0,
    0, 0, 1, 0,
    0, 0, 1, 0), ncol = 4, byrow = TRUE)
  CC <- matrix(c(0, 0, 1, 0,
    1, 1, 0, 0,
    1, 0, 0, 1), ncol = 4, byrow = TRUE) - A
  rat <- ratmak(A, CC)

  time <- seq(0, 1, by = 0.01)
  x0 <- c(X = 1, Y = 2, Z = 0, W = 1)

  z <- numsolve(o = rat, time = time, x0 = x0,
    param = list(theta1 = as.vector(diag(1:3)), theta2 = rep(0, rat$rs[[3]]$p)))
  expect_equal(numsolve(o = mak(A, CC + A), time = time, x0 = x0, param = 1:3), z)

})
