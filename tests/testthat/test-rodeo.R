###---------------###
###  Test rodeo   ###
###---------------###

library(episode)
context("rodeo object")

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


test_that("bull", {
  ro <- opt(y_safe1, lambda = exp(seq(log(1), log(0.0001), length.out = 24)))

  x0s <- array(x0, dim = c(length(x0), 1, 1))
  foo <- bull(m, ro, sc_ = list(k = NULL), params_ = list(list(k = k + .25)), x0s, FALSE)

  d <- foo$dfs
  expect_true(all(d >= 0 & d <= nrow(A)))
  expect_equal(colSums(( apply(foo$params[[1]][[2]], 2, numsolve,
    o = m, time = times_1,
    x0 = foo$params[[1]][[1]][, 1]) - as.vector(y_safe1))^2) / (2 * (nrow(y_safe1) - 1)),
    as.vector(foo$losses))


  m$rs$rate$a <- 3
  m$rs$rate$reg_type <- "scad"
  x0s <- array(x0, dim = c(length(x0), 2, 1))
  ro <- opt(y_safe2, lambda = exp(seq(log(1), log(0.0001), length.out = 24)))
  foo <- bull(m, ro, sc_ = list(k = context), params_ = list(list(k = k + .25)), x0s, FALSE)

  d <- foo$dfs
  expect_true(all(d >= 0 & d <= nrow(A)))


  # three initialisations
  Matrix::Matrix()
  x0s <- array(x0 + 1, dim = c(length(x0), 2, 4))
  m$rs[[1]]$fixed <- TRUE
  ro <- opt(y_safe2, lambda = exp(seq(log(1), log(0.0001), length.out = 24)))
  rodeo_ret <- bull(m, ro, sc_ = list(k = context), params_ = list(list(k = k + .25),
    list(k = rep(0, length(k))),
    list(k = rep(.1, length(k))),
    list(k = rep(.1, length(k)))), x0s, FALSE)
  d <- rodeo_ret$dfs

  # print(rodeo_ret$params[[3]])
  # print(rodeo_ret$params[[4]])
  # print(identical(rodeo_ret$params[[3]], rodeo_ret$params[[4]]))
  # expect_identical(rodeo_ret$params[[3]], rodeo_ret$params[[4]])
  expect_true(all(d >= 0 & d <= nrow(A)))
})



test_that("rodeo.ode", {
  w <- matrix(runif(prod(dim(y_safe2) - c(0, 1)), min = .9, max = 1.1), nrow = nrow(y_safe2))
  v <- runif(length(k), min = .9, max = 1.1)

  ro <- opt(y_safe2, weights = w, nlambda = 14, lambda_min_ratio = 0.0012, tol = 1e-6)
  m <- mak(A, B,
    r = reg(reg_type = "l1", step_max = 200, a = 3.25,
      upper = seq_along(k) / 1.1, backtrack_max = 25,
      lambda_factor = 1, exclude = c(6),
      screen = TRUE, penalty_factor = v,
      contexts = cbind(1, rep(c(1, 2), length.out = 7)),
      scales = c(10.1, 1.2, 1.2, 1.2, rep(1.2, 2), 1.2),
      exact_gradient = TRUE),
    rx0 = reg(reg_type = "none", fixed = TRUE))

  rod <- rodeo(m, ro, x0 = NULL, param = list(k = rep(0, length(k))))

  if (rod$o$rs[[2]]$reg_type == "l1") {
    expect_equal(
      rod$op$lambda * m$rs[[2]]$lambda_factor * as.vector(m$rs[[2]]$penalty_factor %*% (abs(rod$params$k) / m$rs[[2]]$scales)),
      rod$penalties[2, ])
  }

  expect_true(all(rod$x0s == as.vector(t(y_safe2[c(TRUE, diff(y_safe2[, 1]) < 0), -1]))))

  expect_true(all(rod$dfs >= 0 & rod$dfs <= nrow(A)))


  ## mukltiple initialisations
  m <- mak(A, B,
    r = reg(reg_type = "scad", step_max = 200, a = 3.25,
      upper = seq_along(k) / 1.1, backtrack_max = 25,
      lambda_factor = .2, exclude = c(6, 7),
      screen = TRUE, penalty_factor = v,
      contexts = cbind(1, rep(c(1, 2), length.out = 7)),
      scales = c(10.1, 1.2, 1.2, 1.2, rep(1.2, 2), 1.2),
      exact_gradient = TRUE),
    rx0 = reg(reg_type = "none", fixed = TRUE))

  rod <- rodeo(m, ro, x0 = cbind(rep(x0, 2), rep(x0, 2) - .1), param = list(k = cbind(rep(0, length(k)), 1)))



  ## plk
  w <- matrix(runif(prod(dim(y_safe2) - c(0, 1)), min = .9, max = 1.1), nrow = nrow(y_safe2))
  v <- runif(length(A), min = .9, max = 1.1)

  p <- plk(A,
    r = reg(reg_type = "scad", step_max = 200, a = 3.25,
      upper = seq_along(A) / 1.1, backtrack_max = 25,
      lambda_factor = .8,
      screen = TRUE, penalty_factor = v,
      contexts = cbind(1, rep(c(1, 2), length.out = prod(dim(A)))),
      #scales = rep(10, 7),
      scales = c(10.1, rep(1.2, length(A) - 1)),
      exact_gradient = TRUE),
    rx0 = reg(reg_type = "none", fixed = TRUE))
  ro <- opt(y_safe2, weights = w, nlambda = 14, lambda_min_ratio = 0.0012, tol = 1e-6)
  rod <- rodeo(p, ro, x0 = cbind(rep(x0, 2), rep(x0, 2) - .1),
    param = list(k = cbind(rep(0, length(A)), c(1, -1))))


  ## rlk
  w <- matrix(runif(prod(dim(y_safe2) - c(0, 1)), min = .9, max = 1.1), nrow = nrow(y_safe2))
  v <- runif(length(A), min = .9, max = 1.1)

  rr <- rlk(A, B,
           r = reg(reg_type = "scad", step_max = 200, a = 3.25,
                   upper = seq_along(A) / 1.1, backtrack_max = 25,
                   lambda_factor = 1,
                   screen = TRUE, penalty_factor = v,
                   contexts = cbind(1, rep(c(1, 2), length.out = prod(dim(A)))),
                   #scales = rep(10, 7),
                   scales = c(10.1, rep(1.2, length(A) - 1)),
                   exact_gradient = TRUE),
           rx0 = reg(reg_type = "none", fixed = TRUE))
  ro <- opt(y_safe2, weights = w, nlambda = 14, lambda_min_ratio = 0.0012, tol = 1e-6)
  rod <- rodeo(rr, ro, x0 = cbind(rep(x0, 2), rep(x0, 2) - .1),
    param = list(k = cbind(rep(0, length(A)), c(1, -1))))



  ## ratmak
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
  ro <- opt(y_safe2, weights = w, nlambda = 44, lambda_min_ratio = 0.000012, tol = 1e-6)


  rod <- rodeo(rat, ro, x0 = rep(x0, 2),
    param = list(theta1 = rep(0, length(v)),
      theta2 = rep(0, length(v))))

  rod$params
  rod$steps
  rod$dfs
  rod$codes
  rod$jerr

})


