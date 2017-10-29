###----------------###
###  Test mapping  ###
###----------------###

library(episode)
context("mapping")




test_that("power_integral", {
  A <- matrix(c(1, 1, 0, 0,
    0, 1, 1, 0,
    0, 0, 1, 1,
    2, 1, 0, 0,
    2, 0, 0, 1,
    0, 2, 2, 1,
    0, 1, 1, 2), ncol = 4, byrow = TRUE)
  B <- A[c(4, 3, 2, 1, 5, 7, 6), c(1, 3, 2, 4)]
  m <- mak(A, B)
  k <- c(k1 = 1, k2 = 2, k3 = .5, k4 = 0, k5 = 1, k6 = 0, k7 = .1)
  x0 <- c(X = 1, Y = 2, Z = .1, W = 2)

  ti <- seq(0, 1, by = .01)
  y <- numsolve(o = m, time = ti, x0 = x0, param = list(rate = k))
  x <- numsolve(o = m, time = seq(0, 1, by = .001), x0 = x0, param = list(rate = k))
  z <- numsolve(o = m, time = seq(0, 1, by = .1), x0 = x0, param = list(rate = k))

  # if x = y, then easy to check the first row is correct
  ss <- numint(time = ti, x = y, type = "power", A = A, B = A)
  expect_equal(as.vector((y[2, 1] - y[1, 1]) * (exp(A %*% log(y[1, -1])) + exp(A %*% log(y[2, -1])) + 4 * exp(A %*% log((y[1, -1] + y[2, -1]) / 2))) / 6),
    ss[1, ])
  expect_equal(as.vector((y[3, 1] - y[2, 1]) * (exp(A %*% log(y[2, -1])) + exp(A %*% log(y[3, -1])) + 4 * exp(A %*% log((y[2, -1] + y[3, -1]) / 2))) / 6),
    ss[2, ])
})

# ss <- numint(time = ti, x = y, type = "fracpower", A = A, B = A[1, , drop = FALSE])

