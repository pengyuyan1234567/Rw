test_that("RemoveSeasonal removes a known monthly cycle", {
  t_mmdd <- rep(1:12, each = 5)
  mu_true <- 10 * (1:12)
  noise <- rep(c(-2, -1, 0, 1, 2), times = 12)
  y <- mu_true[t_mmdd] + noise

  out <- RemoveSeasonal(y, t_mmdd)

  expect_equal(unname(out$mu), as.numeric(mu_true))
  expect_equal(names(out$mu), as.character(1:12))
  expect_equal(out$anomaly, noise)
})

test_that("RemoveSeasonal is invariant to the order of the observations", {
  t_mmdd <- c(101, 101, 102, 102, 103, 103)
  y <- c(1, 3, 10, 20, -5, 5)

  out <- RemoveSeasonal(y, t_mmdd)
  shuffled <- c(6, 4, 3, 1, 2, 5)
  out_shuffled <- RemoveSeasonal(y[shuffled], t_mmdd[shuffled])

  expect_equal(out$mu, out_shuffled$mu)
  expect_equal(out$anomaly[shuffled], out_shuffled$anomaly)
})

test_that("LS_OneBreak reproduces lm() with a step dummy", {
  X <- seq_len(20)
  Y <- c(seq_len(10), seq_len(10) + 5)
  Ic <- 10L

  out <- LS_OneBreak(Y, X, Ic)
  ref <- stats::lm(Y ~ X + I(X > Ic))

  expect_equal(unname(out$coef), unname(stats::coef(ref)))
  expect_equal(out$fitted, unname(stats::fitted(ref)))
  expect_equal(out$resi, unname(stats::residuals(ref)))
  expect_equal(out$SSE, sum(stats::residuals(ref)^2))
})

test_that("LS_OneBreak without a break reduces to a simple linear trend", {
  X <- seq_len(20)
  Y <- 3 + 2 * X + rep(c(-0.5, 0.5), 10)

  out <- LS_OneBreak(Y, X, Ic = NA)
  ref <- stats::lm(Y ~ X)

  expect_length(out$coef, 2L)
  expect_equal(unname(out$coef), unname(stats::coef(ref)))
  expect_equal(out$SSE, sum(stats::residuals(ref)^2))
})

test_that("LS_MultiBreaks reproduces lm() with one dummy per extra segment", {
  X <- seq_len(30)
  Y <- c(seq_len(10), seq_len(10) + 5, seq_len(10) + 10)
  Ips <- c(10L, 20L, 30L)

  out <- LS_MultiBreaks(Y, X, Ips)
  ref <- stats::lm(Y ~ I(X - mean(X)) + I(X > 10 & X <= 20) + I(X > 20))

  expect_length(out$coef, 4L)
  expect_equal(unname(out$coef), unname(stats::coef(ref)))
  expect_equal(out$fitted, unname(stats::fitted(ref)))
  expect_equal(out$SSE, sum(stats::residuals(ref)^2))
})

test_that("only.SSE returns the same SSE with less output", {
  X <- seq_len(30)
  Y <- c(seq_len(10), seq_len(10) + 5, seq_len(10) + 10)
  Ips <- c(10L, 20L, 30L)

  full <- LS_MultiBreaks(Y, X, Ips)
  fast <- LS_MultiBreaks(Y, X, Ips, only.SSE = TRUE)

  expect_named(fast, "SSE")
  expect_equal(fast$SSE, full$SSE)
  expect_equal(full$SSE, sum(full$resi^2))
})

test_that("the break point must be inside the series", {
  X <- seq_len(10)
  Y <- X

  expect_error(LS_MultiBreaks(Y, X, Ips = c(5, 11)))
  expect_error(LS_MultiBreaks(Y, X[-1], Ips = c(10)))
})

test_that("LS_MultiBreaks without a break reduces to a simple linear trend", {
  X <- seq_len(15)
  Y <- 1 + 0.5 * X + rep(c(-1, 0, 1), 5)

  out <- LS_MultiBreaks(Y, X, Ips = length(Y))
  ref <- stats::lm(Y ~ I(X - mean(X)))

  expect_length(out$coef, 2L)
  expect_equal(unname(out$coef), unname(stats::coef(ref)))
  expect_equal(out$SSE, sum(stats::residuals(ref)^2))
})

test_that("RemoveSeasonal drops NA when building the climatology", {
  t_mmdd <- c(1, 1, 1, 2, 2)
  y <- c(1, 2, NA, 10, 12)

  out <- RemoveSeasonal(y, t_mmdd)

  expect_equal(unname(out$mu), c(1.5, 11))
  expect_equal(out$anomaly, c(-0.5, 0.5, NA, -1, 1))
})
