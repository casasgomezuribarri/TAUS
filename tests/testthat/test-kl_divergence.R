test_that("kl_divergence(p, p) is zero", {
  p <- c(0.2, 0.3, 0.5)
  expect_equal(kl_divergence(p, p), 0, tolerance = 1e-9)
})

test_that("kl_divergence is non-negative", {
  tru <- c(0.25, 0.25, 0.5)
  est <- c(0.2, 0.3, 0.5)
  expect_true(kl_divergence(tru, est) >= 0)
})

test_that("kl_divergence with epsilon avoids log(0)", {
  tru <- c(1, 0, 0)
  est <- c(0.5, 0.5, 0)
  expect_true(is.finite(kl_divergence(tru, est)))
  expect_true(kl_divergence(tru, est) >= 0)
})

test_that("kl_divergence longer vectors", {
  tru <- rep(1/10, 10)
  est <- tru + runif(10, -0.05, 0.05)
  est <- est / sum(est)
  expect_true(kl_divergence(tru, est) >= 0)
  expect_true(is.finite(kl_divergence(tru, est)))
})
