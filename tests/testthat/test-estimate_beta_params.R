test_that("estimate_beta_params returns list with alpha and beta", {
  out <- estimate_beta_params(0.5, 0.3, 0.7)
  expect_type(out, "list")
  expect_named(out, c("alpha", "beta"))
  expect_true(out$alpha > 0)
  expect_true(out$beta > 0)
})

test_that("estimate_beta_params mean 0.5 symmetric CI gives positive alpha, beta", {
  out <- estimate_beta_params(0.5, 0.4, 0.6)
  expect_true(out$alpha > 0)
  expect_true(out$beta > 0)
  expect_equal(out$alpha, out$beta, tolerance = 0.01)
})

test_that("estimate_beta_params edge mean 1 (narrow CI) returns finite alpha, beta", {
  # var must be < 1e-3 to hit the mean=1 edge case in the implementation
  out <- estimate_beta_params(1, 0.9999, 1)
  expect_true(is.finite(out$alpha))
  expect_true(is.finite(out$beta))
  expect_true(out$alpha > 0)
  expect_true(out$beta >= 0)
})

test_that("estimate_beta_params edge mean 0 (narrow CI) returns finite alpha, beta", {
  # var must be < 1e-3 to hit the mean=0 edge case
  out <- estimate_beta_params(0, 0, 0.0005)
  expect_true(is.finite(out$alpha))
  expect_true(is.finite(out$beta))
  expect_true(out$beta > 0)
})

test_that("estimate_beta_params Beta(alpha, beta) has correct mean", {
  out <- estimate_beta_params(0.6, 0.5, 0.7)
  # E[Beta(alpha, beta)] = alpha / (alpha + beta)
  mean_beta <- out$alpha / (out$alpha + out$beta)
  expect_equal(mean_beta, 0.6, tolerance = 0.01)
})
