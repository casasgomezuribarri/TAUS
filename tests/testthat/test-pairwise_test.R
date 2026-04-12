# Reuse same minimal data and cond_surv_mat output
.make_surv_data <- function() {
  tibble::tibble(
    time   = c(2, 4, 6, 5,  3, 5, 4, 6),
    event  = c(1, 1, 1, 0,  1, 1, 0, 0),
    group  = c(rep("A", 4), rep("B", 4))
  )
}

test_that("pairwise_test returns data frame with expected columns", {
  dat <- .make_surv_data()
  cs <- cond_surv_mat(dat, var = "group", time_var = "time", event_var = "event", res = 1)
  out <- suppressMessages(pairwise_test(cs, var_values = c("A", "B"), tau_values = c(4, 4), B = 50))
  expect_s3_class(out, "data.frame")
  expect_true("group_1" %in% names(out))
  expect_true("group_2" %in% names(out))
  expect_true("p_value" %in% names(out))
  expect_true("ks_stat" %in% names(out))
  expect_true("O_tau_1" %in% names(out))
  expect_true("O_tau_2" %in% names(out))
})

test_that("pairwise_test p_values are in [0, 1]", {
  dat <- .make_surv_data()
  cs <- cond_surv_mat(dat, var = "group", time_var = "time", event_var = "event", res = 1)
  out <- suppressMessages(pairwise_test(cs, var_values = c("A", "B"), tau_values = c(4, 4), B = 50))
  expect_true(all(out$p_value >= 0 & out$p_value <= 1))
})

test_that("pairwise_test single tau for all groups", {
  dat <- .make_surv_data()
  cs <- cond_surv_mat(dat, var = "group", time_var = "time", event_var = "event", res = 1)
  out <- suppressMessages(pairwise_test(cs, var_values = NULL, tau_values = 4, B = 50))
  expect_gte(nrow(out), 1)
  expect_true(all(out$p_value >= 0 & out$p_value <= 1))
})

test_that("pairwise_test same group gives p_value near 1 when B is large", {
  # Two "groups" with same data -> same O_tau -> null should not be rejected
  dat <- .make_surv_data()
  dat2 <- dplyr::bind_rows(dat %>% dplyr::mutate(group = "A"), dat %>% dplyr::mutate(group = "A2"))
  cs <- cond_surv_mat(dat2, var = "group", time_var = "time", event_var = "event", res = 1)
  out <- suppressMessages(pairwise_test(cs, var_values = c("A", "A2"), tau_values = c(4, 4), B = 200))
  # Same data so O_tau should be very similar; p-value should be high (not reject null)
  expect_true(out$p_value >= 0.1)
})
