# Minimal survival data: two groups, same time/event columns
# Group A: 3 events at 2,4,6; 1 censored at 5
# Group B: 2 events at 3,5; 2 censored at 4,6
.make_surv_data <- function() {
  tibble::tibble(
    time   = c(2, 4, 6, 5,  3, 5, 4, 6),
    event  = c(1, 1, 1, 0,  1, 1, 0, 0),
    group  = c(rep("A", 4), rep("B", 4))
  )
}

test_that("cond_surv_mat returns list with cond_surv and args", {
  dat <- .make_surv_data()
  out <- cond_surv_mat(dat, var = "group", time_var = "time", event_var = "event", res = 2)
  expect_type(out, "list")
  expect_named(out, c("cond_surv", "args"))
  expect_s3_class(out$cond_surv, "tbl_df")
  expect_type(out$args, "list")
})

test_that("cond_surv_mat output has expected columns", {
  dat <- .make_surv_data()
  out <- cond_surv_mat(dat, var = "group", time_var = "time", event_var = "event", res = 2)
  cs <- out$cond_surv
  expect_true("time" %in% names(cs))
  expect_true("tau" %in% names(cs))
  expect_true("Stau_St" %in% names(cs))
  expect_true("O_tau" %in% names(cs))
  expect_true("O_tau_lo" %in% names(cs))
  expect_true("O_tau_up" %in% names(cs))
  expect_true("unique_label" %in% names(cs))
})

test_that("O_tau and CI bounds are in [0, 1]", {
  dat <- .make_surv_data()
  out <- cond_surv_mat(dat, var = "group", time_var = "time", event_var = "event", res = 2)
  cs <- out$cond_surv
  # One row per (group, time, tau) combination; O_tau is repeated, take unique per group-tau
  o_tau_df <- dplyr::distinct(cs, unique_label, tau, O_tau, O_tau_lo, O_tau_up)
  expect_true(all(o_tau_df$O_tau >= 0 & o_tau_df$O_tau <= 1, na.rm = TRUE))
  expect_true(all(o_tau_df$O_tau_lo >= 0 & o_tau_df$O_tau_lo <= 1, na.rm = TRUE))
  expect_true(all(o_tau_df$O_tau_up >= 0 & o_tau_df$O_tau_up <= 1, na.rm = TRUE))
})

test_that("cond_surv_mat conditional survival Stau_St is in [0, 1]", {
  dat <- .make_surv_data()
  out <- cond_surv_mat(dat, var = "group", time_var = "time", event_var = "event", res = 2)
  cs <- out$cond_surv
  valid <- cs$Stau_St[!is.na(cs$Stau_St)]
  expect_true(all(valid >= 0 & valid <= 1))
})

test_that("cond_surv_mat args reflect input", {
  dat <- .make_surv_data()
  out <- cond_surv_mat(dat, var = "group", time_var = "time", event_var = "event", res = 2)
  expect_equal(out$args$time_var, "time")
  expect_equal(out$args$event_var, "event")
  expect_equal(out$args$var, "group")
  expect_equal(out$args$conf_int_level, 0.95)
})
