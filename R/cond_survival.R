#' Conditional survival matrix and O_tau
#'
#' TAUS (Target-Age Unified Survival): conditional survival P(T>tau|T>t) and
#' the probability O_tau that a randomly sampled individual of unknown age
#' will outlive tau.
#'
#' @name cond_survival
#' @keywords internal
NULL

#' Conditional survival matrix and O_tau
#'
#' Computes the conditional survival matrix P(T>tau|T>t) for all t, tau and
#' the marginal probability O_tau that a randomly sampled individual of
#' unknown age will outlive tau, with confidence intervals.
#'
#' @param data A life table (data frame) with time and event columns.
#' @param var Character; name of the grouping variable (categorical).
#' @param cats Character vector or NULL; additional categorical variables for grouping.
#' @param time_var Character; name of the time variable.
#' @param event_var Character; name of the event variable (0/1 or TRUE/FALSE).
#' @param res Numeric; resolution (step size) for the grid of t and tau.
#' @param conf_int_level Numeric; confidence level (default 0.95).
#' @param aggregate_by_cats Logical; if TRUE, separate by each combination of var and cats.
#' @return A list with \code{cond_surv} (tibble of conditional survival and O_tau) and \code{args} (call arguments).
#' @importFrom magrittr "%>%"
#' @importFrom dplyr mutate group_by ungroup distinct all_of
#' @importFrom tidyr unnest complete nesting unite
#' @importFrom purrr map map_dfr
#' @importFrom rlang sym syms
#' @importFrom tibble tibble
#' @importFrom survival Surv survfit
#' @importFrom stats approx
#' @export
cond_surv_mat <- function(data,
                          var,
                          cats = NULL,
                          time_var,
                          event_var,
                          res,
                          conf_int_level = 0.95,
                          aggregate_by_cats = FALSE) {
  if (is.null(cats)) cats <- character(0)
  # unique_label in the original data to track grouping variable
  if (aggregate_by_cats && length(cats) > 0) {
    data <- data %>%
      tidyr::unite("unique_label", dplyr::all_of(c(var, cats)), sep = "_", remove = FALSE)
  } else {
    data <- data %>%
      dplyr::mutate(unique_label = paste0(!!rlang::sym(var)))
  }

  cond_surv <- data %>%
    dplyr::distinct(dplyr::across(dplyr::all_of(c(var, cats))))

  if (aggregate_by_cats && length(cats) > 0) {
    cond_surv <- cond_surv %>%
      tidyr::unite("unique_label", dplyr::all_of(c(var, cats)), sep = "_", remove = FALSE)
  } else {
    cond_surv <- cond_surv %>%
      dplyr::mutate(unique_label = paste0(!!rlang::sym(var)))
  }

  taus <- seq(1, max(data[[time_var]], na.rm = TRUE), res)

  cond_surv <- cond_surv %>%
    dplyr::mutate(
      conditional_survival = purrr::map(!!rlang::sym("unique_label"), ~ compute_conditional_survival(
        data = data,
        column = "unique_label",
        value = .x,
        time_var = time_var,
        event_var = event_var,
        taus = taus,
        conf_int_level = conf_int_level
      ))
    ) %>%
    tidyr::unnest(conditional_survival) # %>%
  # tidyr::complete(tidyr::nesting(!!rlang::sym(var), !!!rlang::syms(cats)), time = taus, tau = taus)

  if (aggregate_by_cats && length(cats) > 0) {
    cond_surv <- cond_surv %>%
      tidyr::unite("unique_label", dplyr::all_of(c(var, cats)), sep = "_", remove = FALSE)
  } else {
    cond_surv <- cond_surv %>%
      dplyr::mutate(unique_label = paste0(!!rlang::sym(var)))
  }

  cond_surv <- cond_surv %>%
    dplyr::group_by(!!rlang::sym(var), tau, !!!rlang::syms(cats)) %>%
    dplyr::mutate(O_tau = sum(Stau_St * age_dist, na.rm = TRUE)) %>%
    dplyr::mutate(O_tau_lo = sum(Stau_St_lo * age_dist, na.rm = TRUE)) %>%
    dplyr::mutate(O_tau_up = sum(Stau_St_up * age_dist, na.rm = TRUE)) %>%
    dplyr::ungroup()

  list(
    cond_surv = cond_surv,
    args = as.list(environment())
  )
}

#' Compute conditional survival for one group
#'
#' P(T>tau|T>t) for a specific group. Used by \code{cond_surv_mat}.
#'
#' @param data Data frame (life table).
#' @param column Character; grouping column name.
#' @param value Value of the group.
#' @param time_var,event_var Character; time and event column names.
#' @param taus Numeric vector; tau values.
#' @param conf_int_level Numeric; confidence level.
#' @return A tibble of conditional survival for that group.
#' @keywords internal
compute_conditional_survival <- function(data,
                                         column,
                                         value,
                                         time_var,
                                         event_var,
                                         taus,
                                         conf_int_level = 0.95) {
  data_subset <- data[data[[column]] == value, ]

  ci <- 1 - 2 * sqrt((1 - conf_int_level) / 2)

  kmfit <- survival::survfit(
    survival::Surv(data_subset[[time_var]], data_subset[[event_var]]) ~ 1,
    data = data_subset,
    conf.int = ci
  )

  max_time <- max(data_subset[[time_var]], na.rm = TRUE)
  time_points <- taus[taus <= max_time]

  S_t <- approx(x = kmfit$time, y = kmfit$surv, xout = time_points, yleft = 1)$y
  delta_t <- c(diff(time_points), 0)
  mu <- sum(S_t * delta_t)
  A_t <- (S_t * delta_t) / mu

  results <- purrr::map_dfr(time_points, function(tau) {
    S_tau <- summary(kmfit, times = tau)$surv
    S_tau_lo <- summary(kmfit, times = tau)$lower
    S_tau_up <- summary(kmfit, times = tau)$upper

    S_t <- summary(kmfit, times = time_points)$surv
    S_t_lo <- summary(kmfit, times = time_points)$lower
    S_t_up <- summary(kmfit, times = time_points)$upper

    Stau_St <- rep(NA, length(time_points))
    Stau_St_lo <- rep(NA, length(time_points))
    Stau_St_up <- rep(NA, length(time_points))

    valid_idx <- which(S_t > 0 & tau <= max_time)

    Stau_St[valid_idx] <- ifelse(time_points[valid_idx] >= tau, 1, S_tau / S_t[valid_idx])
    Stau_St_lo[valid_idx] <- ifelse(time_points[valid_idx] >= tau, 1, S_tau_lo / S_t_up[valid_idx])
    Stau_St_up[valid_idx] <- ifelse(time_points[valid_idx] >= tau, 1, S_tau_up / S_t_up[valid_idx])

    sum_St <- sum(S_t * delta_t)
    sum_S_tau_end <- sum(S_t[time_points >= tau] * delta_t[time_points >= tau])

    tibble::tibble(
      time = time_points,
      age_dist = A_t,
      tau = tau,
      Stau = S_tau,
      sum_St = sum_St,
      sum_S_tau_end = sum_S_tau_end,
      Stau_St = Stau_St,
      Stau_St_lo = Stau_St_lo,
      Stau_St_up = Stau_St_up
    )
  })

  results
}
