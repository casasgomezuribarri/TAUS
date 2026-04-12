#' Beta parameters and pairwise comparison of O_tau
#'
#' @name stats
NULL

#' Estimate Beta shape parameters from mean and CI limits
#'
#' Used to model O_tau as Beta(alpha, beta) from point estimate and CI width.
#'
#' @param mean Numeric; point estimate (e.g. O_tau).
#' @param lower Numeric; lower CI limit.
#' @param upper Numeric; upper CI limit.
#' @return A list with \code{alpha} and \code{beta}.
#' @export
estimate_beta_params <- function(mean, lower, upper) {
  if (any(is.na(c(mean, lower, upper)))) {
    return(list(alpha = NA_real_, beta = NA_real_))
  }
  var <- upper - lower
  if (var < 1e-3 && abs(mean - 1) < 1e-8) {
    alpha <- 1e6
    beta <- 1e-6
  } else if (var < 1e-3 && abs(mean) < 1e-8) {
    alpha <- 1e-6
    beta <- 1e6
  } else {
    alpha <- ((1 - mean) / var^2 - 1 / mean) * mean^2
    beta <- alpha * (1 / mean - 1)
  }
  list(alpha = alpha, beta = beta)
}

#' Pairwise comparison of O_tau across groups
#'
#' Fits Beta distributions to each O_tau, computes KS distance, and bootstrap p-value
#' under the null that both come from the same Beta.
#'
#' @param cond_surv Output of \code{\link{cond_surv_mat}}.
#' @param var_values Character vector or NULL; groups to compare (NULL = all at same tau).
#' @param tau_values Numeric; if \code{var_values} is NULL, a single tau; else vector of tau per group.
#' @param seed Integer; random seed for bootstrap.
#' @param B Integer; number of bootstrap replicates.
#' @param n_draw Integer; draws per bootstrap replicate for KS.
#' @return A data frame of pairwise comparisons (group_1, group_2, tau, O_tau, ks_stat, p_value, etc.), invisibly.
#' @importFrom magrittr "%>%"
#' @importFrom dplyr filter group_by slice mutate semi_join
#' @importFrom tibble tibble
#' @importFrom stats ks.test optimize pbeta rbeta
#' @export
pairwise_test <- function(cond_surv,
                          var_values = NULL,
                          tau_values = NULL,
                          seed = 314,
                          B = 5000,
                          n_draw = 10) {
  set.seed(seed)

  analytic_ks <- function(a1, b1, a2, b2) {
    optimize(
      function(x) abs(pbeta(x, a1, b1) - pbeta(x, a2, b2)),
      interval = c(0, 1),
      maximum = TRUE
    )$objective
  }

  bootstrap_p <- function(ks_deterministic, a1, b1, a2, b2) {
    a0 <- mean(c(a1, a2))
    b0 <- mean(c(b1, b2))
    null_stats <- replicate(B, {
      x <- rbeta(n_draw, a0, b0)
      y <- rbeta(n_draw, a0, b0)
      suppressWarnings(ks.test(x, y)$statistic)
    })
    mean(null_stats >= ks_deterministic)
  }

  if (is.null(var_values)) {
    if (length(tau_values) != 1) stop("if var_values is not specified, tau_values must be a single number.")
    filtered_cond_surv <- cond_surv$args$cond_surv %>%
      dplyr::filter(tau == tau_values) %>%
      dplyr::group_by(unique_label) %>%
      dplyr::slice(1) %>%
      dplyr::mutate(
        alpha = estimate_beta_params(O_tau, O_tau_lo, O_tau_up)$alpha,
        beta = estimate_beta_params(O_tau, O_tau_lo, O_tau_up)$beta
      )
  } else {
    if (length(tau_values) != length(var_values)) stop("tau_values and var_values must have the same length.")
    pairs <- tibble::tibble(tau = tau_values, unique_label = var_values)
    filtered_cond_surv <- cond_surv$cond_surv %>%
      dplyr::semi_join(pairs, by = c("tau", "unique_label")) %>%
      dplyr::group_by(unique_label) %>%
      dplyr::slice(1) %>%
      dplyr::mutate(
        alpha = estimate_beta_params(O_tau, O_tau_lo, O_tau_up)$alpha,
        beta = estimate_beta_params(O_tau, O_tau_lo, O_tau_up)$beta
      )
  }

  results <- data.frame(
    group_1 = character(),
    tau_1 = numeric(),
    O_tau_1 = numeric(),
    group_2 = character(),
    tau_2 = numeric(),
    O_tau_2 = numeric(),
    ks_stat = numeric(),
    effect_size = numeric(),
    effect_ratio = numeric(),
    p_value = numeric()
  )

  for (i in 1:(nrow(filtered_cond_surv) - 1)) {
    for (j in (i + 1):nrow(filtered_cond_surv)) {
      a1 <- filtered_cond_surv$alpha[i]
      b1 <- filtered_cond_surv$beta[i]
      a2 <- filtered_cond_surv$alpha[j]
      b2 <- filtered_cond_surv$beta[j]

      ks_stat <- analytic_ks(a1, b1, a2, b2)
      p_val <- bootstrap_p(ks_stat, a1, b1, a2, b2)

      results <- rbind(results, data.frame(
        group_1 = filtered_cond_surv$unique_label[i],
        tau_1 = filtered_cond_surv$tau[i],
        O_tau_1 = filtered_cond_surv$O_tau[i],
        group_2 = filtered_cond_surv$unique_label[j],
        tau_2 = filtered_cond_surv$tau[j],
        O_tau_2 = filtered_cond_surv$O_tau[j],
        ks_stat = ks_stat,
        effect_size = filtered_cond_surv$O_tau[i] - filtered_cond_surv$O_tau[j],
        effect_ratio = filtered_cond_surv$O_tau[i] / filtered_cond_surv$O_tau[j],
        p_value = p_val
      ))
    }
  }

  if (requireNamespace("knitr", quietly = TRUE)) {
    print(knitr::kable(results, digits = 3))
  } else {
    print(results)
  }
  invisible(results)
}

#' KL divergence between two distributions
#'
#' Used e.g. for power analysis (estimated vs true age distribution).
#'
#' @param tru,est Numeric vectors; true and estimated probabilities (same length).
#' @param epsilon Small constant to avoid log(0).
#' @return Scalar; KL divergence.
#' @export
kl_divergence <- function(tru, est, epsilon = 1e-10) {
  tru <- tru + epsilon
  est <- est + epsilon
  trh <- tru / sum(tru)
  est <- est / sum(est)
  sum(tru * log(tru / est))
}
