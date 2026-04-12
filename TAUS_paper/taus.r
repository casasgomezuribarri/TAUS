# this script has the functions necessary to apply TAUS - Time-Agnostic Unified Survival

# This approach calculates the probability that a randomly sampled individual from a population
# will die after a certain age.

# Equivalently: the probability that a randomly sampled individual from a population will survive
# to (or have already survived past) a certain age.

# to do:
# - most functions need a lot of arguments just for the computation of cond_surv. instead, make cond_surv the argument.
# - check ps_heatmap function - a bit outdated
##################################################################################################################################
# environment set up
##################################################################################################################################

# load packages
packages <- c(
  "tidyverse",
  "tibble",
  "dplyr",
  "tidyr",
  "rstudioapi",
  "survival",
  "ggplot2",
  #  "ComplexHeatmap",
  "rlang",
  "circlize"
)

for (i in packages) {
  if (!require(i, character.only = TRUE)) install.packages(i)
  library(i, character.only = TRUE)
}

################################################################################################
# conditional survival
################################################################################################

# generates a matrix that is a dataset of conditional survival values - I'm calling it a conditional survival life table
# to do: control over unique_label would be nice - not a priority though

cond_surv_mat <- function(data, # the dataset - a typical life table
                          var, # variable of interest (string, must be categorical)
                          cats = NULL, # other categorical variables that matter (vector of strings, all categorical variables)
                          time_var, # the time variable (string)
                          event_var, # the event variable (string)
                          res, # resolution of the conditional survival grids for t and tau (step size of their ranges)
                          conf_int_level = 0.95, # confidence interval level (default 0.95, i.e. 95% CI)
                          aggregate_by_cats = FALSE # if TRUE, will separate by each unique combination of values in var and cats
) {
  # unique_label in the original data to track grouping variable
  if (aggregate_by_cats) {
    data <- data %>%
      tidyr::unite("unique_label", all_of(c(var, cats)), sep = "_", remove = FALSE)
  } else {
    data <- data %>%
      mutate(unique_label = paste0(!!sym(var)))
  }

  # create new dataframe, one row per unique combination of categorical variables
  cond_surv <- data %>%
    distinct(across(all_of(c(var, cats))))

  # unique_label in cond_surv too, also to track grouping variable
  if (aggregate_by_cats) {
    cond_surv <- cond_surv %>%
      tidyr::unite("unique_label", all_of(c(var, cats)), sep = "_", remove = FALSE)
  } else {
    cond_surv <- cond_surv %>%
      mutate(unique_label = paste0(!!sym(var)))
  }


  # define the space for taus
  taus <- seq(1, max(data[[time_var]], na.rm = TRUE), res) # all values of t, with a step size of res

  # compute conditional survivals
  cond_surv <- cond_surv %>%
    mutate(
      conditional_survival = map(!!sym("unique_label"), ~ compute_conditional_survival( # also a custom function (the next definition), most arguments are inherited
        data = data, # original life table data
        column = "unique_label", # variable to group by (this keeps track of aggregation)
        value = .x, # syntax of map() function (similar to a for loop, but more efficient) - .x is each value in var
        time_var = time_var, # time variable in data
        event_var = event_var, # event variable in data
        taus = taus # taus to compute conditional survival for
      ))
    ) %>%
    unnest(conditional_survival) # %>%
  # complete(nesting(!!sym(var), !!!syms(cats)), time = taus, tau = taus) # add missing rows to standardise plots

  # and expand unique_label values accordingly to avoid NAs
  if (aggregate_by_cats) {
    cond_surv <- cond_surv %>%
      tidyr::unite("unique_label", all_of(c(var, cats)), sep = "_", remove = FALSE)
  } else {
    cond_surv <- cond_surv %>%
      mutate(unique_label = paste0(!!sym(var)))
  }

  # add O_tau - the probability that a randomly sampled individual will outlive tau
  cond_surv <- cond_surv %>%
    group_by(!!sym(var), tau, !!!syms(cats)) %>% # for each unique combination of var, cats, and tau
    mutate(O_tau = sum(Stau_St * age_dist, na.rm = TRUE)) %>% # compute Σ(0-inf) [P(T>tau|T>t) * A(t)]
    mutate(O_tau_lo = sum(Stau_St_lo * age_dist, na.rm = TRUE)) %>% # propagate lower bound of 95%CI
    mutate(O_tau_up = sum(Stau_St_up * age_dist, na.rm = TRUE)) %>% # propagate upper bound of 95%CI
    ungroup() # go back to original format

  return(list(
    cond_surv = cond_surv,
    args = as.list(environment())
  ))
}

# compute conditional survival for a specific group - P(T>tau | T>t) ∀ tau ∈ taus, ∀ t ∈ [1, max_time]
compute_conditional_survival <- function(data, # the dataset
                                         column, # variable of interest (string)
                                         value, # the value within column to compute conditional survival for
                                         time_var, # the time variable (string)
                                         event_var, # the event variable (string)
                                         taus, # the taus to compute conditional survival for
                                         conf_int_level = 0.95 # confidence interval level (default 0.95, i.e. 95% CI)
) {
  data_subset <- data[data[[column]] == value, ] # subset data based on column value

  # since this metric will be a ratio, we need to make sure that the propagated total uncertainty corresponds to the 95%CI
  ci <- 1 - 2 * sqrt((1 - conf_int_level) / 2) # setting ci to this value ensures that the final ratio will have 95%CI

  # first fit km curve for our subset
  kmfit <- survfit(Surv(data_subset[[time_var]], data_subset[[event_var]]) ~ 1,
    data = data_subset,
    conf.int = ci
  )

  # get the maximum age in the subset - conditional survival won't be defined beyond this.
  max_time <- max(data_subset[[time_var]], na.rm = TRUE)

  # define relevant time range
  time_points <- taus[taus <= max_time] # all the possible taus under current condition.
  # A_t = approx(x = kmfit$time, y = kmfit$surv, xout = time_points, yleft = 1)$y # Age pyramid, A(t) - survival function, for ℕ ∈ [1, max_time]
  # A_t = A_t / sum(A_t) # scaled to sum 1 (PMF)

  # better way of doing the above (explicit consideration of bin width) - I'll keep the old version above for reference or in case this breaks anything downstream
  S_t <- approx(x = kmfit$time, y = kmfit$surv, xout = time_points, yleft = 1)$y # Age pyramid - survival function, for all ℕ ∈ [1, max_time]
  delta_t <- c(diff(time_points), 0) # compute interval widths Δt between successive ages - Δt tells us how long each survival value "lasts".
  mu <- sum(S_t * delta_t) # normalization constant μ = ∫ S(t) dt (this is how the integral should be approximated, not just with the sum!)
  A_t <- (S_t * delta_t) / mu # turn S(t) into a proper probability mass function for ages: A(t) = S(t) * Δt / μ


  # compute conditional survival - here's the 'hing!
  results <- map_dfr(time_points, function(tau) { # mapping > for loop (vectorised ops > rowwise ops). We'll do the following ∀ tau ∈ taus
    # tau = taus[10]

    # extract empirical P(T>tau), and its { 1-2*sqrt((1-conf_int_level)/2) } % CI limits from the km curve
    S_tau <- summary(kmfit, times = tau)$surv
    S_tau_lo <- summary(kmfit, times = tau)$lower
    S_tau_up <- summary(kmfit, times = tau)$upper

    # extract empirical P(T>t) ∀ t ∈ [1, max_time], and their { 1-2*sqrt((1-conf_int_level)/2) } % CI limits - also from the km curve
    S_t <- summary(kmfit, times = time_points)$surv
    S_t_lo <- summary(kmfit, times = time_points)$lower
    S_t_up <- summary(kmfit, times = time_points)$upper


    # initialise vectors to store results for this tau (default value is undefined)
    Stau_St <- rep(NA, length(time_points))
    Stau_St_lo <- rep(NA, length(time_points))
    Stau_St_up <- rep(NA, length(time_points))

    # we'll only populate results for valid indices (S(t) > 0 and tau <= max_time) - this shouldn't be needed but doesn't hurt to have it
    valid_idx <- which(S_t > 0 & tau <= max_time)

    # compute P(T > tau | T > t) only for valid indices
    Stau_St[valid_idx] <- ifelse(time_points[valid_idx] >= tau, 1, S_tau / S_t[valid_idx])
    Stau_St_lo[valid_idx] <- ifelse(time_points[valid_idx] >= tau, 1, S_tau_lo / S_t_up[valid_idx])
    Stau_St_up[valid_idx] <- ifelse(time_points[valid_idx] >= tau, 1, S_tau_up / S_t_up[valid_idx])

    # approximate aea under S(t) for all t and for t > tau
    sum_St <- sum(S_t * delta_t)
    sum_S_tau_end <- sum(S_t[time_points >= tau] * delta_t[time_points >= tau])

    # return a dataframe with results for this tau
    tibble(
      time = time_points,
      age_dist = A_t,
      tau = tau,
      Stau = S_tau, # this is useful for a little experiment
      sum_St = sum_St, # this is useful for a little experiment
      sum_S_tau_end = sum_S_tau_end, # this is useful for a little experiment
      Stau_St = Stau_St,
      Stau_St_lo = Stau_St_lo,
      Stau_St_up = Stau_St_up
    )
  })

  return(results)
}

################################################################################################
# stats
################################################################################################


# function to estimate Beta parameters from mean and CI limits
estimate_beta_params <- function(mean, lower, upper) {
  var <- upper - lower # approximate variance from CI width
  if (var < 1e-3 && abs(mean - 1) < 1e-8) { # handle the case of (P = 1, no uncertainty)
    alpha <- 1e6
    beta <- 1e-6
  } else if (var < 1e-3 && abs(mean) < 1e-8) { # handle the case of (P = 0, no uncertainty)
    alpha <- 1e-6
    beta <- 1e6
  } else {
    alpha <- ((1 - mean) / var^2 - 1 / mean) * mean^2
    beta <- alpha * (1 / mean - 1)
  }
  return(list(alpha = alpha, beta = beta))
}

# compare O_taus
pairwise_test <- function(cond_surv, # the output of cond_surv_mat
                          var_values = NULL, # groups to show (NULL = all)
                          tau_values = NULL, # if the above is specified, this must be a vector of the desired tau for each group, in the same order. If a number, applies to all
                          seed = 314,
                          B = 5000, # bootsrtap reps - should be high to fight noise
                          n_draw = 10 # draws per boostrap - should be low to fight KS overdoing its thing
) {
  set.seed(seed)

  # helper f: analytic KS distance (a deterministic metric of distance between two dsitrbutions)
  analytic_ks <- function(a1, b1, a2, b2) {
    # 'optimize()' finds the max or min value of f(x) for all x in interval
    optimize(function(x) abs(pbeta(x, a1, b1) - pbeta(x, a2, b2)), # f(x) is the difference between the two cdfs at x
      interval = c(0, 1), # the whole beta interval
      maximum = TRUE # look for the max value (this is the KS metric)
    )$objective
  }

  # helper f: compute p-value by bootstrap.
  bootstrap_p <- function(ks_deterministic, a1, b1, a2, b2) {
    # first assume a null beta from which both O(tau) values come from (avg of the two alphas and betas)
    a0 <- mean(c(a1, a2))
    b0 <- mean(c(b1, b2))

    # draw B pairs of samples from this null beta, and compute their KS statistic
    null_stats <- replicate(B, {
      x <- rbeta(n_draw, a0, b0)
      y <- rbeta(n_draw, a0, b0)
      suppressWarnings(ks.test(x, y)$statistic)
    })

    # how often is the KS statistic from the null pairs greater than or equal to the deterministic one?
    # or, how likely is it to see a deterministic distance this big if they actually come from the same generative process?
    # or, how likely is it that these two O(tau) values are fundamentally different?
    mean(null_stats >= ks_deterministic)
  }

  # fit beta distributions to each O(tau)
  if (is.null(var_values)) {
    if (length(tau_values) != 1) stop("if var_values is not specified, tau_values must be a single number.")

    filtered_cond_surv <- cond_surv$args$cond_surv %>%
      filter(tau == tau_values) %>%
      group_by(unique_label) %>%
      slice(1) %>% # one row per value (everything else is repeated or irrelevant)
      mutate( # fit beta distributions for each
        alpha = estimate_beta_params(O_tau, O_tau_lo, O_tau_up)$alpha,
        beta = estimate_beta_params(O_tau, O_tau_lo, O_tau_up)$beta
      )
  } else {
    if (length(tau_values) != length(var_values)) stop("tau_values and var_values must have the same length.")

    pairs <- tibble(tau = tau_values, unique_label = var_values) # match tau_values and var_values

    filtered_cond_surv <- cond_surv$cond_surv %>%
      semi_join(pairs, by = c("tau", "unique_label")) %>%
      group_by(unique_label) %>%
      slice(1) %>% # one row per value (everything else is repeated or irrelevant)
      mutate( # fit a beta dist to each
        alpha = estimate_beta_params(O_tau, O_tau_lo, O_tau_up)$alpha,
        beta = estimate_beta_params(O_tau, O_tau_lo, O_tau_up)$beta
      )
  }

  # compare the beta distributions

  # initialise results dataframe
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

  # pairwise comparisons
  for (i in 1:(nrow(filtered_cond_surv) - 1)) {
    for (j in (i + 1):nrow(filtered_cond_surv)) {
      # extract alphas and betas
      a1 <- filtered_cond_surv$alpha[i]
      b1 <- filtered_cond_surv$beta[i]
      a2 <- filtered_cond_surv$alpha[j]
      b2 <- filtered_cond_surv$beta[j]

      ks_stat <- analytic_ks(a1, b1, a2, b2) # deterministic KS distance between the two betas
      p_val <- bootstrap_p(ks_stat, a1, b1, a2, b2) # bootstrap-based p‑value - how likely are we to see this KS_stat if they're actually from the same beta?

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
  print(knitr::kable(results, digits = 3)) # print a pretty table
  invisible(results) # also allow assingment to a variable
}

compare_parametric_fits <- function(data # a lifetable
                                    , time_var # the variable with times
                                    , event_var # variable with status  (conventions from survival package)
                                    , dists = c("genf", "gengamma", "weibull", "gompertz", "gamma", "llogis", "lnorm", "exp") # distributions to fit, default is all of them
                                    , funcs = c("survival", "hazard", "cumhaz") # types of plots to produce, default is all of them
                                    , plot_title = "Parametric Survival Model Fits") {
  # useful for later
  pretty_dict <- c(
    genf     = "Generalised F",
    gengamma = "Generalised Gamma",
    weibull  = "Weibull",
    gompertz = "Gompertz",
    gamma    = "Gamma",
    llogis   = "Log-Logistic",
    lnorm    = "Log-Normal",
    exp      = "Exponential"
  )

  # fit all models
  fit_list <- list()
  for (dist in dists) {
    fit_list[[dist]] <- tryCatch(
      flexsurvreg(Surv(data[[time_var]], data[[event_var]]) ~ 1, data = data, dist = dist),
      error = function(e) {
        message(paste("Skipping", dist, "due to error:", e$message))
        NULL
      }
    )
  }

  # remove NULLS
  fit_list <- fit_list[!sapply(fit_list, is.null)]
  if (length(fit_list) == 0) { # dont continue if all NULL
    warning("No models were successfully fitted.")
    return(list(comparison = data.frame(), fits = NULL, plot = NULL))
  }

  # extract AIC, loglik, and pretty name for each fit
  aics <- sapply(fit_list, AIC)
  logliks <- sapply(fit_list, function(fit) fit$loglik)

  # generating a pvalue based on coxsnell residuals
  ks_pvals <- sapply(fit_list, function(fit) { # for each fitted model
    cs <- coxsnell_flexsurvreg(fit) # get the coxsnell residuals
    km_fit <- survfit(Surv(cs$est, cs$status) ~ 1) # fit a km curve to those residuals (*reason below)

    # data from km_fit
    obs_times <- km_fit$time
    emp_surv <- km_fit$surv

    # data from a theoretical population with exponential hazard
    theo_surv <- exp(-obs_times)

    # are they the same thing?
    ks_pval <- ks.test(emp_surv, theo_surv)$p.value
  })
  # *reason for that km fit
  # - CS residuals should resemble data from an Exp(1) without censoring
  # - fitting a KM to it effectively handles the censoring present
  # - check https://search.r-project.org/CRAN/refmans/flexsurv/html/coxsnell_flexsurvreg.html

  # use the dicitonary from above
  dist_names <- names(aics) # this is a named vector and we can use it
  dist_labels <- pretty_dict[dist_names]
  dist_labels[is.na(dist_labels)] <- dist_names[is.na(dist_labels)]

  # compile in a table, then sort it
  comparison <- data.frame(
    Distribution = dist_labels,
    AIC = aics,
    LogLik = logliks,
    KS_p = ks_pvals,
    row.names = NULL
  )
  comparison <- comparison[order(comparison$AIC), ]

  # reorder everything else (useful later)
  fit_list_ordered <- fit_list[comparison$Distribution %>% match(dist_labels)]
  # dist_ordered <- dist_names[order(aics)]
  # dist_labels_ordered <- dist_labels[order(aics)]


  # now let's produce plots in a table as well:
  n_dists <- length(fit_list_ordered)
  n_funcs <- length(funcs)
  op <- par(no.readonly = TRUE) # save old par settings
  par(mfrow = c(n_funcs, n_dists), mar = c(2, 2, 6, 1), oma = c(4, 4, 6, 2)) # margins and layout

  # # titles (top row)
  # for (label in dist_labels_ordered) {
  #   plot.new()
  #   title(main = label, cex.main = 1.2)
  # }

  # the actual plots
  for (func in funcs) {
    for (dist in dist_names) {
      fit <- fit_list[[dist]]
      if (func == funcs[1]) {
        plot(fit,
          type = func, xlab = "", ylab = func,
          main = paste0(pretty_dict[[dist]], "\nAIC: ", round(aics[dist], 2), "\nfit p-val: ", round(ks_pvals[dist], 3))
        )
      } else {
        plot(fit, type = func, xlab = "", ylab = func, main = NULL)
      }
    }
  }

  # extras
  mtext("Time", side = 1, outer = TRUE, line = 2)
  mtext("", side = 2, outer = TRUE, line = 2)
  mtext(plot_title, side = 3, outer = TRUE, line = 4, cex = 1.5)

  plot <- recordPlot()

  par(op) # reset plotting parameters

  invisible(list(fits = fit_list, comparison = comparison, plot = plot))
}

# relative entropy, information loss metric. Used in power analysis
kl_divergence <- function(tru, est, epsilon = 1e-10) {
  # Add a tiny epsilon to avoid log(0) and division by zero
  tru <- tru + epsilon
  est <- est + epsilon

  # Normalize again to ensure they sum to 1
  trh <- tru / sum(tru)
  est <- est / sum(est)

  sum(tru * log(tru / est))
}

################################################################################################
# plotting
################################################################################################

# km  curves - check survival_KM.r and functions.r for extended usage involving
# color control and annotations of median survival

create_survival_curve <- function(data,
                                  title,
                                  colors = NULL,
                                  linetypes = NULL,
                                  title_size = 20,
                                  conf_int = TRUE,
                                  group,
                                  time_var,
                                  event_var,
                                  xmax,
                                  legend = "none",
                                  median_surv = TRUE) {
  if (!(is.na(group))) {
    # convert group to a factor, define formula, fit KM model
    data[[group]] <- as.factor(data[[group]])
    surv_formula <- as.formula(paste0("Surv(", time_var, ", ", event_var, ") ~ ", group))
    fit <- survfit(surv_formula, data = data)

    # extract median survival for each level in data[[group]]
    surv_data <- summary(fit)$table
    annotations <- data.frame(
      group = rownames(surv_data),
      x = surv_data[, "median"]
    )
    labels <- levels(data[[group]]) # labels for the annotation
  } else {
    # convert group to a factor, define formula, fit KM model
    surv_formula <- as.formula("data$event ~ 1")
    fit <- survfit2(surv_formula)

    # extract median survival for each level in data[[group]]
    surv_data <- summary(fit)$table
    annotations <- data.frame(
      group = c("All"),
      x = surv_data["median"]
    )
    labels <- c("All")
  }

  use_linetype_aes <- !is.null(linetypes)

  # Plot the survival curves
  p <- survfit2(surv_formula, data = data) %>%
    ggsurvfit(type = "survival", linetype_aes = use_linetype_aes) +
    labs(
      x = "Time (days)",
      y = "Proportion alive"
    ) +
    xlim(0, xmax) + # Make sure all plots have the same axis limits
    ylim(0, 1) + # Make sure all plots have the same axis limits
    (if (!is.null(colors)) scale_fill_manual(values = colors) else NULL) + # Set fill colors based on the palette
    (if (!is.null(colors)) scale_color_manual(values = colors) else NULL) + # Set line colors based on the palette
    (if (!is.null(linetypes)) scale_linetype_manual(values = linetypes) else NULL) +
    theme(
      axis.text = element_text(size = 13), # Font size for axis ticks
      axis.title = element_text(size = 16), # Adjust font of labels
      legend.position = legend, # Try "bottom" or "none"
      legend.text = element_text(size = 13),
      legend.title = element_text(size = 16),
      plot.margin = margin(8, 8, 8, 8), # Plot margins (t, r, b, l)
      plot.title = element_text(size = title_size, hjust = 0.5) # Title settings
    ) +
    ggtitle(title)

  if (conf_int == TRUE) p <- p + add_confidence_interval(alpha = 0.1)


  if (median_surv == TRUE) {
    if (is.na(colors)) {
      stop("Median survival annotations need color specifications")
    }

    # Add a single "Median survival" annotation
    p <- p + annotate("text",
      x = 0, y = 0.25,
      label = "Median survival:",
      size = 8, color = "black", fontface = "bold", hjust = 0
    )

    # add the median survival for each level, and the lines
    for (i in 1:nrow(annotations)) {
      # add horizontal line
      p <- p + geom_segment(
        x = 0, xend = annotations$x[i], y = 0.5, yend = 0.5,
        color = colors[i], linetype = "solid"
      )

      # add vertical line
      p <- p + geom_segment(
        x = annotations$x[i], xend = annotations$x[i], y = 0.5, yend = 0,
        color = colors[i], linetype = "solid"
      )

      # Add annotation text
      p <- p + annotate("text",
        x = 0, y = 0.25 - i * 0.05, # Adjust y position for each label
        label = paste(labels[i], ": ", round(annotations$x[i], 1)),
        size = 8, color = colors[i], hjust = 0
      )
    }
  }
  print(p)
  invisible(p)
}


create_legend_plot <- function(
  data, colors, ltypes,
  group = "exposed"
) {
  # convert the group to a factor and define formula
  data[[group]] <- as.factor(data[[group]])
  surv_formula <- as.formula(paste("data$event ~ data$", group, sep = ""))

  # define labels
  if (group == "exposed") {
    labels <- levels(data[[group]]) # to check the order of the labels
    l_title <- "Exposure:"
  } else if (group == "species") {
    labels <- levels(data[[group]]) # to check the order of the labels
    l_title <- "Species:"
  } else if (group == "temp_mean") {
    labels <- levels(data[[group]]) # to check the order of the labels
    l_title <- "Mean temperature:"
  } else if (group == "temp_range") {
    labels <- levels(data[[group]]) # to check the order of the labels
    l_title <- "Temperature oscillations:"
  } else {
    labels <- levels(data[[group]]) # to check the order of the labels
    l_title <- "Group:"
  }
  # maybe double check in which order each should be... but it seems fine to me like this

  # first make a plot
  p <- survfit2(surv_formula) %>%
    ggsurvfit() +
    scale_fill_manual(values = colors, labels = labels) + # Set fill colors and labels
    scale_color_manual(values = colors, labels = labels) + # Set line colors and labels
    add_confidence_interval() +
    theme(
      legend.text = element_text(size = 10), # Adjust legend text size
      legend.position = "bottom", # Position the legend at the bottom
      legend.title = element_text(size = 11), # Adjust legend title size
    ) +
    guides(
      fill = guide_legend(title = l_title),
      color = guide_legend(title = l_title)
    )

  # then xtract the legend
  grob <- ggplotGrob(p) # convert the plot to a grob
  # print(grob$layout$name) # find where the legend is (guide-box-bottom)
  positions <- which(grob$layout$name == "guide-box-bottom")
  legend <- grob$grobs[[positions[1]]] # extract the  legend
  # as_ggplot(legend) # plot the legend
  return(legend)
}


# heatmap for many pairwise p-value comparisons
ps_heatmap <- function(results, # heatmap values
                       surv, # for margin annotations
                       title, # main title of the whole plot
                       sig_threshold = 0.05, # significance threshold
                       ex_colors = c("Exposed" = "skyblue", "Control" = "tomato"),
                       tm_colors = c("21" = "darkgreen", "27" = "gold"),
                       tr_colors = c("0" = "purple", "6" = "orange")) {
  # dictionary for margin annotations
  ann <- surv[, (names(surv) %in% c("pot", "temp_mean", "temp_range", "exposed"))] %>%
    unique()
  rownames(ann) <- ann$pot
  ann$pot <- NULL # we don't ned this nomore

  # prepare gambiae and coluzzii matrices - two heatmaps, two matrices
  gam <- filter(results, species == "An. gambiae")

  mat_gam <- expand_grid(Group1 = unique(gam$Group1), Group2 = unique(gam$Group2)) %>% # create grid of combinations
    left_join(gam, by = c("Group1", "Group2")) %>% # expand original with NAs
    select(c(Group1, Group2, KS_p_log10)) %>% # get rid of useless columns
    pivot_wider(names_from = Group2, values_from = KS_p_log10) %>% # convert to square matrix
    column_to_rownames("Group1") %>%
    as.matrix()


  col <- filter(results, species == "An. coluzzii")

  mat_col <- expand_grid(Group1 = unique(col$Group1), Group2 = unique(col$Group2)) %>% # create grid of combinations
    left_join(col, by = c("Group1", "Group2")) %>% # expand original with NAs
    select(c(Group1, Group2, KS_p_log10)) %>% # get rid of useless columns
    pivot_wider(names_from = Group2, values_from = KS_p_log10) %>% # convert to square matrix
    column_to_rownames("Group1") %>%
    as.matrix()

  # colormap for annotations (change to match those in KMs... or the other way around?)
  colormap <- list(
    exposed = ex_colors,
    temp_mean = tm_colors,
    temp_range = tr_colors
  )

  # subset and reorder annotation dictionary for each matrix
  ann_for_gam <- ann[rownames(mat_gam), , drop = FALSE]
  ann_for_col <- ann[rownames(mat_col), , drop = FALSE]
  ann_for_col_reordered <- ann[rownames(mat_col), c("temp_range", "temp_mean", "exposed"), drop = FALSE] # reorder columns for right annotations

  # annotation objects
  row_annot_gam <- rowAnnotation(
    df = ann_for_gam,
    col = colormap,
    show_annotation_name = FALSE
  )

  col_annot_gam <- HeatmapAnnotation(
    df = ann_for_gam,
    col = colormap,
    show_annotation_name = FALSE
  )

  row_annot_col <- rowAnnotation(
    df = ann_for_col_reordered,
    col = colormap,
    show_annotation_name = FALSE
  )

  col_annot_col <- HeatmapAnnotation(
    df = ann_for_col,
    col = colormap,
    show_annotation_name = FALSE
  )

  # colorbar
  max <- max(c(mat_gam, mat_col), na.rm = TRUE) # maximum value in either dataset (shared colorbar)
  mid <- -log10(sig_threshold) # change of color at the significance threshold
  min <- -max(c(mat_gam, mat_col), na.rm = TRUE) # symmetry around midpoint keeps color intensity informative of magnitude

  color_bar <- colorRamp2(c(min, mid, max), c("blue", "white", "red"))

  # the heatmaps
  ht1 <- Heatmap(mat_gam,
    name = "Significance (-log(p))",
    col = color_bar,
    top_annotation = col_annot_gam,
    left_annotation = row_annot_gam,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    column_title = "An. gambiae",
    row_names_side = "left",
    column_title_gp = gpar(fontsize = 18),
    row_title_gp = gpar(fontsize = 18),
    row_names_gp = gpar(fontsize = 18),
    column_names_gp = gpar(fontsize = 18)
  )

  ht2 <- Heatmap(mat_col,
    name = "Significance (-log(p))",
    col = color_bar,
    top_annotation = col_annot_col,
    right_annotation = row_annot_col,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    column_title = "An. coluzzii",
    column_title_gp = gpar(fontsize = 18),
    row_title_gp = gpar(fontsize = 18),
    row_names_gp = gpar(fontsize = 18),
    column_names_gp = gpar(fontsize = 18)
  )

  # combine
  draw(ht1 + ht2,
    heatmap_legend_side = "bottom",
    merge_legend = TRUE,
    padding = unit(c(2, 2, 20, 2), "mm")
  )

  # title
  grid.text(title,
    y = unit(1, "npc") - unit(6, "mm"),
    gp = gpar(fontsize = 22, fontface = "bold")
  )
}

# f(t) with color for censoring separated by some variable of interest
# to do: allow a third colour to show left censoring

f_t_census <- function(data, # the dataset
                       time_var, # name of the time variable (string)
                       event_var, # name of the event variable (string) - values will show in the legend
                       var, # variable to facet by (string)
                       nbins = max(data[[time_var]]), # max(t) by default
                       ncol = 1) {
  ft <- ggplot(data, aes(x = !!sym(time_var), fill = !!sym(event_var))) +
    geom_histogram(position = "dodge", bins = nbins, alpha = 1) +
    scale_fill_manual(values = c("Dead" = "red", "Censored" = "green")) +
    labs(
      title = "Histogram of survival times",
      x = "Time (days)",
      y = "Frequency",
      fill = "Exit reason"
    ) +
    facet_wrap(as.formula(paste("~", var)), ncol = ncol) + # parse var
    # theme_minimal() +
    theme(
      plot.title = element_text(size = 20, hjust = 0.5),
      axis.title = element_text(size = 16),
      axis.text = element_text(size = 13),
      strip.text = element_text(size = 16),
      legend.title = element_text(size = 16),
      legend.text = element_text(size = 13)
    )
  print(ft) # display the plot
  invisible(ft) # return without printing - so it can be saved in a variable
}

# visualise conditional survival (heatmaps and line plots)

cond_surv_plot <- function(cond_surv, # the output of cond_surv_mat
                           var_value = NULL, # optional filter - must be one of the values of var
                           tau_value = NULL, # optional filter - must be a valid value within the range taus
                           color_var = NULL, # variable to color by (if NULL, use var)
                           linetype_var = NULL, # variable to linetype by (if NULL, no linetype)
                           markershape_var = NULL, # variable to use for marker shape (if NULL, no shape)
                           ncol_hm = 1, # number of columns if doing a heatmap (no filters)
                           collapse_time = FALSE # if TRUE, plot P(T>tau) instead of P(T>tau|T>t) for all t - collapsing t using A(t) as weights
) {
  # aggregate if appropriate
  if (cond_surv$args$aggregate_by_cats) {
    # compute maximum death observed in each group (levels in var)
    max_values <- cond_surv$args$data %>%
      filter(!!sym(cond_surv$args$event_var) == 1) %>% # only deaths
      group_by(!!sym(cond_surv$args$var), !!!syms(cond_surv$args$cats)) %>% # for each group of interest
      summarise(max_t = max(!!sym(cond_surv$args$time_var), na.rm = TRUE)) %>% # max_t is the latest death observed
      tidyr::unite("unique_label", all_of(c(cond_surv$args$var, cond_surv$args$cats)),
        sep = "_", remove = FALSE # in the future this could be made customisable for better legends
      )
  } else {
    # compute maximum death observed in each group (levels combinations in var and cats)
    max_values <- cond_surv$args$data %>%
      filter(!!(sym(cond_surv$args$event_var)) == 1) %>% # only deaths
      group_by(!!sym(cond_surv$args$var)) %>% # for each group of interest
      summarise(max_t = max(!!sym(cond_surv$args$time_var), na.rm = TRUE)) %>% # max_t is the latest death observed
      mutate(
        unique_label = paste0(!!sym(cond_surv$args$var)) # in the future this could be made customisable for better legends
      )
  }
  # useful for later
  max_time <- cond_surv$args$data %>%
    summarise(max_value = max(!!sym(cond_surv$args$time_var[[1]]), na.rm = TRUE))

  # if either filter is specified
  if (!is.null(tau_value) && collapse_time == FALSE) { # if tau_value is specified and collapse_time is FALSE the plot is P(T>tau|T>t) against t

    filtered_cond_surv <- cond_surv$cond_surv %>%
      filter(tau == tau_value)

    # for later
    title <- bquote("Time-dependant probability of outliving " ~ .(tau_value))

    if (!is.null(var_value)) { # if var_value is specified, filter by it too
      filtered_cond_surv <- filtered_cond_surv %>%
        filter(unique_label == var_value)

      # update
      title <- bquote("Time-dependant probability of outliving " ~ .(tau_value))
    }

    if (is.null(color_var)) { # if color_var isn't specified, use unique_var
      color_var <- "unique_label"
    }

    # plot the curve. If linetype_var is NULL:
    if (is.null(linetype_var)) { # if linetype_var isn't specified, all are solid
      plot <- ggplot(filtered_cond_surv, aes(
        x = time,
        y = Stau_St,
        group = unique_label,
        color = !!sym(color_var)
      ))
    } else { # but if it is, then use it
      plot <- ggplot(filtered_cond_surv, aes(
        x = time,
        y = Stau_St,
        group = unique_label,
        color = !!sym(color_var),
        linetype = !!sym(linetype_var)
      ))
    }

    # plot the curve
    plot <- plot +
      geom_line() +
      geom_ribbon(aes(ymin = Stau_St_lo, ymax = Stau_St_up, group = unique_label, fill = !!sym(color_var)), alpha = 0.2) +
      labs(x = "Time", y = bquote(O[tau == .(tau_value)](t)), title = title) +
      # scale_x_continuous(breaks = 1:30) +
      ylim(0, 1) +
      theme_minimal() +
      theme(
        plot.title = element_text(size = 18),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 13),
        strip.text = element_text(size = 16),
        legend.text = element_text(size = 13),
        legend.title = element_blank(),
        # legend.position = "bottom"
      )

    print(plot) # display the plot
    invisible(plot) # return without printing - so it can be saved in a variable
  } else if (!is.null(tau_value) && collapse_time == TRUE) { # if tau_value is specified and collapse_time is TRUE the plot is P(T>tau) for each var level

    # if var_values is NULL, tau_value is a number and is the same for all values of the variable of interest
    if (var_value %>% is.null()) {
      # tau_value must be a number
      if (!is.numeric(tau_value) || length(tau_value) != 1) {
        stop("tau_value must be a single numeric value or a vector of the same length as var_value")
      }
      # and this number can be in the plot title
      title <- bquote("Probability of outliving " ~ .(tau_value) ~ " when age is unknown")

      # filter cond_surv
      filtered_cond_surv <- cond_surv$cond_surv %>%
        filter(tau == tau_value)

      if (is.null(color_var)) { # if color_var isn't specified, use unique_var
        color_var <- "unique_label"
      }

      if (is.null(markershape_var)) { # if markershape_var isn't specified, no shape
        shape_mapping <- NULL
      } else {
        shape_mapping <- aes(shape = !!(sym(markershape_var)))
      }

      # plot it
      plot <- ggplot(filtered_cond_surv, aes(x = unique_label, y = O_tau, color = !!(sym(color_var)))) +
        geom_point(shape_mapping, position = position_dodge(width = 0.5), size = 3) +
        geom_errorbar(aes(ymin = O_tau_lo, ymax = O_tau_up), width = 0.2, position = position_dodge(width = 0.5)) +
        labs(x = "Group", y = bquote(O[tau == .(tau_value)]), title = title) +
        ylim(0, 1) +
        theme_minimal() +
        theme(
          plot.title = element_text(size = 18),
          axis.title = element_text(size = 16),
          axis.text.x = element_text(angle = 13, hjust = 1),
          axis.text = element_text(size = 13),
          strip.text = element_text(size = 16),
          strip.background = element_rect(fill = "#bdbdef", color = "white", size = 2),
          legend.text = element_text(size = 13),
          legend.title = element_blank(),
          # legend.position = "bottom"
        )

      print(plot) # display the plot
      invisible(plot) # return without printing - so it can be saved in a variable

      # if var_value has been specified, tau_values should be a vector of taus for each level, in the same order
    } else {
      # tau_value must be the same length
      if (length(tau_value) != length(var_value)) {
        stop("tau_value must be a single numeric value or a vector of the same length as var_value")
      }
      # and the title is a bit more generic
      title <- bquote("Probability of outliving " * tau * " when age is unknown")

      # first make a tibble with the pairs of tau and var
      pairs <- tibble(
        tau = tau_value,
        unique_label = var_value # rlang's !!sym()
      )

      # now filter cond_surv
      filtered_cond_surv <- cond_surv$cond_surv %>%
        semi_join(pairs, by = c("tau", "unique_label")) # nice and neat

      # create a variable that joiins the group and the tau choice for labelling, with proper handling of characters for downstream printing of the greek symbol tau
      filtered_cond_surv <- filtered_cond_surv %>%
        mutate(
          label = paste0("plain('", unique_label, "') ~ (tau==", tau, ")")
        )

      if (is.null(color_var)) { # if color_var isn't specified, use unique_var
        color_var <- "unique_label"
      }

      if (is.null(markershape_var)) { # if markershape_var isn't specified, no shape
        shape_mapping <- NULL
      } else {
        shape_mapping <- aes(shape = !!(sym(markershape_var)))
      }

      # plot it
      plot <- ggplot(filtered_cond_surv, aes(x = label, y = O_tau, color = !!(sym(color_var)))) +
        geom_point(shape_mapping, position = position_dodge(width = 0.5), size = 3) +
        geom_errorbar(aes(ymin = O_tau_lo, ymax = O_tau_up), width = 0.2, position = position_dodge(width = 0.5)) +
        labs(x = "Group", y = bquote(O[tau]), title = title) +
        scale_x_discrete(labels = function(x) parse(text = x)) + # to print the tau symbol, not the "tau" string
        ylim(0, 1) +
        theme_minimal() +
        theme(
          plot.title = element_text(size = 18),
          axis.title = element_text(size = 16),
          axis.text.x = element_text(angle = 13, hjust = 1),
          axis.text = element_text(size = 13),
          strip.text = element_text(size = 16),
          strip.background = element_rect(fill = "#bdbdef", color = "white", size = 2),
          legend.text = element_text(size = 13),
          legend.title = element_blank(),
          # legend.position = "bottom"
        )

      print(plot) # display the plot
      invisible(plot) # return without printing - so it can be saved in a variable
    }
  } else if (is.null(tau_value) && collapse_time == FALSE) { # if tau_value is not specified and collapse_time is FALSE, the plot is a heatmap.

    filtered_cond_surv <- cond_surv$cond_surv
    title <- bquote("Conditional survival per group")

    if (!is.null(var_value)) { # if var_value is specified, filter by it too
      filtered_cond_surv <- filtered_cond_surv %>%
        filter(unique_label == var_value)

      # update title
      title <- bquote("Conditional survival. Group: " ~ filtered_cond_surv$unique_label)
    }
    n_groups <- length(unique(filtered_cond_surv$unique_label)) # number of groups to plot (useful later)
    n_rows <- ceiling(n_groups / ncol_hm) # number of rows in the grid, useful later

    # heatmap
    plot <- ggplot(filtered_cond_surv, aes(x = time, y = tau, fill = Stau_St)) +
      geom_tile(color = NA) +
      geom_vline(data = max_values, aes(xintercept = max_t), color = "#00f891", linetype = "solid", linewidth = 1) + # adds a vertical line at the time of the last death observed
      scale_fill_gradientn(colors = brewer.pal(9, "YlOrRd"), na.value = "grey") +
      labs(
        title = title,
        x = "Age (t)",
        y = bquote("Age to outlive (" * tau * ")"),
        fill = bquote("P(T > " * tau * " | T > t)")
      ) +
      facet_wrap(~unique_label, scales = "free", ncol = ncol_hm) +
      theme_minimal() +
      theme(
        aspect.ratio = 1,
        plot.title = element_blank(),
        axis.title = element_text(size = max(6, 11 - 1.5 * n_rows) + 1), # depends on number of groups
        axis.text = element_text(size = max(6, 11 - 1.5 * n_rows)), # depends on number of groups
        strip.text = element_text(size = max(6, 11 - 1.5 * n_rows)), # depends on number of groups
        legend.text = element_text(size = 10),
        legend.title = element_blank(),
      ) +
      guides(
        fill = guide_colorbar(
          barwidth = 0.7,
          barheight = n_rows * 4, # depends on number of rows
          title.hjust = 0.5
        )
      )

    print(plot) # display the plot
    invisible(plot) # return without printing - so it can be saved in a variable
  } else if (is.null(tau_value) && collapse_time == TRUE) { # if tau_value is not specified and collapse_time is TRUE, the plot is P(T>tau) against tau

    filtered_cond_surv <- cond_surv$cond_surv
    title <- bquote("Probability of outliving " * tau * " when age is unknown")

    if (!is.null(var_value)) { # if var_value is specified, filter by it too
      filtered_cond_surv <- filtered_cond_surv %>%
        filter(unique_label == var_value)

      # update title
      title <- bquote("P(T>" ~ tau ~ ") when age is unknown. Group: " ~ filtered_cond_surv$unique_label)
    }

    if (is.null(color_var)) { # if color_var isn't specified, use unique_var
      color_var <- "unique_label"
    }

    # plot the curve. If linetype_var is NULL:
    if (is.null(linetype_var)) { # if linetype_var isn't specified, all are solid
      plot <- ggplot(filtered_cond_surv, aes(
        x = tau,
        y = O_tau,
        group = unique_label,
        color = !!sym(color_var)
      ))
    } else { # but if it is, then use it
      plot <- ggplot(filtered_cond_surv, aes(
        x = tau,
        y = O_tau,
        group = unique_label,
        color = !!sym(color_var),
        linetype = !!sym(linetype_var)
      ))
    }

    # plot the curve
    plot <- plot +
      geom_line() +
      geom_ribbon(aes(ymin = O_tau_lo, ymax = O_tau_up, group = unique_label, fill = !!sym(color_var)), alpha = 0.2) +
      labs(x = expression(tau), y = expression(O(tau)), title = title) +
      xlim(0, max_time[[1]]) +
      ylim(0, 1) +
      theme_minimal() +
      theme(
        plot.title = element_text(size = 18),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 13),
        strip.text = element_text(size = 16),
        strip.background = element_rect(fill = "#bdbdef", color = "white", size = 2),
        legend.text = element_text(size = 13),
        legend.title = element_blank(),
        # legend.position = "bottom"
      )

    print(plot) # display the plot
    invisible(plot) # return without printing - so it can be saved in a variable
  }
}

# visualise O_tau computation
o_tau_computation <- function(
  cond_surv # the output of cond_surv_mat
  , var_value,
  tau_value
) {
  # filter cond_surv
  filtered_cond_surv <- cond_surv$cond_surv %>%
    filter(cond_surv$cond_surv$unique_label == var_value & tau == tau_value)

  # chose the scale for the age_dist (otherwise it will be too small)
  scale <- 1 / filtered_cond_surv$age_dist[1]

  # draw the O_tau(t) curve alone first
  otau <- ggplot(filtered_cond_surv, aes(x = time)) +
    geom_line(aes(x = time, y = Stau_St, color = "P(T>τ|T>t)"), size = 2) + # conditional survival curve
    geom_ribbon(aes(ymin = Stau_St_lo, ymax = Stau_St_up, fill = "95%CI"), alpha = 0.5) + # 95% CI of conditional survival curve
    labs(x = "t", y = paste0("P(T>", tau_value, "|T>t)"), title = paste0(filtered_cond_surv$unique_label, ", τ = ", tau_value)) +
    scale_color_manual(
      name = "",
      values = c("P(T>τ|T>t)" = "black")
    ) +
    scale_fill_manual(
      name = "",
      values = c("95%CI" = "gray")
    ) +
    ylim(c(0, 1)) +
    theme_minimal(base_family = "serif") +
    theme(
      plot.title = element_text(size = 28, face = "bold"),
      axis.title = element_text(size = 24),
      axis.text = element_text(size = 22),
      legend.text = element_text(size = 20),
      axis.line = element_line(size = 1, color = "black"),
      axis.ticks = element_line(size = 1, color = "black")
    )


  # add the age distribution
  true_dist <- otau +
    stat_summary(fun = mean, geom = "bar", aes(y = age_dist * scale, fill = "A(t)"), alpha = 0.2) + # average over each time value, in case our group var is not exhaustive of all factors
    scale_fill_manual(
      name = "Legend",
      values = c("A(t)" = "blue")
    ) +
    guides(
      fill = guide_legend(order = 1),
      color = guide_legend(order = 2)
    ) +
    scale_y_continuous(sec.axis = sec_axis(~ . / scale, name = "A(t)")) # add secondary y-axis for O_tau(t)

  # add the product of the two
  true_dist_otau_t <- true_dist +
    geom_line(aes(x = time, y = Stau_St * age_dist * scale, color = "P(T>τ|T>t) * A(t)"), size = 2) +
    scale_color_manual(
      name = "",
      values = c("P(T>τ|T>t)" = "black", "P(T>τ|T>t) * A(t)" = "red")
    ) +
    scale_fill_manual(
      name = "Legend",
      values = c("A(t)" = "blue")
    ) +
    guides(
      fill = guide_legend(order = 1),
      color = guide_legend(order = 2)
    )

  # add the shaded area
  O_tau_plot <- true_dist_otau_t +
    geom_area(aes(x = time, y = Stau_St * age_dist * scale, fill = "O_tau"), alpha = 0.2) +
    scale_fill_manual(
      name = "Legend",
      values = c("A(t)" = "blue", "O_tau" = "red"),
      labels = c("A(t)", expression(O[tau]))
    ) +
    scale_color_manual(
      name = "",
      values = c("P(T>τ|T>t)" = "black", "P(T>τ|T>t) * A(t)" = "red")
    ) +
    guides(
      fill = guide_legend(order = 1),
      color = guide_legend(order = 2)
    ) +
    annotate(
      "text",
      x = Inf, y = 0.6,
      label = expression(O[tau] == " red area"),
      hjust = 1.1, vjust = 2,
      size = 12,
    ) +
    theme_minimal(base_family = "serif") +
    theme(
      plot.title = element_text(size = 28, face = "bold"),
      axis.title = element_text(size = 24),
      axis.text = element_text(size = 22),
      legend.text = element_text(size = 20),
      legend.title = element_text(size = 22),
      axis.line = element_line(size = 1, color = "black"),
      axis.ticks = element_line(size = 1, color = "black")
    )


  print(O_tau_plot) # display the last plot
  invisible(list(otau, true_dist, true_dist_otau_t, O_tau_plot)) # return the three plots
}

# visualise O_tau for all taus
# Deprecated, this is now a case of cond_surv_plot (tau_values = NULL, collapse_time = TRUE)
o_tau_curve <- function(cond_surv # the output of cond_surv_mat
                        , facet = FALSE # switch for facet_wrap
                        , facet_var = "unique_label" # default makes sense but can be changed (string)
                        , ncol = 1 # for facet_wrap
) {
  max_time <- cond_surv$args$data %>%
    summarise(max_value = max(!!sym(cond_surv$args$time_var[[1]]), na.rm = TRUE))


  otau_curve <- ggplot(cond_surv$cond_surv, aes(x = tau, y = O_tau, color = unique_label)) +
    geom_line(size = 1.5) +
    geom_ribbon(aes(ymin = O_tau_lo, ymax = O_tau_up, fill = !!(sym(facet_var))), alpha = 0.2) +
    labs(x = "tau", y = paste("P(T>tau) when age is unknown"), title = paste("Prob. random individual will outlive tau")) +
    xlim(0, max_time[[1]]) +
    ylim(0, 1) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 16, face = "bold"),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 11),
      strip.text = element_text(size = 11),
      strip.background = element_rect(fill = "#bdbdef", color = "white", size = 2),
      legend.text = element_text(size = 10),
      legend.title = element_text(size = 11),
      # legend.position = "bottom"
    )

  if (facet) {
    p <- p + facet_wrap(as.formula(paste("~", facet_var)), scales = "free", ncol = ncol)
  }

  print(otau_curve) # display the plot
  invisible(otau_curve) # return without printing - so it can be saved in a variable
}

# visualise O_tau for a specific tau
o_tau_scatter <- function(cond_surv # the output of cond_surv_mat
                          , tau_value # the tau value we want
                          , var_values = NULL # if !NULL, must be a vector of all possible values in the variable of interest
                          , tau_values = NULL # if the above is specified, this must be a vector of the desired tau for each value, in the same order
                          , facet = FALSE # switch for facet_wrap
                          , facet_var = "unique_label" # default makes sense but can be changed (string)
                          , ncol = 1 # for facet_wrap
) {
  # if var_values is NULL, tau_value is a number and is the same for all values of the variable of interest
  if (var_values %>% is.null()) {
    # filter cond_surv
    filtered_cond_surv <- cond_surv$cond_surv %>%
      filter(tau == tau_value)

    # plot it
    otau_scatter <- ggplot(filtered_cond_surv, aes(x = !!sym(facet_var), y = O_tau, color = !!(sym(cond_surv$args$var)))) +
      geom_point(position = position_dodge(width = 0.5), size = 4) +
      geom_errorbar(aes(ymin = O_tau_lo, ymax = O_tau_up), width = 0.2, position = position_dodge(width = 0.5)) +
      labs(x = "Group", y = paste0("P(T>", tau_value, ") when age is unknown"), title = paste0("Prob. random individual will outlive ", tau_value, "d")) +
      ylim(0, 1) +
      theme_minimal() +
      theme(
        plot.title = element_text(size = 16, face = "bold"),
        axis.title = element_text(size = 12),
        axis.text.x = element_text(angle = 11, hjust = 1),
        axis.text = element_text(size = 11),
        strip.text = element_text(size = 11),
        strip.background = element_rect(fill = "#bdbdef", color = "white", size = 2),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 11),
        # legend.position = "bottom"
      )

    if (facet) {
      p <- p + facet_wrap(as.formula(paste("~", facet_var)), scales = "free", ncol = ncol)
    }


    print(otau_scatter) # display the plot
    invisible(otau_scatter) # return without printing - so it can be saved in a variable


    # if var_values has been specified, tau_values should be a vector of taus for each level, in the same order
  } else {
    # first make a tibble with the pairs of tau and var
    pairs <- tibble(
      tau = tau_values,
      unique_label = var_values # rlang's !!sym()
    )

    # now filter cond_surv
    filtered_cond_surv <- cond_surv$cond_surv %>%
      semi_join(pairs, by = c("tau", "unique_label")) # nice and neat

    # create a variable that joiins the group and the tau choice for labelling
    filtered_cond_surv <- filtered_cond_surv %>%
      mutate(label = paste0(!!sym(facet_var), " (tau = ", tau, ")"))

    facet_var <- "label" # update facet_var

    # plot it
    otau_scatter <- ggplot(filtered_cond_surv, aes(x = !!sym(facet_var), y = O_tau, color = !!(sym(cond_surv$args$var)))) +
      geom_point(position = position_dodge(width = 0.5), size = 4) +
      geom_errorbar(aes(ymin = O_tau_lo, ymax = O_tau_up), width = 0.2, position = position_dodge(width = 0.5)) +
      labs(x = "Group", y = paste0("P(T>tau) when age is unknown"), title = paste0("Prob. random individual will outlive age tau")) +
      ylim(0, 1) +
      theme_minimal() +
      theme(
        plot.title = element_text(size = 16, face = "bold"),
        axis.title = element_text(size = 12),
        axis.text.x = element_text(angle = 11, hjust = 1),
        axis.text = element_text(size = 11),
        strip.text = element_text(size = 11),
        strip.background = element_rect(fill = "#bdbdef", color = "white", size = 2),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 11),
        # legend.position = "bottom"
      )

    if (facet) {
      p <- p + facet_wrap(as.formula(paste("~", facet_var)), scales = "free", ncol = ncol)
    }


    print(otau_scatter) # display the plot
    invisible(otau_scatter) # return without printing - so it can be saved in a variable
  }
}
