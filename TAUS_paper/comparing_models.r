# this script compares taus cox and parametric survival models

##################################################################################################################################
# Environment
##################################################################################################################################
# a few packages
packages <- c(
  "tidyverse",
  "tibble",
  "dplyr",
  "rstudioapi",
  "remotes" # allows us to install TAUS from gh
)

for (i in packages) {
  if (!require(i, character.only = TRUE)) install.packages(i)
  library(i, character.only = TRUE)
}

# a couple useful fucntions for setting up the right working directory:

# returns full path to parent folder of current script
whereami <- function() {
  this_file <- commandArgs() %>%
    tibble::enframe(name = NULL) %>%
    tidyr::separate(col = value, into = c("key", "value"), sep = "=", fill = "right") %>%
    dplyr::filter(key == "--file") %>%
    dplyr::pull(value)
  if (length(this_file) == 0) {
    this_file <- rstudioapi::getSourceEditorContext()$path
  }
  return(dirname(this_file))
}

# set working directory to be the parent of the above
setwd_grandparent <- function() {
  this_script <- whereami()
  setwd(dirname(this_script))
  # print working directory
  print("Current working directory:")
  print(getwd())
}

# apply it
setwd_grandparent() # working directory should be the folder of where TAUS is

# install TAUS fro GH
library(TAUS) # doesnt work at first (ensure we don't have it at this point)
# remove.packages("TAUS") # if the above works, remove and reinstall
remotes::install_github("casasgomezuribarri/TAUS") # then we install from GH
library(TAUS)

# some more packages:
packages <- c(
  "RColorBrewer",
  "betareg",
  "effectsize",
  "survminer",
  "nortest",
  "scales",
  "ggsurvfit",
  "cramer",
  "ggtext",
  "survival",
  "flexsurv",
  "dplyr",
  "plotly"
)

for (i in packages) {
  if (!require(i, character.only = TRUE)) install.packages(i)
  library(i, character.only = TRUE)
}

##################################################################################################################################
# load the data
##################################################################################################################################
set.seed(1984) # reproducibility of random processes

# choose the dataset (simulation or realworld excmples - comment/uncomment as appropriate)

dataset <- "simu" # simulated dataset
# dataset <- "dros1" # example 2: survival by temperature
# dataset <- "dros2" # example 1: survival by country


# datasets were named in chronological order of analysis
# examples were named in order of appearance in the manuscript (simpler to more complex)

# Yes, I am aware that the names are confusing.
# No, I'm not going to rename everything (bound to make mistakes...)

# from now on everything should run without manual input:

if (dataset == "simu") {
  surv <- read.csv("TAUS/Data/synth/simu.csv") # dataset

  # realistic conditions
  p <- 0.4 # proporiton censored
  n <- 1000 # sample size
  f <- 200 # follow-up period

  # format dataset into appropriate lifetable
  surv <- surv %>%
    mutate(
      span = round(death_time) - birth_time, # calculate lifespan
      habitat = as.factor(surv$habitat),
      mutation = as.factor(surv$mutation),
      event = sample(c(1), size = n(), replace = TRUE) # all deaths observed
    )

  # censor the data
  surv_tmp <- surv %>%
    mutate(
      event = sample(c(0, 1), size = n(), replace = TRUE, prob = c(p, 1 - p)) # 0 - censored, 1 - death observed. p% censored
    ) %>%
    mutate(
      time = if_else(event == 0, round(span * runif(n())), span), # if censored, span observed is a random fraction of the total lifespan
      unique_var = paste0(habitat, "_", mutation) # create a unique variable for each combination of habitat and mutation, might be useful
    )

  # subsample dataset
  surv_tmp <- surv_tmp[sample(nrow(surv), n), ]

  # limit follow-up period
  surv_tmp <- surv_tmp %>%
    mutate(
      event = if_else(time > f, 0L, event), # change event to 0 if time > f
      time = if_else(time > f, f, time) # change time to f it time > f
    )

  # define a few useful things
  var <- "unique_var" # variable of interest
  cats <- c("habitat") # other categorical variables in the dataset
  time_var <- "time"
  event_var <- "event"
  surv_tmp$event_factor <- factor(surv_tmp[[event_var]], levels = c(0, 1), labels = c("Censored", "Dead"))
  event_factor <- "event_factor" # event variable as a factor
  max_time <- max(surv_tmp[[time_var]])
  plot_names <- "simu"
} else if (dataset == "dros1") {
  surv_tmp <- read.csv("TAUS/Data/online/dros1.csv") # dataset

  surv_tmp <- surv_tmp %>%
    mutate( # extract replicate and phenotype from selection line
      rep = str_sub(SelectionLine, -1, -1),
      phen = str_sub(SelectionLine, 1, -2),
      temp = paste0(Treatment, "°C"),
      group = paste0(temp, "_", phen)
    ) %>%
    group_by(SelectionLine) %>% # for each specific selectionline group...
    mutate(
      event = if_else(Longevity_Days == -9, 0, 1), # if age is -9, the fly outlived the experiment (censored, event = 0). Otherwise death was observed (event = 1)
      time = if_else(Longevity_Days == -9, max(Longevity_Days, na.rm = TRUE), Longevity_Days) # time variable should be like this, not like what they had
    ) %>%
    ungroup() %>% # go back to original format
    subset(phen == "Control") %>% # remove genetic variation
    as.data.frame() # group_by() converted it to a tibble

  # define a few useful things

  var <- "temp" # variable of interest
  cats <- c("rep", "phen") # other categorical variables in the dataset
  time_var <- "time"
  event_var <- "event"
  surv_tmp$event_factor <- factor(surv_tmp[[event_var]], levels = c(0, 1), labels = c("Censored", "Dead"))
  event_factor <- "event_factor" # event variable as a factor
  max_time <- max(surv_tmp[[time_var]])
  plot_names <- "dros1"
} else if (dataset == "dros2") {
  surv_tmp <- read.csv("TAUS/Data/online/dros2.csv") # dataset

  surv_tmp <- surv_tmp %>% # remove unnecessary columns
    select(
      -starts_with("Fec"),
      -starts_with("Hatch"),
      -starts_with("Ovariole"),
      -starts_with("X")
    ) %>%
    mutate(
      event = 1, # all files were followed until death
      time = Female.Age.at.Death,
      event_factor = factor(event, levels = c(0, 1), labels = c("Censored", "Dead")) # in factor format mkaes plotting easier
    )

  # define a few useful things

  var <- "Population" # variable of interest
  cats <- c() # other categorical variables in the dataset
  time_var <- "time"
  event_var <- "event"
  surv_tmp$event_factor <- factor(surv_tmp[[event_var]], levels = c(0, 1), labels = c("Censored", "Dead"))
  event_factor <- "event_factor" # event variable as a factor
  max_time <- max(surv_tmp[[time_var]])
  plot_names <- "dros2"
}

##################################################################################################################################
# describe the data
##################################################################################################################################

# f(t) with color for censoring separated by some variable of interest
f_t <- f_t_census(
  data = surv_tmp,
  time_var = time_var,
  event_var = event_factor,
  var = var,
  nbins = round(max_time / 5),
  ncol = 1
)
# save plot
# ggsave(plot = f_t, paste0("Figures/", plot_names, "_histograms.png"), width = 6.5, height = 4.5, dpi = 400)

# formatting is dataset-specific

if (dataset == "simu") {
  colors <- c("#F8766D", "#F8766D", "#00BFC4", "#00BFC4") # color by habitat
  linetypes <- c("solid", "dashed", "solid", "dashed") # linetype by mutation
} else if (dataset == "dros1") {
  colors <- NULL # automatic colours (one per curve)
  linetypes <- NULL # automatic linetypes (all solid)
} else if (dataset == "dros2") {
  colors <- NULL # automatic colours (one per curve)
  linetypes <- NULL # automatic linetypes (all solid)
}

# KM curves
km_curves_emp <- create_survival_curve(
  data = surv_tmp, title = "Empirical survival curves",
  group = var,
  time_var = time_var,
  colors = colors,
  linetypes = linetypes,
  event_var = event_var,
  xmax = max_time # x-axis max
  , legend = "right" # legend position
  , median_surv = FALSE # add median survival
)

# save plot
# ggsave(plot = km_curves_emp, paste0("Figures/", plot_names, "_KMcurves.png"), width = 6.5, height = 4.5, dpi = 400)

# KM curves - no confidence intervals (cleaner plos it useful later)
km_curves_emp_norib <- create_survival_curve(
  data = surv_tmp, title = "Empirical survival curves",
  group = var,
  time_var = time_var,
  colors = colors # manually colour lines
  , linetypes = linetypes # manually change line
  , event_var = event_var,
  conf_int = FALSE,
  xmax = max_time # x-axis max
  , legend = "right" # legend position
  , median_surv = FALSE # add median survival
)

# save plot
# ggsave(plot = km_curves_emp_norib, paste0("Figures/", plot_names, "_KMcurves_noribbon.png"), width = 6.5, height = 4.5, dpi = 400)

##################################################################################################################################
# cox
##################################################################################################################################

if (dataset == "simu") {
  coxmodel <- coxph(Surv(time, event) ~ mutation * habitat, data = surv_tmp)
  ph <- cox.zph(coxmodel)
  su <- summary(coxmodel)
} else if (dataset == "dros1") {
  coxmodel <- coxph(Surv(time, event) ~ temp, data = surv_tmp)
  ph <- cox.zph(coxmodel)
  su <- summary(coxmodel)
} else if (dataset == "dros2") {
  coxmodel <- coxph(Surv(time, event) ~ Population, data = surv_tmp)
  ph <- cox.zph(coxmodel)
  su <- summary(coxmodel)
}

ph # checking ph assumption
su # model's output

# plot the schoenfield residuals

if (dataset == "simu") {
  sch_mut <- plot(cox.zph(coxmodel), var = "mutation")
  sch_hab <- plot(cox.zph(coxmodel), var = "habitat")
  sch_int <- plot(cox.zph(coxmodel), var = "mutation:habitat")

  # save them all
  #  for (i in 1:3) {
  #    if (i == 1) var_sch <- "mutation"
  #    if (i == 2) var_sch <- "habitat"
  #    if (i == 3) var_sch <- "interaction"
  #    png(
  #      filename = paste0("Figures/", plot_names, "_SchRes_", var_sch, ".png"),
  #      width = 2000, height = 2000, res = 400
  #    )
  #    plot(cox.zph(coxmodel)[i])
  #    dev.off()
  #  }
} else if (dataset == "dros1") {
  sch_pop <- plot(cox.zph(coxmodel), var = "temp")

  # save
  #  var_sch <- "population"
  #  png(
  #    filename = paste0("Figures/", plot_names, "_SchRes_", var_sch, ".png"),
  #    width = 2000, height = 2000, res = 400
  #  )
  #  plot(cox.zph(coxmodel)[1])
  #  dev.off()
} else if (dataset == "dros2") {
  sch_pop <- plot(cox.zph(coxmodel), var = "Population")

  # save
  #  var_sch <- "population"
  #  png(
  #    filename = paste0("Figures/", plot_names, "_SchRes_", var_sch, ".png"),
  #    width = 2000, height = 2000, res = 400
  #  )
  #  plot(cox.zph(coxmodel)[1])
  #  dev.off()
}

##################################################################################################################################
# parametric
##################################################################################################################################

surv_par <- filter(surv_tmp, time > 0) # this package doesn't like survival times of 0

# we need to define n p and f for the dros datasets
if (dataset == "simu") {
  # n, p and f are alredy defined
} else if (dataset == "dros1") {
  n <- nrow(surv_par) # sample size
  p <- 1 - (sum(surv_par$event) / nrow(surv_par)) # proportion censored
  f <- max_time
} else if (dataset == "dros2") {
  n <- nrow(surv_par) # sample size
  p <- 1 - (sum(surv_par$event) / nrow(surv_par)) # proportion censored
  f <- max_time
}

# choose a distribution
par_fit <- compare_parametric_fits(
  data = surv_par # this custom function returns a big plot and a table for comparison
  , time_var = time_var,
  event_var = event_var,
  plot_title = paste0("Sample Size: ", n, ", Follow-up period: ", f, " Censored: ", p)
)

par_fit$comparison

# we'll choose the best fit that is not GenF or Genγ (they're not very parsimonious options...)

# fit it
if (dataset == "simu") {
  parmodel <- flexsurvreg(Surv(time, event) ~ mutation * habitat, data = surv_par, dist = "weibull")
} else if (dataset == "dros1") {
  parmodel <- flexsurvreg(Surv(time, event) ~ temp, data = surv_par, dist = "llogis")
} else if (dataset == "dros2") {
  parmodel <- flexsurvreg(Surv(time, event) ~ Population, data = surv_par, dist = "weibull")
}

tidy(parmodel)
plot(parmodel)

# save it
# png(
#  filename = paste0("Figures/", plot_names, "_parametric_fit_simpleton.png"),
#  width = 2000, height = 2000, res = 400
# )
plot(parmodel)
dev.off()


# right. Now with colours

if (dataset == "simu") {
  # first take out the model fits into a dataframe
  fits <- summary(parmodel, type = "survival")
  fit_df <- bind_rows(lapply(seq_along(fits), function(i) {
    data.frame(fits[[i]], strata = names(fits)[i])
  })) %>% # change the variable levels to match the km curves plot
    mutate(group = ifelse(strata == "mutation=no,habitat=Island B",
      "Island B_no", ifelse(strata == "mutation=no,habitat=Island A",
        "Island A_no", ifelse(strata == "mutation=yes,habitat=Island B",
          "Island B_yes", "Island A_yes"
        )
      )
    ))
  linetypes <- c("solid", "dashed", "solid", "dashed") # linetype by mutation
} else if (dataset == "dros1") {
  fits <- summary(parmodel, type = "survival")
  fit_df <- bind_rows(lapply(seq_along(fits), function(i) {
    data.frame(fits[[i]], strata = names(fits)[i])
  })) %>% # change the variable levels to match the km curves plot
    mutate(group = ifelse(strata == "temp=14°C",
      "14°C", ifelse(strata == "temp=18°C",
        "18°C", ifelse(strata == "temp=25°C",
          "25°C", ifelse(strata == "temp=30°C",
            "30°C", ifelse(strata == "temp=33°C",
              "33°C", "I got left out!"
            )
          )
        )
      )
    ))
  linetypes <- NULL
} else if (dataset == "dros2") {
  # first take out the model fits into a dataframe
  fits <- summary(parmodel, type = "survival")
  fit_df <- bind_rows(lapply(seq_along(fits), function(i) {
    data.frame(fits[[i]], strata = names(fits)[i])
  })) %>% # change the variable levels to match the km curves plot
    mutate(group = ifelse(strata == "Population=Zambia",
      "Zambia", ifelse(strata == "Population=South Africa",
        "South Africa", ifelse(strata == "Population=Austria",
          "Austria", "I got left out!"
        )
      )
    ))
  linetypes <- NULL
}

# now add them to our plot with the KM curves
plot <- km_curves_emp +
  geom_line(
    data = fit_df,
    aes(time, est, color = group, linetype = if (!is.null(linetypes)) group else NULL),
    linewidth = 1,
    show.legend = FALSE
  ) +
  geom_ribbon(
    data = fit_df,
    aes(x = time, ymin = lcl, ymax = ucl, fill = group),
    alpha = 0.1,
    show.legend = FALSE
  ) +
  theme(plot.title = element_text(size = 18)) +
  (if (!is.null(linetypes)) scale_linetype_manual(values = linetypes) else NULL) +
  labs(y = "Proportion alive", x = "Time", title = "Empirical survival curves and their parametric fits")
plot

# save plot
# ggsave(plot = plot, filename = paste0("Figures/", plot_names, "_KMcurves_with_parametric.png"), width = 6.5, height = 4.5, dpi = 400)


# right. Now without ribbons
plot_clean <- km_curves_emp_norib +
  geom_line(data = fit_df, aes(time, est, color = group, linetype = if (!is.null(linetypes)) group else NULL), size = 1, show.legend = FALSE) +
  labs(y = "Proportion alive", x = "Time", color = "Group", fill = "Group", title = "Empirical survival curves and their parametric fits")

plot_clean

# save plot
# ggsave(plot = plot_clean, filename = paste0("Figures/", plot_names, "_KMcurves_with_parametric_cleaner.png"), width = 6.5, height = 4.5, dpi = 400)

##################################################################################################################################
# taus
##################################################################################################################################

if (dataset == "simu") {
  categories <- c("mutation", "habitat")
} else if (dataset == "dros1") {
  categories <- NULL
} else if (dataset == "dros2") {
  categories <- NULL
}


# conditional survival matrix:
cond_surv <- cond_surv_mat(
  data = surv_tmp, # the dataset - a typical life table
  var = var, # variable of interest (string, must be categorical)
  cats = categories,
  time_var = time_var, # the time variable (string)
  event_var = event_var, # the event variable (string)
  res = 1, # resolution of the conditional survival grids for t and tau
  conf_int_level = 0.95, # confidence interval level (default 0.95)
  aggregate_by_cats = FALSE
)

# add power analysis here: assuming that A(t) is true, when did we get it?
# but this is only relevant for the simulation study...

if (dataset == "simu") {
  # let's check our ability to estimate A(t). First fit a km to the whole simu dataset
  km_fit_all <- survfit(Surv(span, event) ~ 1, data = surv) # km curve
  time_points <- seq(1, max(surv_tmp$time), by = 1) # useful to have in hand
  A_t_theo <- approx(x = km_fit_all$time, y = km_fit_all$surv, xout = time_points, yleft = 1)$y # estimate A(t) as S(t)
  A_t_theo <- A_t_theo / sum(A_t_theo, na.rm = TRUE) # normalise into PDF for comparison


  # KMs for subsets
  pa <- tibble(
    sample_size = numeric(),
    kl_div = numeric()
  )

  # now compare how well subsets can estimate A(t)
  for (i in seq(10, 2000, 1)) {
    # ss = i * nrow(surv_tmp) # sample size
    subset <- surv[sample(nrow(surv), i), ] # subset data
    km_fit_sub <- survfit(Surv(span, event) ~ 1, data = subset) # fit km
    A_t_empi <- approx(
      x = km_fit_sub$time,
      y = km_fit_sub$surv,
      xout = time_points,
      yleft = 1, rule = 2
    )$y # extend with 0s till the end
    A_t_empi <- A_t_empi / sum(A_t_empi, na.rm = TRUE) # PDF for comparison

    # compare with KL (information loss when approximating a distribution)
    loss <- kl_divergence(tru = A_t_theo, est = A_t_empi)
    pa <- bind_rows(pa, tibble(
      sample_size = i,
      kl_div = loss
    ))
  }


  # Create ggplot
  p <- ggplot(pa, aes(x = sample_size, y = kl_div)) +
    geom_line() +
    geom_point(size = 1) +
    # x axis log
    ylim(c(0, 2)) +
    labs(
      title = "Error in A(t) estimation as a function of n",
      x = "Sample Size (n)",
      y = "A(t) error (KL divergence)"
    ) +
    theme(
      axis.text = element_text(size = 11), # Font size for axis ticks
      axis.title = element_text(size = 12), # Adjust font of labels
      legend.text = element_text(size = 10),
      legend.title = element_text(size = 11),
      plot.margin = margin(8, 8, 8, 8), # Plot margins (t, r, b, l)
      plot.title = element_text(size = 16, hjust = 0.5) # Title settings
    )

  p
  # ggsave("Figures/PowerAnalysisAt.png", plot = p, width = 6.5, height = 4.5, dpi = 400)
}

# stats:
if (dataset == "simu") {
  tau_value <- 50
} else if (dataset == "dros1") {
  tau_value <- 125 # for these, do 43 nd 125
} else if (dataset == "dros2") {
  tau_value <- 53 # for these, do 31 and 53
}

stats <- pairwise_test(cond_surv, # the output of cond_surv_mat
  var_values = NULL, # which groups to show. If NULL, will show all groups
  tau_values = tau_value # can also be a vector of values, see below
)


# alternatively:

if (dataset == "simu") {
  var_values <- unique(cond_surv$cond_surv$unique_label) # specify which groups are of interest
  tau_values <- c(20, 20, 40, 40) # which tau is of interest for each group? (in the same order)
} else if (dataset == "dros1") {
  var_values <- unique(cond_surv$cond_surv$unique_label) # specify which groups are of interest
  tau_values <- c(125, 125, 125, 125, 125) # which tau is of interest for each group? (in the same order)
} else if (dataset == "dros2") {
  var_values <- unique(cond_surv$cond_surv$unique_label) # specify which groups are of interest
  tau_values <- c(53, 53, 53)
}

stats <- pairwise_test(cond_surv # the output of cond_surv_mat
  ,
  var_values = var_values # if !NULL, must be a vector of all possible values in the variable of interest
  , tau_values = tau_values # the age beyond which survival is interesting
)


# plotting with cond_surv_plot:

# if tau_value == NULL && collapse_time == FALSE (default values for both)
# heatmap of P(T>tau|T>t) for all t and tau
cond_surv_heatmap <- cond_surv_plot(cond_surv,
  ncol_hm = 2
)
# save plot
# ggsave(plot = cond_surv_heatmap, paste0("Figures/", plot_names, "_heatmaps.png"), width = 6.5, height = 4.5, dpi = 400)

# if tau_value is specified
# P(T>tau|T>t) ~ t for the specified tau value
# optional: select a specific group to plot. Default is NULL (all)

if (dataset == "simu") {
  var_value <- unique(cond_surv$cond_surv$unique_label) # must be one of unique(cond_surv$cond_surv$unique_label)
  tau_value <- 50 # choose a tau value
  color_var <- "habitat"
  linetype_var <- "mutation"
  markershape_var <- "mutation"
} else if (dataset == "dros1") {
  var_value <- unique(cond_surv$cond_surv$unique_label) # specify which groups are of interest
  tau_value <- 125 # for these, do 43 nd 125
  color_var <- NULL
  linetype_var <- NULL
  markershape_var <- NULL
} else if (dataset == "dros2") {
  var_value <- unique(cond_surv$cond_surv$unique_label)
  tau_value <- 53 # for this one, do 31 and 53
  color_var <- NULL
  linetype_var <- NULL
  markershape_var <- NULL
}

cond_surv_line <- cond_surv_plot(cond_surv,
  var_value = var_value,
  tau_value = tau_value,
  color_var = color_var,
  linetype_var = linetype_var
)
# save plot
# ggsave(plot = cond_surv_line, paste0("Figures/", plot_names, "_P(T>", tau_value, "|T>t).png"), width = 6.5, height = 4.5, dpi = 400)

# if tau_value == NULL && collapse_time == TRUE
# P(T>tau) ~ tau
cond_surv_curve <- cond_surv_plot(cond_surv,
  collapse_time = TRUE,
  color_var = color_var,
  linetype_var = linetype_var
)

# save plot
# ggsave(plot = cond_surv_curve, paste0("Figures/", plot_names, "_P(T>tau).png"), width = 6.5, height = 4.5, dpi = 400)


# if tau_value != NULL && collapse_time == TRUE
# P(T>tau) when age is unknown
if (dataset == "simu") {
  tau_value <- 50
} else if (dataset == "dros1") {
  tau_value <- 125 # for this one, do 43 and 125
} else if (dataset == "dros2") {
  tau_value <- 53 # for this one, do 31 and 53
}


cond_surv_scatter <- cond_surv_plot(cond_surv,
  tau_value = tau_value,
  collapse_time = TRUE,
  color_var = color_var,
  markershape_var = markershape_var
)
# ggsave(plot = cond_surv_scatter, paste0("Figures/", plot_names, "_P(T>", tau_value, ").png"), width = 6.5, height = 4.5, dpi = 400)

# if each group has a different tau of interest (only for simu):

if (dataset == "simu") {
  var_value <- unique(cond_surv$cond_surv$unique_label) # must be within unique(cond_surv$cond_surv$unique_label)
  tau_value <- c(20, 20, 40, 40) # same order (only really applicatble to the simu case)

  cond_surv_scatter2 <- cond_surv_plot(cond_surv,
    var_value = var_value,
    tau_value = tau_value,
    collapse_time = TRUE,
    color_var = color_var,
    markershape_var = markershape_var
  )

  # ggsave(plot = cond_surv_scatter2, paste0("Figures/", plot_names, "_P(T>taus).png"), width = 6.5, height = 4.5, dpi = 400)
}
