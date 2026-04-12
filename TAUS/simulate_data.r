# this script creates a dataset of individuals on two islands, A and B, with different traits and lifespans.

# a few packages 
packages = c("tidyverse"
             , "tibble"
             , "dplyr"
             , "rstudioapi"
)

for (i in packages) {
  if (!require(i, character.only = TRUE)) install.packages(i)
  library(i, character.only = TRUE)
}

# a couple useful fucntions for setting up the right working directory:

# returns full path to parent folder of current script
whereami <-  function()
{
  this_file <- commandArgs() %>% 
    tibble::enframe(name = NULL) %>%
    tidyr::separate(col=value, into=c("key", "value"), sep="=", fill='right') %>%
    dplyr::filter(key == "--file") %>%
    dplyr::pull(value)
  if (length(this_file)==0)
  {
    this_file <- rstudioapi::getSourceEditorContext()$path
  }
  return(dirname(this_file))
}

# set working directory to be the parent of the above
setwd_grandparent <- function()
{
  this_script <- whereami()
  setwd(dirname(this_script))
  # print working directory
  print("Current working directory:")
  print(getwd())
}

# apply it
setwd_grandparent() # working directory should be the folder of where TAUS is

source("TAUS/taus.r") # loads dependencies and custom functions for the TAUS project

# load packages
packages <- c(
    "FamilyRank",
    "flexsurv"
)

for (i in packages) {
    if (!require(i, character.only = TRUE)) install.packages(i)
    library(i, character.only = TRUE)
}



# parameters
initial_pop <- 1000
max_time <- 10000

# initialise
individuals <- list()
id_counter <- 1
set.seed(1999)

# simulate data
for (island in c("Island A", "Island B")) { # in each island
    for (t in 1:max_time) { # each day there are a few births
        if (t == 1) {
            n_births <- initial_pop # hard-coded for the first day in both islands
        } else if (island == "Island A") {
            n_births <- rpois(1, lambda = 7) # otherwise, a poisson process
        } else if (island == "Island B") {
            n_births <- rpois(1, lambda = 12) # there's chocolate on island b
        }

        for (i in 1:n_births) { # for each newborn...
            birth_time <- t # record birthtime
            habitat <- island # record island population
            mutation <- sample(c("yes", "no"), 1) # random mutation

            # compute lifespan - first, a baseline lifespan assigned by gOD:
            # a bimodal distribution wth some chance of early mortality - island specific
            if (island == "Island A") { # high infant mortality on island A
                for_gODS_sake <- rbinorm(n = 1, mean1 = 10, mean2 = 60, sd1 = 3, sd2 = 20, prop = .34)
            } else if (island == "Island B") { # low infant mortality on island B
                for_gODS_sake <- rbinorm(n = 1, mean1 = 10, mean2 = 60, sd1 = 3, sd2 = 20, prop = .18)
            }

            # can check below what for_gODS_sake looks like
            # hist(rbinorm(n=1000, mean1=10, mean2=60, sd1=3, sd2=20, prop=.1), breaks = 80)

            # force it to be positive
            base_lifespan <- max(for_gODS_sake, 1)

            if (base_lifespan < 15) { # for those with low lifespan, the mutation is good
                lifespan <- base_lifespan +
                    ifelse(mutation == "yes" # if mutated
                        , rnorm(1, 10, 5), 0
                    ) + # extend lifespan by rnorm(1, 10, 5)
                    rnorm(1, 0, 1) # plus individual-level noise
            } else if (base_lifespan > 25) { # for those with high lifespan
                if (island == "Island A") {
                    lifespan <- base_lifespan + # mutation is good on island A
                        ifelse(mutation == "yes" # if mutated:
                            , rnorm(1, 35, 20), 0
                        ) +
                        rnorm(1, 0, 1) # plus individual-level noise
                } else if (island == "Island B") {
                    lifespan <- base_lifespan + # mutation is bad on island B
                        ifelse(mutation == "yes" # if mutated:
                            , rnorm(1, -15, 5), 0
                        ) +
                        rnorm(1, 0, 1) # plus individual-level noise
                }
            } else {
                lifespan <- base_lifespan + rnorm(1, 0, 1) # individual-level noise
            }

            death_time <- birth_time + max(lifespan, 0.5) # ensure lifespan > 0

            # record in dataset
            individuals[[id_counter]] <- data.frame(
                id = id_counter,
                birth_time = birth_time,
                death_time = death_time,
                habitat = island,
                mutation = mutation
            )

            id_counter <- id_counter + 1
        }
    }
}
# combine all individuals into a single data frame
pop_df <- do.call(rbind, individuals)

# polish the dataframe
simu_surv <- pop_df %>%
    mutate(
        span = round(death_time) - birth_time, # calculate lifespan
        event = sample(1, size = n(), replace = TRUE) # all deaths observed
    ) %>%
    mutate(
        time = if_else(event == 0, round(span * runif(n())), span) # if censored, span observed is a random fraction of the total lifespan
    )

# survival curves
km <- survfit(Surv(span, event) ~ habitat + mutation, data = simu_surv)
plot(km, conf.int = TRUE)

# try parametric fits
simu_surv <- simu_surv %>% filter(time > 0) # parametric survival package doesn't like survival times of 0

par_fits <- compare_parametric_fits( # a custom function to visually compare all fits. See TAUS/taus.r
    data = simu_surv,
    time_var = "time",
    event_var = "event"
)


# try cox fit
full <- coxph(Surv(time, event) ~ mutation * habitat, data = simu_surv)
summary(full)

cox.zph(full) # if p<0.05, PH assumption is violated and cox model is not valid

# save
write.csv(pop_df, file = "Data/synth/simu_toy.csv", row.names = FALSE)
