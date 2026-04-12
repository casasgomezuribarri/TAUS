# TAUS: Target-Age Unified Survival (R package)

**TAUS** computes conditional survival P(T>τ|T>t) and the probability O_τ that a randomly sampled individual of unknown age will outlive a target age τ. It does not assume proportional hazards or a parametric survival form. Pairwise comparison of O_τ across groups is supported via a Beta–Kolmogorov–Smirnov bootstrap.

## Installation

From source (from the repo root, with R ≥ 4.0):

```r
install.packages("remotes")
remotes::install_github("casasgomezuribarri/Survival", subdir = ".", ref = "package")
```

Or clone the repo, then in R:

```r
install.packages(".", repos = NULL, type = "source")
```

## Quick example

```r
library(TAUS)
library(dplyr)
# Minimal survival data: time, event, group
dat <- tibble(
  time  = c(2, 4, 6, 5, 3, 5, 4, 6),
  event = c(1, 1, 1, 0, 1, 1, 0, 0),
  group = c(rep("A", 4), rep("B", 4))
)
# Conditional survival matrix and O_tau
out <- cond_surv_mat(dat, var = "group", time_var = "time", event_var = "event", res = 2)
# Pairwise comparison at tau = 4
pairwise_test(out, var_values = c("A", "B"), tau_values = c(4, 4), B = 500)
```

## Main functions

| Function | Description |
|----------|-------------|
| `cond_surv_mat()` | Conditional survival matrix and O_τ with confidence intervals |
| `estimate_beta_params()` | Beta shape parameters from mean and CI limits |
| `kl_divergence()` | Kullback–Leibler divergence between two probability vectors |
| `pairwise_test()` | Pairwise comparison of O_τ across groups (Beta–KS bootstrap) |

## Tests

```r
devtools::test()   # or: testthat::test_local(".")
```

## Authors

Iván Casas Gomez-Uribarri (maintainer), Simon A. Babayan, Fredros Okumu, Francesco Baldini, Mauro Pazmiño Betancourth.

Methodology: see the TAUS manuscript (this repository and associated thesis/paper materials).

## Licence

GPL-3.

---

# This repository is part of a group
This repository belongs to a group of repositories with all my PhD work (data, code, figures, reports...). All of them are private, but they might become partially public as needed (e.g. for revision, collaboration or publication). The general structure is as follows:

```
.
├── Eggs
├── EIP
├── MIRS
├── Oocysts
├── qPCR
├── Quarto
└── Survival  # <- you are here
```

Each repository deals with with a part of the project, but they're all related.

Once all the work is done, these may be restructured and published as one.

# Current repo: Survival

This repo analyses the effect of exposure to _P. falciparum_ on mosquito survival considering the following covariates: mean temp, temp range and mosquito species.

It also develops a whole new model for survival analysis, and tests it with simulated and real-world data

#### Repo structure:

```
|Survival
│ └── Code
│   └── age.r
│   └── compile_data_survival.r
│   └── conditional_survival_simulation.r
│   └── conditional_survival.r
│   └── considering_causality.r
│   └── functions.r
│   └── set_up_env.r
│   └── survival_analysis.r
│   └── survival_cox_weighted.r
│   └── survival_cox.r
│   └── survival_KM.r
│   └── survival_parametric.r
│   └── survival_taus_old.r
│   └── survival_taus.r
│   └── survival_timeless.r
│ └── Data
│   └── online
│     └── dros1.csv
│     └── dros2.csv
│   └── synth
│     └── simu.csv
│   └── surv_pot_summary.csv
│   └── surv.xlsx
│   └── survival.csv
│   └── survival.xlsx
│ └── Figures 
│   └── ...
│ └── logs 
│   └── ...
│ └── models 
│   └── ...
│ └── TAUS
│   └── comparing_models.r
│   └── simulate_data.r
│   └── taus.r
```

The **Code** folder contains the scripts for the analysis of mosquito survival.

- `age.r` - compares the ages at the time of exposure, showing they are not significantly different across replicates.
- `compile_data_survival.r` - fecthes data from Excel and outputs it in the right format. Shouldn't need to be opened because it is run from other scripts
- `conditional_survival_simulation.r` - ...
- `conditional_survival.r` - ...
- `considering_causality.r` - calculates propensity scores, draws dags, and other useful tools for consieration of causality in the analysis
- `functions.r` - defines useful custom functions
- `set_up_env.r` - defines a useful function for setting up environment right
- `survival_analysis.r` - survival analysis - here the results of the scripts below will be eventually joined for a nice story.
- `survival_cox_weighted.r` - coxphw analysis (includes augmented backward elimination)
- `survival_cox.r` - coxph analysis (includes tests for PH assumption)
- `survival_KM.r` - KM (descriptive) survival analysis
- `survival_parametric.r` - parametric survival analysis (does this really make sense?)
- `survival_taus_old.r` - archived code just in case - not in use for now (I think)
- `survival_taus.r` - taus survival analysis
- `survival_timeless.r` - ...

The **Data** folder contains the data files (some csv files are created by `compile_data_survival.r`)
- online
  - `dros1.csv`
  - `dros2.csv`
- synth
  - `simu.csv`
- `surv_pot_summary.csv`
- `surv.csv`
- `survival.csv`
- `survival.xlsx`

The **Figures** folder stores all plots generated during this analysis. Not all of them will be published, but they were all useful.

The **logs** folder stores all logs generated during augmented backaward elimination (currently empty) - deprecated method (ignore)

The **models** folder stores all models generated during augmented backaward elimination - deprecated method (ignore)

The **TAUS** folder contains scripts for the TAUS manuscript (simulations, comparisons). The installable R package lives in the repo root: `R/`, `man/`, `tests/`, and `DESCRIPTION`. See the TAUS package section at the top of this readme for installation and usage.
