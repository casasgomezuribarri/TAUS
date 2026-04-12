# TAUS: Target-Age Unified Survival (R package)

**TAUS** computes conditional survival P(T>τ|T>t) and the probability O_τ that a randomly sampled individual of unknown age will outlive a target age τ. It does not assume proportional hazards or a parametric survival form. Pairwise comparison of O_τ across groups is supported via a Beta–Kolmogorov–Smirnov bootstrap.

## Installation

From source (from the repo root, with R ≥ 4.0):

```r
install.packages("remotes")
remotes::install_github("casasgomezuribarri/TAUS")
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

None

---
