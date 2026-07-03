# Component-additive multilevel network meta-regression (ML-NMR)

The Bayesian flagship of cpaic. The relative effect of every treatment
is the sum of its component effects (`theta = C beta`), estimated
jointly from individual patient data (IPD) and aggregate data (AgD).
Aggregate arms are fitted by integrating the individual-level model over
each study's covariate distribution. Because disconnected sub-networks
share component parameters, the network is connected by construction.

## Usage

``` r
cmlnmr(
  ipd,
  agd,
  effect_modifiers,
  inactive = NULL,
  sep.comps = "+",
  family = "binomial",
  study = ".study",
  trt = ".trt",
  outcome = ".y",
  time = ".time",
  exposure = ".exposure",
  r = "r",
  n = "n",
  E = "E",
  se = "se",
  cut_points = NULL,
  interval = ".interval",
  baseline = c("piecewise", "mspline"),
  n_basis = 6L,
  cor = NULL,
  n_int = 64L,
  prior_intercept_sd = 10,
  prior_beta_sd = 10,
  prior_reg_sd = 2.5,
  chains = 4L,
  iter_warmup = 500L,
  iter_sampling = 500L,
  seed = NULL,
  ...
)
```

## Arguments

- ipd:

  Individual patient data (one row per patient).

- agd:

  Aggregate data (one row per arm) with the per-study covariate
  summaries `x_mean`, `x_sd` for each effect modifier `x`.

- effect_modifiers:

  Character vector of effect-modifier names.

- inactive, sep.comps:

  Component coding (see
  [`cpaic_network()`](https://choxos.github.io/cpaic/reference/cpaic_network.md)).

- family:

  One of `"binomial"`, `"gaussian"`, `"poisson"`, `"survival"`.

- study, trt:

  Column names (in both `ipd` and `agd`).

- outcome:

  IPD outcome column: 0/1 (binomial), numeric (gaussian), count
  (poisson), or the event indicator for survival.

- time, exposure:

  IPD time / exposure column (survival, poisson).

- r, n, E, se:

  Aggregate columns: events `r`, sample size `n` (binomial), exposure
  `E` (poisson/survival), mean `outcome` and its standard error `se`
  (gaussian).

- cut_points:

  Survival only: interior interval boundaries for the
  piecewise-exponential baseline. `NULL` (default) gives the exponential
  model; e.g. `c(6, 12)` gives three intervals.

- interval:

  Survival only: name of the aggregate interval-index column (values
  `1..K`) required when `cut_points` are supplied.

- baseline:

  Survival baseline hazard: `"piecewise"` (default, step function) or
  `"mspline"` (smooth M-spline over the `cut_points` grid).

- n_basis:

  Number of M-spline basis functions when `baseline = "mspline"` (must
  be `<=` the number of intervals).

- cor:

  Optional covariate correlation matrix for the Gaussian-copula
  integration. Defaults to the IPD correlation; pass a matrix to
  override or the identity to integrate covariates independently.

- n_int:

  Integration points per aggregate arm (ignored for `gaussian`, which is
  exact at the covariate means).

- prior_intercept_sd, prior_beta_sd, prior_reg_sd:

  Prior SDs.

- chains, iter_warmup, iter_sampling, seed:

  Passed to `cmdstanr`.

- ...:

  Reserved.

## Value

An object of class `cpaic_mlnmr` with the `cmdstanr` fit, the component
design, and a tidy table of component effects.

## Details

Supported families: `"binomial"` (logit), `"gaussian"` (identity),
`"poisson"` (log), and `"survival"`. Survival uses a
proportional-hazards model with a flexible baseline:
`baseline = "piecewise"` gives a piecewise-exponential step baseline
(one level per interval defined by `cut_points`; exponential when
`cut_points = NULL`), while `baseline = "mspline"` gives a smooth
baseline hazard from a non-negative M-spline basis (`n_basis` functions,
simplex weights) evaluated at the interval midpoints over the
`cut_points` grid. Individual patient data are split at the cut points
internally; aggregate survival data are supplied as events and
person-time per arm and (when `cut_points` are given) per interval. The
aggregate likelihood approximates the expected events in an arm-interval
by person-time times the integrated hazard; this assumes the person-time
is independent of the covariates within an interval (most accurate when
the intervals are narrow), so it is approximate when effect modifiers
also drive censoring.

Aggregate covariates are integrated with a Gaussian copula: the
correlation is estimated from the IPD (or supplied via `cor`) so the
integration points respect the covariate correlation structure rather
than treating the effect modifiers as independent.

Identifiability note: a component whose effect-modifier interaction is
informed only by aggregate arms is weakly identified, because main and
interaction effects are constrained only through the integrated
population-average outcome. Supply IPD for such components where
possible; `prior_reg_sd` regularizes otherwise.

## See also

[`cmaic()`](https://choxos.github.io/cpaic/reference/cmaic.md),
[`cstc()`](https://choxos.github.io/cpaic/reference/cstc.md),
[`cnma_bridge()`](https://choxos.github.io/cpaic/reference/cnma_bridge.md)

## Examples

``` r
if (FALSE) { # requireNamespace("cmdstanr", quietly = TRUE) && !inherits(try(cmdstanr::cmdstan_path(), silent = TRUE), "try-error")
# \donttest{
ipd <- data.frame(.study = "S1",
                  .trt = rep(c("Placebo", "A"), each = 100),
                  .y = rbinom(200, 1, 0.5), x1 = rnorm(200))
agd <- data.frame(.study = "S2", .trt = c("Placebo", "A+B"),
                  r = c(40, 55), n = c(100, 100),
                  x1_mean = c(0.2, 0.2), x1_sd = c(1, 1))
fit <- cmlnmr(ipd, agd, effect_modifiers = "x1", inactive = "Placebo",
              chains = 2, iter_warmup = 200, iter_sampling = 200)
component_effects(fit)
# }
}
```
