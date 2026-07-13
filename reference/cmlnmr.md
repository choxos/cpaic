# Component-additive multilevel network meta-regression (ML-NMR)

The Bayesian flagship of cpaic. The relative effect of every treatment
is the sum of its component effects, estimated jointly from individual
patient data (IPD) and aggregate data (AgD). Aggregate arms are fitted
by integrating the individual-level model over each study's covariate
distribution, averaging the outcome on its natural scale (not the link
scale). Because disconnected sub-networks share component parameters,
the network is connected by construction.

## Usage

``` r
cmlnmr(
  ipd,
  agd,
  effect_modifiers,
  inactive = NULL,
  sep.comps = "+",
  family = "binomial",
  margins = NULL,
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
  summaries `x_mean` (and `x_sd` for normal margins) for each effect
  modifier `x`.

- effect_modifiers:

  Character vector of effect-modifier names.

- inactive, sep.comps:

  Component coding (see
  [`cpaic_network()`](https://choxos.github.io/cpaic/reference/cpaic_network.md)).
  `inactive = NULL` gives the *unanchored* component parameterization,
  in which every unit receives its own parameter (Wigle & Béliveau
  2022).

- family:

  One of `"binomial"`, `"gaussian"`, `"poisson"`, `"survival"`.

- margins:

  Optional named character vector giving the integration margin of each
  effect modifier: `"normal"` or `"bernoulli"`. Defaults to Bernoulli
  for 0/1 covariates and normal otherwise.

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

  Survival only: interior interval boundaries for the piecewise-constant
  baseline. `NULL` (default) gives the exponential model; e.g.
  `c(6, 12)` gives three intervals.

- interval:

  Survival only: name of the aggregate interval-index column (values
  `1..K`) required when `cut_points` are supplied.

- baseline:

  Survival baseline hazard: `"piecewise"` (default, free step heights)
  or `"mspline"` (step heights smoothed by an M-spline evaluated at
  interval midpoints; see the Survival section).

- n_basis:

  Number of M-spline basis functions when `baseline = "mspline"` (must
  be `<=` the number of intervals).

- cor:

  Optional covariate correlation matrix for the Gaussian-copula
  integration. Must be a positive-definite correlation matrix (unit
  diagonal). Defaults to the within-study IPD correlation.

- n_int:

  Integration points per aggregate arm (ignored for `gaussian`, which is
  exact at the covariate means).

- prior_intercept_sd, prior_beta_sd, prior_reg_sd:

  Prior SDs.

- chains, iter_warmup, iter_sampling, seed:

  Passed to `cmdstanr`.

- ...:

  Passed to the `cmdstanr` sampler (e.g. `adapt_delta`).

## Value

An object of class `cpaic_mlnmr` with the `cmdstanr` fit, the component
design, and a tidy table of component effects.

## Details

The model includes component x effect-modifier interactions `gamma`, so
the treatment effect is **population-specific**: \$\$\theta_t(x) = C_t'
(\beta + \Gamma x).\$\$ The component main effects `beta` are the
effects at the covariate origin (`x = 0`) and are *not* by themselves a
population-adjusted quantity. Use `newdata` in
[`relative_effects()`](https://choxos.github.io/cpaic/reference/relative_effects.md)
/
[`component_effects()`](https://choxos.github.io/cpaic/reference/component_effects.md)
to obtain effects in a named target population.

Supported families: `"binomial"` (logit), `"gaussian"` (identity),
`"poisson"` (log), and `"survival"`.

## Integration

Aggregate covariates are integrated with Sobol' quasi-Monte-Carlo points
coupled by a Gaussian copula, whose correlation is pooled *within* IPD
studies on the Fisher z scale (or supplied via `cor`). Each covariate is
pushed through its own marginal inverse CDF: `margins` may be `"normal"`
(using `x_mean` and `x_sd`) or `"bernoulli"` (using `x_mean` as the
prevalence). Margins default to Bernoulli for covariates that are 0/1 in
the IPD and normal otherwise; a normal margin on a binary covariate
would integrate over a population that cannot occur.

## Survival (approximations, read before use)

Survival uses a proportional-hazards model with a piecewise-constant
baseline log-hazard on the `cut_points` grid (`cut_points = NULL` gives
the exponential model). With `baseline = "mspline"` the interval heights
are *smoothed* by an M-spline basis evaluated at the interval midpoints.
This is a piecewise-exponential model with a smoothed step baseline; it
is **not** the continuous-time integrated M-spline survival likelihood
of `multinma`, which uses both the M-spline hazard basis and its
integrated (I-spline) cumulative hazard and supports left-,
interval-censoring and delayed entry. cpaic handles right-censoring
only.

Aggregate survival data are supplied as events and person-time per arm
and per interval, and the expected events are approximated by
`person-time x mean hazard`. This is an approximation even without
censoring, because higher-hazard individuals leave the risk set earlier,
so the *baseline* covariate distribution does not describe the covariate
distribution of the accumulated person-time. In a two-group example
(hazards 0.1 and 0.4, 50:50, follow-up to t = 10) it overstates expected
events by about 36%. Narrow intervals reduce the bias; supply IPD, or
interval-specific risk-set covariate summaries, where accuracy matters.

## Identifiability

A relative effect is uniquely estimable only if its component contrast
lies in the row space of the within-study component design (Wigle et al.
2026);
[`relative_effects()`](https://choxos.github.io/cpaic/reference/relative_effects.md)
returns `NA` otherwise rather than a prior-driven number. Note this
checks identification of `beta`; a component x effect-modifier
interaction is additionally identified only by covariate variation on
the contrasts that involve it, and interactions informed only by
aggregate arms are weakly identified (`prior_reg_sd` regularizes).

## References

Phillippo DM, Dias S, Ades AE, et al. (2020). Multilevel network
meta-regression for population-adjusted treatment comparisons. *JRSS A*,
183(3), 1189–1210.

Wigle A, Beliveau A, Nikolakopoulou A, Lin L (2026). Creating Treatment
and Component Hierarchies in Component Network Meta-Analysis.

## See also

[`cmaic()`](https://choxos.github.io/cpaic/reference/cmaic.md),
[`cstc()`](https://choxos.github.io/cpaic/reference/cstc.md),
[`cnma_bridge()`](https://choxos.github.io/cpaic/reference/cnma_bridge.md),
[`estimable_effects()`](https://choxos.github.io/cpaic/reference/estimable_effects.md)

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
# Effects in a named target population (x1 = 0.2), not at the origin:
relative_effects(fit, newdata = data.frame(x1 = 0.2))
# }
}
```
