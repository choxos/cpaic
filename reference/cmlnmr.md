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
  start = ".start",
  entry = ".entry",
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
  QR = FALSE,
  trt_effects = c("fixed", "random"),
  re_parameterization = c("noncentered", "centered"),
  prior_intercept_sd = 2.5,
  prior_aux_sd = 1,
  prior_beta_sd = 2.5,
  prior_sigma_sd = 2.5,
  prior_reg_sd = 1,
  prior_gamma_dist = c("normal", "student_t"),
  prior_gamma_scale = 1,
  prior_gamma_df = 4,
  prior_tau_dist = c("half-normal", "half-student-t"),
  prior_tau_scale = 1,
  prior_tau_df = 4,
  prior_predictive = FALSE,
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
  effect modifier: `"normal"`, `"bernoulli"`, `"gamma"`, `"lognormal"`,
  or `"beta"`. The last three are set from the study mean and SD by
  method of moments (`gamma`/`lognormal` need a positive mean; `beta`
  needs a mean in `(0, 1)` and `sd^2 < mean(1 - mean)`). Defaults to
  Bernoulli for 0/1 covariates and normal otherwise.

- study, trt:

  Column names (in both `ipd` and `agd`).

- outcome:

  IPD outcome column: 0/1 (binomial), numeric (gaussian), count
  (poisson), or the event indicator for survival.

- time, exposure:

  Outcome-time column for survival in both IPD and AgD; IPD exposure
  column for Poisson outcomes.

- start, entry:

  Survival columns giving the lower endpoint for interval-censored
  outcomes and the delayed-entry time. Missing columns imply zero.

- r, n, E, se:

  Aggregate columns: events `r`, sample size `n` (binomial), exposure
  `E` (poisson), mean `outcome` and its standard error `se` (gaussian).

- cut_points:

  Survival only: interior interval boundaries for a piecewise baseline.
  `NULL` gives the exponential model. This argument is ignored for a
  continuous M-spline baseline.

- interval:

  Retained for source compatibility; exact survival data do not use
  interval-indexed event counts.

- baseline:

  Survival baseline hazard: `"piecewise"` (default, free step heights)
  or `"mspline"` (a continuous cubic M-spline with its exact integrated
  basis).

- n_basis:

  Number of cubic M-spline basis functions. Must be at least 4.

- cor:

  Optional covariate correlation matrix for the Gaussian-copula
  integration. Must be a positive-definite correlation matrix (unit
  diagonal). Defaults to the within-study IPD correlation. For gamma,
  lognormal, or beta margins the auto-estimated correlation is only an
  approximation to the latent copula correlation; supply `cor` on the
  latent scale to control it exactly.

- n_int:

  Integration points per aggregate arm (ignored for `gaussian`, which is
  exact at the covariate means).

  This is the main cost lever for the survival families. An aggregate
  survival arm is supplied as reconstructed pseudo-IPD, so the aggregate
  likelihood is evaluated once per (aggregate row x integration point):
  the work grows as `nrow(agd) * n_int`, and the default of 64 is
  expensive on a trial with several hundred reconstructed patients.
  Sampling is usually well behaved on the fixed-effects model (no
  divergences in the fixed-effects checks here), though the
  random-effects survival model can still produce a few divergent
  transitions and occasional rejected simplex proposals; inspect the
  diagnostics rather than assuming they are clean. If a survival fit is
  slow, reduce `n_int` before suspecting the geometry, and confirm the
  answer is stable with
  [`plot_integration_error()`](https://choxos.github.io/cpaic/reference/plot_integration_error.md).

- QR:

  Logical scalar. If `TRUE`, apply the scaled thin QR reparameterization
  used by `multinma` to the complete fixed-effects design matrix. This
  is only a reparameterization: it must not change the posterior
  distribution, only the geometry the sampler explores. The default is
  `FALSE`, matching `multinma`.

  Do not turn this on expecting a free improvement. On the component
  networks tested here the fixed-effects design was not badly
  conditioned (a condition number near 19, in a network where every
  active treatment shared a component), and `QR = TRUE` gave *fewer*
  effective samples per second than `QR = FALSE`, with no divergent
  transitions either way. The intuition that a component design must be
  severely collinear, because one component recurs across many
  multi-component treatments, is not borne out: the study intercepts and
  the spread of the integration points keep the conditioning mild. Check
  `Z_cond` on the fit, and reach for `QR = TRUE` only when it is large.

- trt_effects:

  Treatment-effect model: `"fixed"` or `"random"`.

- re_parameterization:

  Random-effects parameterization. The default `"noncentered"` should be
  used for inference; `"centered"` is provided for sampling diagnostics.

- prior_intercept_sd, prior_beta_sd, prior_reg_sd:

  Standard deviations for study-intercept, component-effect, and
  prognostic-regression normal priors.

- prior_aux_sd:

  Scale of the half-normal prior on the baseline-hazard smoothing
  parameter (survival families only). Each study has its own baseline
  hazard, given a first-order random-walk prior on the log spline
  coefficients with this shared smoothing scale. This is a simplified
  relative of the smoothing prior in `multinma`, not the same prior:
  smaller values shrink every study's baseline toward equal spline
  weights, which is a smooth default shape and not a constant hazard.
  The default of 1 follows the Stan recommendation of a
  half-normal(0, 1) prior for a hierarchical scale.

- prior_sigma_sd:

  Scale of the half-normal prior on the Gaussian residual standard
  deviation (gaussian family only), kept separate from `prior_beta_sd`
  so the treatment-effect prior and the residual-noise prior are
  independent.

- prior_gamma_dist, prior_gamma_scale, prior_gamma_df:

  Distribution, scale, and degrees of freedom for interaction priors.
  The Student t option uses the stated degrees of freedom.

- prior_tau_dist, prior_tau_scale, prior_tau_df:

  Distribution, scale, and degrees of freedom for the positive
  heterogeneity prior.

- prior_predictive:

  If `TRUE`, sample from the prior and omit the observed likelihood.
  Replicated outcomes remain available for
  [`prior_predictive_check()`](https://choxos.github.io/cpaic/reference/prior_predictive_check.md).

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

## Random effects

`trt_effects = "random"` adds study-arm deviations around the
component-implied relative effects. Deviations use a non-centered
parameterization by default. Within a multi-arm study, deviations
relative to the study baseline have the standard NMA correlation of 0.5.
The heterogeneity standard deviation `tau` has a half-normal(0, 1) prior
by default. The centered parameterization is available only to reproduce
sampling comparisons.

## Priors

Defaults follow the Stan prior-choice recommendations. Component effects
use normal(0, 2.5), component by effect-modifier interactions use
normal(0, 1), study intercepts use normal(0, 2.5), and `tau` uses
half-normal(0, 1). Interaction priors do real regularization when Gamma
is weakly identified, so every fitted object records the complete prior
specification. Use
[`prior_sensitivity()`](https://choxos.github.io/cpaic/reference/prior_sensitivity.md)
to quantify contrast movement and `prior_predictive = TRUE` with
[`prior_predictive_check()`](https://choxos.github.io/cpaic/reference/prior_predictive_check.md)
to inspect prior implications before fitting the likelihood.

## Survival

Survival uses the exact individual likelihood ported from `multinma`
(Phillippo et al. 2020). The model evaluates a hazard basis and its
integrated cumulative-hazard basis at every outcome, interval start, and
delayed-entry time. It supports observed events, right censoring, left
censoring, interval censoring, and delayed entry.
`baseline = "piecewise"` gives a piecewise-exponential baseline;
`baseline = "mspline"` gives a continuous cubic M-spline baseline.

Aggregate survival input must contain reconstructed event and censoring
rows with the same outcome-time columns as IPD, plus repeated arm-level
covariate summaries. The likelihood of every aggregate row is averaged
over its covariate integration points with `log_sum_exp`. Aggregate
event counts and person-time alone cannot recover this likelihood and
are rejected explicitly.

Two qualifications, so that "exact" is not read more broadly than it
should be.

- **The likelihood is exact; the covariate integration is not.** Every
  individual contribution (event, right, left and interval censoring,
  delayed entry) is the exact analytic expression, verified against
  closed form to machine precision. The *aggregate* likelihood, however,
  averages that exact contribution over a finite quasi-Monte-Carlo grid
  of `n_int` covariate points, so it carries an integration error that
  shrinks with `n_int` but is not zero. Increase `n_int` and confirm
  that the estimates are stable before relying on them. The earlier
  person-time approximation was biased by 36% in a two-group example;
  that particular bias is removed, but a finite integration error
  remains.

- **Each study has its own baseline hazard shape.** Every study carries
  its own set of spline (or step) coefficients, smoothed toward a common
  shape by a shared random-walk scale (`prior_aux_sd`), so the treatment
  effects do not have to absorb baseline misfit. A single global spline
  basis is built from the pooled follow-up range, so a study with much
  shorter follow-up may not inform the coefficients of the latest basis
  functions; those are then determined by the smoothing prior rather
  than by that study's data.

## Scope and current limitations

Two gaps are worth naming for anyone comparing this with `multinma`.

- **Effects are reported as conditional contrasts at a covariate
  value**, `(C_t - C_u)'(beta + Gamma x)`, on the linear-predictor
  scale.
  [`relative_effects()`](https://choxos.github.io/cpaic/reference/relative_effects.md)
  evaluates this at the target in `newdata`. There is no marginal
  (population-standardized) effect path yet: on a non-collapsible scale
  the conditional effect at a point differs from the average effect over
  a population with a distribution of covariates, and only the former is
  returned.

- **Every effect modifier enters both the prognostic terms and the full
  set of component interactions.** There is no prognostic-only covariate
  role (unlike
  [`cstc()`](https://choxos.github.io/cpaic/reference/cstc.md), which
  separates `prognostics`), so a covariate that shifts outcomes without
  modifying any component effect still adds interaction parameters that
  the data must then constrain toward zero.

## Identifiability

A relative effect is uniquely estimable only if its component contrast
lies in the row space of the within-study component design (Wigle et al.
2026);
[`relative_effects()`](https://choxos.github.io/cpaic/reference/relative_effects.md)
returns `NA` otherwise rather than a prior-driven number. Note this
checks identification of `beta`; a component x effect-modifier
interaction is additionally identified only by covariate variation on
the contrasts that involve it, and interactions informed only by
aggregate arms are weakly identified (`prior_gamma_scale` regularizes).

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
