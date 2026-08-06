# Marginal effects from a component ML-NMR fit

Standardizes treatment-specific posterior outcomes over an explicit
target distribution, then forms treatment contrasts within each
posterior draw. This is posterior standardization, not evaluation at a
single covariate profile. Treatment-specific outcomes are averaged
first, and contrasts are calculated second.

## Usage

``` r
marginal_effects(
  object,
  target,
  weights = NULL,
  reference = NULL,
  all_contrasts = FALSE,
  measure = NULL,
  baseline_study = NULL,
  times = NULL,
  random_effect = "population",
  backtransf = TRUE,
  level = 0.95,
  ...
)
```

## Arguments

- object:

  A fitted
  [`cmlnmr()`](https://choxos.github.io/cpaic/reference/cmlnmr.md)
  object.

- target:

  Either a data frame containing exactly the fitted effect modifiers, or
  a list with named `means`, `sds`, `margins`, and optional `cor` and
  `n_int` entries.

- weights:

  Optional nonnegative target-row weights. This is only valid when
  `target` is a data frame.

- reference:

  Reference treatment.

- all_contrasts:

  Return every ordered treatment contrast if `TRUE`.

- measure:

  Marginal contrast measure. Available measures depend on the outcome
  family.

- baseline_study:

  Study whose fitted intercept, and for survival whose fitted baseline
  hazard, is transported to the target population.

- times:

  Positive prediction times for survival measures, restricted to the
  observed follow-up support of `baseline_study`. Predictions begin at
  model time zero rather than a landmark or delayed-entry time. For RMST
  these are the integration horizons. RMST measures require a
  piecewise-exponential donor baseline and are analytic over its
  intervals.

- random_effect:

  Random-effect prediction policy. Only the population mean, which sets
  study-arm deviations to zero, is currently available.

- backtransf:

  Report ratio measures on their natural scale if `TRUE`.

- level:

  Credible interval level.

- ...:

  Unused.

## Value

A `cpaic_effects` data frame. `estimate`, `lower`, and `upper` are on
the reporting scale. `estimate_contrast`, `se_contrast`, and
`contrast_scale` describe the draw-level contrast used for inference: a
log contrast for ratios and a natural-scale difference for difference
measures. Attributes record the measure, target nodes and weights, donor
baseline, random-effect policy, and sampler diagnostic status.

## Measures

Available `measure` values are:

- binomial: `"odds_ratio"`, `"risk_ratio"`, and `"risk_difference"`;

- Gaussian: `"mean_difference"`;

- Poisson: `"rate_ratio"` and `"rate_difference"`, for unit exposure;

- survival: `"survival_difference"`, `"survival_ratio"`,
  `"risk_difference"`, `"risk_ratio"`, `"rmst_difference"`,
  `"rmst_ratio"`, and `"time_specific_hazard_ratio"`.

A scalar marginal hazard ratio is not defined. The time-specific
marginal hazard ratio generally changes with time even when the fitted
conditional hazards are proportional.

Survival and risk ratios, restricted mean survival time ratios, and
time-specific marginal hazard ratios are accumulated and contrasted on
the log scale with log-sum-exp calculations before optional
back-transformation. Set `backtransf = FALSE` to retain the log-ratio
reporting scale. RMST measures require a piecewise-exponential donor
baseline and are evaluated analytically over its intervals. M-spline
fits support marginal survival, risk, and time-specific marginal hazard
measures, but not RMST.

## Target distribution

A data-frame target supplies empirical integration rows. `weights` are
normalized, and zero-weight rows are removed. Bernoulli modifiers must
be represented by actual 0/1 rows; a fractional prevalence is not an
empirical pseudo-patient. Continuous modifiers are accepted at their
observed values, subject to the support of the fitted margin.

A list target defines a Gaussian-copula distribution through named
`means`, `sds`, and `margins`. Optional `cor` is a latent-scale
correlation matrix and `n_int` is the deterministic Sobol' integration
size. If `cor` is omitted, the correlation retained by
[`cmlnmr()`](https://choxos.github.io/cpaic/reference/cmlnmr.md) is
used. Each continuous modifier may use a `normal`, `gamma`, `lognormal`,
or `beta` margin; fitted Bernoulli modifiers must remain Bernoulli.
Independent Bernoulli modifiers are enumerated exactly, and continuous
modifiers use deterministic Sobol' nodes. Summary nodes must reproduce
requested Bernoulli strata and non-normal means and SDs within the
implementation's moment-fidelity tolerances: Bernoulli prevalence is
checked to `max(1e-8, 0.1 * min(p, 1 - p))`, non-normal means to 0.1
target SD, and non-normal SDs to 15 percent relative error. A failed
gate stops with a request for a larger `n_int` or an empirical target.
Mixed or non-normal multivariable targets also require resolved pairwise
correlations to agree within 0.05 of a high-resolution deterministic
reference grid. The resolved nodes and weights are retained on the
result. Target-distribution uncertainty is not propagated.

## Baseline transport

Binomial measures, Poisson rate differences, and every survival measure
need `baseline_study`. Its fitted intercept is transported to the
target; survival measures also transport its fitted baseline hazard. The
target covariate distribution is always supplied separately. This is an
explicit transport assumption, not something identified by relative
treatment effects. Gaussian mean differences and Poisson rate ratios do
not need a baseline study because the common intercept cancels.

Survival predictions start at model time zero. They are not landmark
predictions and are not conditioned on delayed entry or a supplied
row-level entry time. Prediction times must be positive and within the
observed follow-up support retained for the selected donor study, which
can be shorter than the global survival support.

## Identification and random effects

Nonlinear marginal contrasts are screened at every positive-weight
target node. A passing result is labelled `"first-order screen"`, never
`"exact"`, because aggregate-data interactions and survival baselines
can still be weakly identified. Gaussian identity-link mean differences
can be labelled `"exact"` when their target-mean contrast is supported
by IPD. Nonlinear results that fail this conservative treatment-surface
check are returned as `NA` with `basis = "first-order screen failed"`.
This screen covers the treatment `beta`/`Gamma` surface, not full
identification of prognostic effects, donor intercepts, or donor
baseline-hazard parameters.

The target checks are support and moment checks, not a formal overlap or
positivity diagnostic. No overlap statistic is estimated, and
extrapolation risk remains a substantive limitation for the analyst to
assess.

Only `random_effect = "population"` is available. It sets study-arm
deviations to zero. The function does not make predictive draws over a
new study's heterogeneity distribution and does not return marginal
component effects or marginal rankings, since nonlinear standardization
does not preserve component additivity. Nonlinear marginal effects
remain treatment-level contrasts.
