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
  backend = c("rstan", "cmdstanr"),
  chains = 4L,
  iter_warmup = 500L,
  iter_sampling = 500L,
  seed = NULL,
  adapt_delta = NULL,
  max_treedepth = NULL,
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
  column for Poisson outcomes. If the Poisson exposure column is absent,
  every patient is given an exposure of 1 (equal follow-up) and a
  message says so, matching
  [`cstc()`](https://choxos.github.io/cpaic/reference/cstc.md) and
  [`cmaic()`](https://choxos.github.io/cpaic/reference/cmaic.md), which
  drop the log offset when no exposure column is named.

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

- backend:

  Sampler engine: `"rstan"` (default) or `"cmdstanr"`. The two fit the
  same Stan models and are interchangeable; see the section below.

- chains, iter_warmup, iter_sampling, seed:

  Sampler settings. `iter_warmup` and `iter_sampling` are counted
  separately whichever backend is used; the translation to rstan's
  combined `iter` is handled internally.

- adapt_delta, max_treedepth:

  Sampler tuning, or `NULL` for the engine default. These are named
  arguments rather than left to `...` because the two backends take them
  in different places.

- ...:

  Further arguments for the sampler. These are passed through untouched
  and are therefore **backend-specific**: they reach
  [`rstan::sampling()`](https://mc-stan.org/rstan/reference/stanmodel-method-sampling.html)
  or the `cmdstanr` `$sample()` method as given. Prefer the named
  arguments above for anything that has one.

  Under `backend = "rstan"` an argument rstan does not accept is
  rejected by name before the fit, rather than being handed to
  [`rstan::sampling()`](https://mc-stan.org/rstan/reference/stanmodel-method-sampling.html)
  to kill every chain with a message that names nothing. About twenty
  cmdstanr sampler arguments have no rstan equivalent (`step_size`,
  `metric`, `inv_metric`, `adapt_engaged`, `parallel_chains`,
  `save_latent_dynamics`, and so on), and a misspelled argument is
  caught the same way.

## Value

An object of class `cpaic_mlnmr` with the fitted Stan object, the
component design, and a tidy table of component effects.

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
  that the estimates are stable before relying on them. (The cruder
  alternative of summarizing an aggregate arm by its event count and
  person-time was biased by 36% in a two-group example, which is why it
  is rejected outright rather than offered as a fallback.)

- **Each study has its own baseline hazard shape.** Every study carries
  its own set of spline (or step) coefficients, smoothed toward a common
  shape by a shared random-walk scale (`prior_aux_sd`), so the treatment
  effects do not have to absorb baseline misfit. A single global spline
  basis is built from the pooled follow-up range, so a study with much
  shorter follow-up may not inform the coefficients of the latest basis
  functions; those are then determined by the smoothing prior rather
  than by that study's data.

With `baseline = "mspline"` the sampler may print, during warmup,
`coefficients[...] is not a valid simplex. sum(...) = nan`. The
smoothing scale is unbounded above, so an early proposal can drive the
random walk on the log spline coefficients past the floating-point range
and the softmax returns `NaN`. Stan rejects that proposal and adaptation
continues; this is the rejection mechanism doing its job, not a
fitted-model problem. Read it as a warning only if it persists past
warmup, which would indicate a genuinely ill-conditioned baseline.

## Within-study versus ecological effect modification

A single `Gamma` multiplies the individual covariates of the IPD **and**
the covariate means of the aggregate arms. These are not the same
parameter. Write the effect as \$\$\alpha + \gamma_W (x - \bar x_s) +
\gamma_B \bar x_s .\$\$ An aggregate contrast depends only on
`alpha + gamma_B xbar_s`, so it carries no information about the
within-study interaction `gamma_W`; fitting one `Gamma` imposes
`gamma_W = gamma_B`. Randomization identifies each study's treatment
effect but does not randomize covariate means *across* studies, so a
between-study gradient is confounded in a way a within-study slope is
not (Berlin et al. 2002; Freeman et al. 2018).

The practical consequence: an interaction supported only by aggregate
arms is an **ecological** association being read as effect modification.
[`estimable_effects_at()`](https://choxos.github.io/cpaic/reference/estimable_effects_at.md)
separates the two in its `identified_by` column (`"IPD"` versus
`"aggregate"`) and marks the latter `basis = "first-order screen"`;
[`cpaic_ranks()`](https://choxos.github.io/cpaic/reference/cpaic_ranks.md)
drops such elements from a hierarchy by default. Treat a
target-population effect that leans on aggregate-identified interactions
as exploratory, and check it with
[`prior_sensitivity()`](https://choxos.github.io/cpaic/reference/prior_sensitivity.md).

## Survival status coding

`cmlnmr()` uses the four-level convention `0` right-censored, `1`
observed event, `2` left-censored, `3` interval-censored. This is
**not** the coding
[`cstc()`](https://choxos.github.io/cpaic/reference/cstc.md) and
[`cmaic()`](https://choxos.github.io/cpaic/reference/cmaic.md) use:
those pass the column straight to
[`survival::Surv()`](https://rdrr.io/pkg/survival/man/Surv.html), which
reads `0`/`1` or `1`/`2`, so a `2` there is an event rather than a
left-censored observation. Do not reuse one status column across the two
layers without recoding it.

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

- **The Gaussian model has one residual standard deviation for the whole
  network.** A single `sigma` covers every individual-level observation
  in every study and arm. Studies whose residual variance genuinely
  differs are then weighted mostly by sample size rather than by
  precision, so their relative contribution to a conflicting component
  effect, and the width of the resulting interval, are not right. Fit
  the families separately, or rescale, if the residual scales are far
  apart.

- **Poisson aggregate arms assume exposure is independent of the
  covariates.** The aggregate mean is `E * mean(exp(eta))` over the
  integration points, which equals the correct `sum_i E_i exp(eta_i)`
  only when individual person-time is unrelated to the effect modifiers.
  When longer-followed patients differ systematically the aggregate
  contribution is biased; the interface has no way to accept
  exposure-weighted covariate moments.

- **The copula correlation for a non-normal margin is approximate.** For
  a Bernoulli margin the observed-to-latent map is multinma's
  closed-form `cor_adjust = "pearson"` adjustment, which does not use
  the prevalences: at a prevalence of 0.1 a requested observed
  correlation of 0.5 comes back as about 0.42. For gamma, lognormal, and
  beta margins the observed correlation is used unadjusted and a warning
  says so. Supply `cor` on the latent scale to set it exactly.

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

## Backends

`backend = "rstan"` is the default. Its models are compiled when cpaic
is installed, so nothing else is needed and the examples and tests run
anywhere. `backend = "cmdstanr"` fits the identical models with CmdStan,
which tracks Stan releases more closely and is often faster, but it
needs the `cmdstanr` package and a separate CmdStan installation.

The two are interchangeable, not identical: they do not share a random
number stream, so the same `seed` gives different draws on each.
Convergence diagnostics are computed the same way for both (through
`posterior`), so `rhat`, `ess_bulk`, and `ess_tail` mean the same thing
whichever produced the fit, and everything downstream of the fit works
on either.

They also differ in how they run chains. cmdstanr runs all `chains` at
once. rstan follows the R convention of taking its core count from
`getOption("mc.cores")`, which is 1 unless you set it, so on the default
backend the chains run one after another until you do:

    options(mc.cores = parallel::detectCores())

Do not reach into `fit$fit` to summarize parameters. That slot holds
whatever the backend returned, an S4 `stanfit` or an R6 CmdStan object,
and the two share no accessors, so code written against one fails on the
other. Use
[`posterior_summary()`](https://choxos.github.io/cpaic/reference/posterior_summary.md),
which returns the same table either way.

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
# \donttest{
ipd <- data.frame(.study = "S1",
                  .trt = rep(c("Placebo", "A"), each = 100),
                  .y = rbinom(200, 1, 0.5), x1 = rnorm(200))
agd <- data.frame(.study = "S2", .trt = c("Placebo", "A+B"),
                  r = c(40, 55), n = c(100, 100),
                  x1_mean = c(0.2, 0.2), x1_sd = c(1, 1))
fit <- cmlnmr(ipd, agd, effect_modifiers = "x1", inactive = "Placebo",
              chains = 2, iter_warmup = 200, iter_sampling = 200)
#> Warning: Bulk Effective Samples Size (ESS) is too low, indicating posterior means and medians may be unreliable.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#bulk-ess
#> Warning: Tail Effective Samples Size (ESS) is too low, indicating posterior variances and tail quantiles may be unreliable.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#tail-ess
# Effects in a named target population (x1 = 0.2), not at the origin:
relative_effects(fit, newdata = data.frame(x1 = 0.2))
#> Relative effects (OR, back-transformed)
#>   Conditional effect at covariate profile: x1 = 0.2
#>  treatment comparator estimate    se lower upper pr_gt0              basis
#>          A    Placebo    2.121 0.286 1.141 3.557  0.995              exact
#>        A+B    Placebo    2.020 0.312 1.037 3.438  0.975 first-order screen
#>   `se` is on the link (log) scale; the interval is back-transformed.
#>   basis "first-order screen" = estimable by the row-space criterion but leaning
#>   on aggregate arms or a survival baseline, so it can be optimistic; check
#>   with prior_sensitivity() / estimable_effects_at().
# }
```
