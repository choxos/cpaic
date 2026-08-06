# cpaic

**Component-Based Population-Adjusted Indirect Comparison**

`cpaic` is experimental research software for combining component
network meta-analysis (cNMA) with population-adjusted indirect
comparison (PAIC) in a **disconnected** treatment network. Shared
treatment components can identify a bridge across the disconnection.
IPD-bearing comparisons can then be adjusted for differences in measured
effect modifiers, subject to the restrictions below.

> **Research use only.** The methodology and implementation have not
> been validated for clinical, regulatory, reimbursement, or other
> decision use. A successful fit and clean sampler diagnostics do not
> validate the cross-sub-network additivity or transportability
> assumptions.

Standard network meta-analysis needs a connected network. When the
network is disconnected (no common comparator links two sub-networks),
it cannot be analyzed directly. Two ideas each solve half of the
problem:

- **Component NMA** (Rücker et al. 2020) decomposes multi-component
  treatments into additive component effects. If sub-networks *share
  components*, those component effects bridge the gap; but the method
  uses aggregate data only and ignores effect-modifier imbalance.
- **PAIC** (STC, MAIC, ML-NMR) adjusts for effect-modifier imbalance
  using individual patient data (IPD), but assumes the network is
  already connected.

`cpaic` implements experimental ways to combine these ideas. The
one-stage
[`cmlnmr()`](https://choxos.github.io/cpaic/reference/cmlnmr.md) model
and the two-stage
[`cstc()`](https://choxos.github.io/cpaic/reference/cstc.md) /
[`cmaic()`](https://choxos.github.io/cpaic/reference/cmaic.md) bridge
have different estimands and different failure modes. They should not be
described as interchangeable whole-network population adjustments.

## Two senses of “anchored”

These are different words that happen to be spelled the same, and
conflating them causes real confusion:

- **Anchored PAIC** (the NICE sense): the adjusted comparison runs
  through a *common comparator arm* within each trial, rather than
  comparing single arms.
  [`cmaic()`](https://choxos.github.io/cpaic/reference/cmaic.md) and
  [`cstc()`](https://choxos.github.io/cpaic/reference/cstc.md) are
  anchored in this sense.
- **Anchored cNMA** (the Rücker/Wigle sense): an *inactive component* is
  fixed at zero in the component parameterization. Setting
  `inactive = NULL` gives the **unanchored** parameterization of Wigle &
  Béliveau (2022), in which every unit receives its own parameter and no
  anchor can be misspecified.

Bridging a disconnected network through shared components is **not**
anchoring by a common comparator. It is identification through
additivity.

## Status

- **Frequentist core.**
  [`cnma_bridge()`](https://choxos.github.io/cpaic/reference/cnma_bridge.md)
  (reconnect via components),
  [`cmaic()`](https://choxos.github.io/cpaic/reference/cmaic.md)
  (component MAIC),
  [`cstc()`](https://choxos.github.io/cpaic/reference/cstc.md)
  (component anchored STC). Binary, continuous, count, and survival
  outcomes.
- **Bayesian
  [`cmlnmr()`](https://choxos.github.io/cpaic/reference/cmlnmr.md)**
  (component-additive ML-NMR), fitted with rstan by default or with
  cmdstanr on request. Binary, continuous, count, and survival outcomes,
  with separate conditional and target-standardized marginal reporting
  paths.

The four outcome-family implementations are under active development.
They are not decision validated. Read the estimand and bridge
limitations before interpreting a result.

## Estimability

Reconnecting a network through shared components does **not** guarantee
that the effects you want are estimable. A relative effect is uniquely
estimable exactly when its contrast vector lies in the **row space** of
the component design matrix `X = B C` (Wigle et al. 2026). Full column
rank of `X` is sufficient but *not necessary*, so a rank-deficient
network can still identify useful cross-sub-network contrasts.

This matters because both engines will otherwise hand back a
confident-looking number for a contrast that carries no information: the
frequentist fit through the Moore-Penrose pseudoinverse, and the
Bayesian fit through the prior. `cpaic` checks every contrast and
returns `NA` instead.

``` r

estimable_effects(net)     # which relative effects are identified?
cpaic_connectivity(net)    # rank, null space, bridging components
```

## Installation

``` r

# install.packages("remotes")
remotes::install_github("choxos/cpaic")
```

`cpaic` builds on `netmeta` for the component-NMA engine and `maicplus`
for the MAIC weights, both on CRAN. The component-additive ML-NMR models
and their quasi-Monte-Carlo integration are implemented in the package
itself, following Phillippo et al. (2020).

Those models are compiled when cpaic is installed, so
[`cmlnmr()`](https://choxos.github.io/cpaic/reference/cmlnmr.md) works
with no further setup. It fits through **rstan** by default; pass
`backend = "cmdstanr"` to use CmdStan instead, which tracks Stan
releases more closely and is often faster, and which needs the
`cmdstanr` package from <https://stan-dev.r-universe.dev> plus a CmdStan
installation. The two fit the same models and agree up to Monte Carlo
error; they do not share a random number stream, so the same `seed`
gives different draws on each. Under rstan the chains run one after
another unless you set `options(mc.cores = ...)`, as usual for R.

Summarize a fit with
[`posterior_summary()`](https://choxos.github.io/cpaic/reference/posterior_summary.md)
rather than the sampler object under `fit$fit`: that slot is an S4
`stanfit` on one backend and an R6 object on the other, and they share
no accessors.

## Quick start

``` r

library(cpaic)

# A disconnected network: sub-network {Placebo, A, B} and isolated
# sub-network {A+B, A+B+C, A+B+D}, bridged by the shared components A and B.
net <- cpaic_network(cpaic_bin_agd, ipd = cpaic_bin_ipd, sm = "OR",
                     family = "binomial", ipd_covariates = "x1",
                     inactive = "Placebo")

# Which effects can this network actually identify?
cpaic_connectivity(net)
estimable_effects(net)

# 1. Connect only (aggregate component NMA)
cnma_bridge(net)

# 2. Experimental two-stage STC bridge. This network retains AgD edges, so
#    explicit opt-in is required because the bridge mixes populations.
cstc(net, target = c(x1 = 0), effect_modifiers = "x1",
     allow_experimental_bridge = TRUE)

# 3. Experimental two-stage MAIC bridge. On the OR scale, marginal weighted
#    contrasts also need not obey component additivity.
cmaic(net, target = c(x1 = 0), effect_modifiers = "x1",
      allow_experimental_bridge = TRUE)
```

By default, [`cstc()`](https://choxos.github.io/cpaic/reference/cstc.md)
and [`cmaic()`](https://choxos.github.io/cpaic/reference/cmaic.md) stop
when retained aggregate contrasts would be combined with target-specific
adjusted contrasts.
[`cmaic()`](https://choxos.github.io/cpaic/reference/cmaic.md) also
stops by default for non-Gaussian families because a marginal weighted
contrast is not generally additive across components on a nonlinear
scale. `allow_experimental_bridge = TRUE` records the override in the
fitted object and emits a warning. It does not make the estimands
coherent.

For the Bayesian model, conditional relative effects vary with the
effect modifiers because `theta_t(x) = C_t' (beta + Gamma x)`. A one-row
`newdata` value is interpreted as a vector of target covariate means.
Because this contrast is linear in `x` on the link scale, it is an
average conditional link-scale effect at those means.

``` r

fit <- cmlnmr(ipd, agd, effect_modifiers = "x1", inactive = "Placebo")
# Average conditional link-scale effects at target mean x1 = 0.3
relative_effects(fit, newdata = data.frame(x1 = 0.3))
component_effects(fit, newdata = data.frame(x1 = 0.3))
```

Marginal reporting uses a separate target-distribution interface. It
averages treatment-specific outcomes within every posterior draw, then
forms the requested contrast. For example, an empirical binary target
distribution can be supplied as actual 0/1 rows with weights:

``` r

target_population <- data.frame(x1 = c(0, 1))
target_weights <- c(0.7, 0.3)

# Binomial marginal risk ratio. Replace "study_id" with a study retained by fit.
marginal_effects(
  fit, target = target_population, weights = target_weights,
  measure = "risk_ratio", baseline_study = "study_id"
)
```

The donor study supplies the fitted intercept needed for absolute risks.
The target rows remain a separate input. Gaussian mean differences and
Poisson rate ratios do not need a donor baseline because the common
intercept cancels.

An empirical target is a data frame with exactly the fitted effect
modifiers. Its finite, nonnegative weights are normalized and
zero-weight rows are discarded. Bernoulli modifiers must be represented
by actual 0/1 rows. A summary target is a named list of `means`, `sds`,
and `margins`, with optional latent-scale `cor` and deterministic Sobol’
`n_int`. Continuous modifiers can use `normal`, `gamma`, `lognormal`, or
`beta` margins independently; fitted Bernoulli modifiers must remain
Bernoulli. Independent Bernoulli strata are enumerated exactly, while
continuous margins use deterministic integration nodes. Summary targets
have moment-fidelity gates for Bernoulli strata and non-normal means and
SDs. The Bernoulli prevalence tolerance is
`max(1e-8, 0.1 * min(p, 1 - p))`; non-normal means must be within 0.1
target SD and non-normal SDs within 15 percent relative error. A failed
gate stops with a request for more nodes or an empirical target. Mixed
or non-normal multivariable targets also require resolved pairwise
correlations to agree within 0.05 of a high-resolution deterministic
reference grid. Target-distribution uncertainty is not propagated, and
these checks do not constitute a formal overlap or positivity
diagnostic.

For survival, predictions start at model time zero and are not landmark
or delayed-entry conditioned. `times` must lie within the observed
follow-up support of the selected donor study, which may be shorter than
the global support. Survival and risk ratios, piecewise-exponential RMST
ratios, and time-specific marginal hazard ratios are calculated on the
log scale with stable log-sum-exp operations before optional
back-transformation. RMST measures require a piecewise-exponential donor
baseline and are evaluated analytically over its intervals. M-spline
fits support marginal survival, risk, and time-specific marginal hazard
measures, but not RMST. Nonlinear marginal effects are treatment-level
contrasts; the package does not report marginal component effects or
rankings.

[`relative_effects()`](https://choxos.github.io/cpaic/reference/relative_effects.md),
[`marginal_effects()`](https://choxos.github.io/cpaic/reference/marginal_effects.md),
[`league_table()`](https://choxos.github.io/cpaic/reference/league_table.md),
[`component_effects()`](https://choxos.github.io/cpaic/reference/component_effects.md),
[`estimable_effects()`](https://choxos.github.io/cpaic/reference/estimable_effects.md),
[`additivity_test()`](https://choxos.github.io/cpaic/reference/additivity_test.md),
[`effective_sample_size()`](https://choxos.github.io/cpaic/reference/effective_sample_size.md),
[`forest()`](https://choxos.github.io/cpaic/reference/forest.md) and
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) summarize and
visualize a fit.

## Validation status

The package does not yet have a tracked, independently reproducible
validation suite that supports decision use. Local exploratory
simulations are useful for finding failure modes, but they are not
release evidence and are not summarized as performance claims here. At
minimum, future validation must cover overlap, model misspecification,
partial IPD coverage, nonlinear estimands, multi-arm studies, survival
reconstruction, and violations of cross-sub-network additivity.

## Limitations

The following apply to any result produced by this package.

- **The bridging assumption is untestable.** Reconnecting a disconnected
  network requires component effects (and, under population adjustment,
  component by effect-modifier interactions) to be *constant across
  sub-networks*. There is by construction no cross-gap evidence to test
  this against.
  [`additivity_test()`](https://choxos.github.io/cpaic/reference/additivity_test.md)
  reports the fit of the additive model *within* the observed evidence;
  a large p-value is not a license to bridge.
- **A two-stage bridge can mix populations.**
  [`cmaic()`](https://choxos.github.io/cpaic/reference/cmaic.md) and
  [`cstc()`](https://choxos.github.io/cpaic/reference/cstc.md) adjust
  the edges that have IPD. Retained aggregate contrasts remain on their
  original study-population estimands. Combining both types does not
  produce a coherent single-target network merely because the component
  model can solve the algebra. The default is to stop;
  `allow_experimental_bridge = TRUE` is an explicit research override,
  not a validity correction.
- **Estimands differ across the methods.**
  [`cstc()`](https://choxos.github.io/cpaic/reference/cstc.md) returns a
  conditional effect on the link scale at the target covariate means;
  [`cmaic()`](https://choxos.github.io/cpaic/reference/cmaic.md) returns
  a marginal weighted contrast. On a nonlinear scale a marginal contrast
  is not generally component-additive, so feeding cMAIC edges into an
  additive bridge can be methodologically incoherent even with perfect
  matching. Non-Gaussian cMAIC bridges therefore require the
  experimental override.
- **Conditional and marginal cML-NMR summaries are separate estimands.**
  A one-row `newdata` profile reports an average conditional link-scale
  contrast.
  [`marginal_effects()`](https://choxos.github.io/cpaic/reference/marginal_effects.md)
  requires a target distribution and standardizes absolute outcomes
  before forming treatment-level contrasts. Nonlinear marginal effects
  are not component-additive, so marginal component effects and marginal
  treatment rankings are not reported.
- **Some marginal measures borrow a study baseline.** Binomial effects,
  Poisson rate differences, and every survival measure transport a
  selected fitted study intercept or baseline hazard to the target
  population. This is an explicit assumption. Target-distribution
  uncertainty is not propagated, and study-arm random effects are set to
  zero for population-mean predictions. Survival predictions begin at
  model time zero, not at a landmark or delayed entry time, and
  prediction times are restricted to the selected donor study’s observed
  follow-up support.
- **Target checks are not overlap diagnostics.** Empirical and summary
  targets are checked for finite values, fitted-margin support,
  Bernoulli strata, and summary moment fidelity. The package does not
  estimate a formal overlap or positivity statistic, so extrapolation
  risk remains an analyst-level responsibility.
- **Survival conditions on the supplied rows.**
  [`cmlnmr()`](https://choxos.github.io/cpaic/reference/cmlnmr.md)
  evaluates the exact analytic individual survival likelihood (events
  plus right, left, and interval censoring, and delayed entry),
  integrated over each aggregate arm’s covariate distribution, with a
  study-specific piecewise-exponential or continuous cubic M-spline
  baseline. It is exact only *conditional on* the supplied (or
  reconstructed) pseudo-individual rows, the chosen baseline basis, and
  the finite numerical integration; uncertainty from reconstructing
  pseudo-IPD out of a published Kaplan-Meier curve is not propagated,
  and proportional hazards is assumed.
  [`marginal_effects()`](https://choxos.github.io/cpaic/reference/marginal_effects.md)
  can report standardized survival and risk contrasts and time-specific
  marginal hazard ratios within the fitted follow-up support. Restricted
  mean survival time contrasts are available for piecewise-exponential
  baselines. M-spline fits retain the pointwise survival, risk, and
  time-specific hazard measures. These outputs additionally require
  transport of a selected study’s baseline hazard.
- **Survival interfaces have different status contracts.** The two-stage
  [`cstc()`](https://choxos.github.io/cpaic/reference/cstc.md) and
  [`cmaic()`](https://choxos.github.io/cpaic/reference/cmaic.md) network
  input accepts only `0` for right censoring and `1` for an event.
  [`cmlnmr()`](https://choxos.github.io/cpaic/reference/cmlnmr.md)
  separately uses `0` right-censored, `1` event, `2` left-censored, and
  `3` interval-censored rows, with reconstructed rows required for
  aggregate survival arms. The package reports proportional-hazard
  contrasts through
  [`relative_effects()`](https://choxos.github.io/cpaic/reference/relative_effects.md).
  A scalar marginal hazard ratio is not reported because the
  standardized hazard ratio generally varies with time.
- **[`cmlnmr()`](https://choxos.github.io/cpaic/reference/cmlnmr.md)
  treatment effects may be fixed or random.** `trt_effects = "random"`
  adds study-arm heterogeneity with a single common `tau`, which is one
  exchangeability assumption, not a component-, sub-network-, or
  era-specific heterogeneity structure.

## Documentation

The introductory and methods vignettes ship with the package. Five
archived outcome examples are excluded until they are regenerated
against the current estimand labels, safeguards, and API. The standalone
manual in the development-only `documentation/` tree is not public
package guidance.

## References

- Rücker G, Petropoulou M, Schwarzer G (2020). Network meta-analysis of
  multicomponent interventions. *Biometrical Journal* 62(3):808-821.
- Wigle A, Béliveau A (2022). Bayesian unanchored additive models for
  component network meta-analysis. *Statistics in Medicine*
  41(22):4444-4466.
- Wigle A, Béliveau A, Nikolakopoulou A, Lin L (2026). Creating
  treatment and component hierarchies in component network
  meta-analysis.
- Efthimiou O, et al. (2022). A Bayesian model for combining aggregate
  and individual participant data in component network meta-analysis.
  *Statistics in Medicine* 41(14):2586-2606.
- Phillippo DM, et al. (2020). Multilevel network meta-regression for
  population-adjusted treatment comparisons. *JRSS A* 183(3):1189-1210.

## License

GPL-3
