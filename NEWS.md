# cpaic 0.1.0

First release. Experimental research software. The methodology and
implementation are not validated for clinical, regulatory, reimbursement,
guideline, or other decision use.

## Methodology safeguards

* `cstc()` and `cmaic()` stop by default when adjusted IPD edges would be
  mixed with retained aggregate edges from other study populations. Non-Gaussian
  cMAIC bridges also stop because marginal effects are not generally
  component-additive on nonlinear scales. `allow_experimental_bridge = TRUE`
  records a research-only override and its reasons.
* Partial-IPD multi-arm replacements are rejected. An IPD-only study edge
  requires `allow_ipd_only_studies = TRUE` and is recorded in the result.
* cML-NMR reporting and ranking methods identify their supported estimand as
  `average_conditional_link` at target effect-modifier means. Marginal
  standardized effects are not available. Bernoulli means in `[0, 1]` are
  valid prevalences.
* Two-stage network construction rejects missing analysis values, malformed
  study labels, invalid family-specific outcomes, and ambiguous survival event
  coding. Two-stage survival requires `0 = censored`, `1 = event`.
* Ranking summaries and rank probabilities share one retained set and one exact
  tie rule. Tied elements split probability across occupied ranks. Target-grid
  failures remain explicit rows and plots break at them.
* Poisson aggregate rates are calculated with log-sum-exp and fitted with the
  log-rate likelihood. Predictive RNG overflow is represented by explicit
  sentinels and counts. Predictive summaries require valid draws and never use
  capped rates. Posterior fitted means average response-scale rates.
* Sampler diagnostic extraction returns unknown values as `NA`, not zero. Fits
  store a `passed`, `failed`, or `unknown` diagnostic status.
* Relative-effect tables separate `estimate_link` and `se_link` from the
  reporting-scale `estimate`, interval, and `scale` field.
* Package and pkgdown builds include the current introduction and methods
  articles. Five archived rendered examples remain outside published builds
  until they are regenerated against the current interfaces.

## What the package does

cpaic combines component network meta-analysis (cNMA) with experimental
population-adjustment methods for a **disconnected** treatment network. Shared
components can identify a bridge. Scientific transportability across that gap
remains an untestable assumption.

* `cpaic_network()` builds a (possibly disconnected) contrast-level network and
  codes multi-component treatment labels into a treatment-by-component matrix
  with `build_C_matrix()`.
* `cnma_bridge()` reconnects the network through the additive component model of
  Rücker et al. (2020), on top of `netmeta::discomb()`.
* `cstc()` and `cmaic()` are the two-stage frequentist routes. They replace
  IPD-bearing edges with target-specific contrasts. A gate prevents automatic
  synthesis with incompatible retained edges or nonlinear marginal cMAIC
  contrasts.
* `cmlnmr()` is the Bayesian flagship: component-additive multilevel network
  meta-regression, fitted with rstan or CmdStan. The treatment effect is `C %*% beta`
  and the model carries component by effect-modifier interactions through the
  whole network, so disconnected sub-networks are connected by construction.
  Reported target summaries are average conditional link-scale effects at
  covariate means, not marginally standardized effects.

Binary, continuous, count, and time-to-event outcomes are supported throughout.

## Estimability is checked, not assumed

Reconnecting a network does not make the effects you want estimable. A relative
effect is uniquely estimable exactly when its contrast lies in the row space of
the component design `X = B C` (Wigle et al. 2026). Both engines would otherwise
return a confident-looking number for a contrast carrying no information: the
frequentist fit through the Moore-Penrose pseudoinverse, the Bayesian fit through
the prior.

* `cpaic_connectivity()` and `estimable_effects()` report the rank, the null
  space, the bridging components, and which relative effects are identified.
* `estimable_effects_at()` extends the criterion to the target-mean average
  conditional estimand and grades each contrast
  `"exact"`, `"first-order screen"`, or `"not identified"`.
* `relative_effects()`, `component_effects()`, and `league_table()` return `NA`
  for a contrast that is not identified; `cpaic_ranks()` drops it from the
  hierarchy rather than ranking a prior.

## Failing closed

* `cstc()` and `cmaic()` reject a study whose regression or weight fit did not
  converge, is separated, is rank deficient, has a degenerate
  treatment-coefficient covariance, or (for `cmaic()`) did not achieve moment
  balance. Each `cmaic()` bootstrap replicate is held to the same weight-validity
  gate as the point estimate, and an edge whose successful replicates fall below
  `min_boot_success` is refused rather than given a fragile standard error. An
  invalid edge is never passed silently into the additive bridge.
* `cmlnmr()` validates its inputs before compiling the Stan model: a study
  present in both `ipd` and `agd`, a single-arm study, a non-numeric or
  incomplete effect modifier, fractional aggregate counts, a malformed seed, a
  non-positive `chains`/`iter_warmup`/`iter_sampling`, and protected sampler
  arguments such as `data` in `...` are all rejected by name. Zero warmup is
  refused explicitly: rstan accepts it and samples without ever adapting the
  step size, which returns draws that look like a fit and are not one.
* Sampler arguments in `...` are checked against what the chosen backend
  actually accepts, so a cmdstanr-only argument (or a misspelled one) is named
  rather than handed to rstan to kill every chain with a message about `@`.
* `cpaic_network()` rejects self-comparisons, duplicate {study, treatment-pair}
  rows, missing identifiers, missing or nonfinite analysis values, invalid
  family-specific outcomes, and malformed survival status; `build_C_matrix()`
  rejects empty component tokens and an `inactive` label matching no component.
* An unset seed is drawn and recorded, so an unseeded fit reproduces.

## Diagnostics

* `additivity_test()` reports the Cochran Q statistics, and says explicitly that
  a saturated model (zero residual degrees of freedom) gives `Q = 0` by
  arithmetic rather than as evidence of fit, and that neither statistic can test
  whether component effects are constant *across* sub-networks.
* `edge_influence()` reports the weight each edge carries on a requested
  contrast, using `1 / (seTE^2 + tau^2)`, and warns when individual patient data
  sit on an edge that cannot affect the answer. The effective sample size from
  `effective_sample_size()` cannot detect that; `weight_diagnostics()` exposes
  weight concentration that the effective sample size also hides.
* `bridge_fragility()` quantifies how much un-testable cross-sub-network drift
  would overturn a conclusion.
* `cmlnmr()` checks divergences, tree depth, E-BFMI, R-hat, and effective sample
  size across every sampled parameter block. It reports `NA` when a diagnostic
  is unavailable and records a formal sampler-validity status.
* `posterior_summary()` summarizes any sampled parameter block of a `cmlnmr()`
  fit, with the same columns whichever sampler backend produced it.
* `prior_sensitivity()`, `prior_predictive_check()`, `dic()`, `loo()`, `waic()`,
  and `redact_fit()` cover prior movement, prior implications, model comparison,
  and sharing a fit without row-level data.
* Plots: network, forest, rankogram, deviance, leverage, prior-versus-posterior,
  integration error, MCMC, and survival curves, plus three specific to cpaic (the
  target-mean rank curve, the estimability map, and edge influence).

## Documented limitations

These are stated in the manual pages, not only here.

* **The bridging assumption is untestable.** There is by construction no
  cross-gap evidence against which to test that component effects are constant
  across sub-networks.
* **Only IPD edges are adjusted by the two-stage routes.** `cstc()` and `cmaic()`
  leave every aggregate-only edge in its own study population. The default gate
  refuses that mixed bridge.
* **Marginal effects do not add.** `cmaic()` targets a marginal effect, and on a
  non-collapsible scale the additive component model is false; the resulting bias
  survives perfect matching and infinite data.
* **`cmlnmr()` reports an average conditional link-scale contrast at target
  means**, not a population-standardized marginal effect.
* **One `Gamma` serves both roles.** It multiplies individual covariates and
  aggregate study means alike, so an interaction supported only by aggregate arms
  is an ecological association read as effect modification.
* **Survival status contracts differ.** `cmlnmr()` accepts 0 right-censored, 1
  event, 2 left-censored, and 3 interval-censored. The two-stage constructor
  accepts only 0 right-censored and 1 event.
* Further approximations recorded in `?cmlnmr`: the Gaussian model has one
  residual standard deviation for the whole network; Poisson aggregate arms
  assume person-time is independent of the effect modifiers; the
  observed-to-latent copula correlation is approximate for a non-normal margin;
  and the aggregate likelihood carries a finite quasi-Monte-Carlo integration
  error.

## Sampler backends

`cmlnmr()` fits through either engine, selected with `backend`.

* `"rstan"` is the default. The Stan models are compiled when cpaic is
  installed, so the Bayesian engine works with no further setup, and the
  examples and tests run anywhere.
* `"cmdstanr"` fits the identical models with CmdStan, which tracks Stan
  releases more closely and is often faster. It needs the `cmdstanr` package and
  a separate CmdStan installation.

The two are intended to fit the same Stan programs but do not share a random
number stream,
so the same `seed` gives different draws on each. Backend parity checks
are not a validation of every family and parameterization. Diagnostics are
computed the same way for both, so `rhat`, `ess_bulk`, and `ess_tail` mean the
same thing whichever produced a fit, and every summary, plot, and diagnostic in
the package works on either. The engine and its version are recorded in the
fit's provenance, and in the arguments `prior_sensitivity()` refits with, so a
sensitivity analysis cannot silently switch engines.

Reach for `posterior_summary()` rather than the sampler object under `fit$fit`.
That object is an S4 `stanfit` under one backend and an R6 object under the
other, and they share no accessors, so code written against one fails on the
other.

`adapt_delta` and `max_treedepth` are arguments of `cmlnmr()` rather than being
passed through `...`, because the two engines take them in different places.
Anything still passed through `...` reaches the chosen sampler unchanged and is
therefore backend-specific.

## Dependencies

Imports are kept to what is load-bearing: `netmeta` for the component-NMA engine,
`maicplus` for the MAIC weights, `randtoolbox` for the Sobol' integration points,
`igraph` for network connectivity, `loo` for the `loo`/`waic` generics,
`posterior` for the convergence quantities (which `loo` already brings in), and
otherwise packages that ship with R. The component-additive ML-NMR models and
their quasi-Monte-Carlo integration are implemented in cpaic itself; `multinma`
is a `Suggests` used only by the test that keeps the random-effects correlation
in step with multinma's `RE_cor()`. `cmlnmr()` fits through `rstan` by default;
`backend = "cmdstanr"` additionally needs `cmdstanr`, from
<https://stan-dev.r-universe.dev>.
