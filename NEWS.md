# cpaic 0.1.0

First release. Research software: read the limitations before using a result for
a decision.

## What the package does

cpaic extends component network meta-analysis (cNMA) to population-adjusted
indirect comparison (PAIC), so a **disconnected** treatment network can be
reconnected through shared treatment components and, at the same time, adjusted
for between-study imbalance in effect modifiers.

* `cpaic_network()` builds a (possibly disconnected) contrast-level network and
  codes multi-component treatment labels into a treatment-by-component matrix
  with `build_C_matrix()`.
* `cnma_bridge()` reconnects the network through the additive component model of
  Rücker et al. (2020), on top of `netmeta::discomb()`.
* `cstc()` and `cmaic()` are the two-stage frequentist routes: anchored simulated
  treatment comparison and anchored matching-adjusted indirect comparison replace
  each IPD-bearing edge with a population-adjusted contrast, which the component
  bridge then combines.
* `cmlnmr()` is the Bayesian flagship: component-additive multilevel network
  meta-regression, fitted with rstan or CmdStan. The treatment effect is `C %*% beta`
  and the model carries component by effect-modifier interactions through the
  whole network, so disconnected sub-networks are connected by construction and
  every edge is adjusted to one target population coherently.

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
* `estimable_effects_at()` extends the criterion to the population-adjusted
  estimand, which depends on the target population, and grades each contrast
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
  incomplete effect modifier, fractional aggregate counts, a malformed seed, and
  protected sampler arguments such as `data` in `...` are all rejected by name.
* `cpaic_network()` rejects self-comparisons, duplicate {study, treatment-pair}
  rows, and missing treatment labels; `build_C_matrix()` rejects empty component
  tokens and an `inactive` label matching no component.
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
  size across every sampled parameter block, and reports `NA` rather than an
  ideal infinity when a diagnostic is unavailable.
* `posterior_summary()` summarizes any sampled parameter block of a `cmlnmr()`
  fit, with the same columns whichever sampler backend produced it.
* `prior_sensitivity()`, `prior_predictive_check()`, `dic()`, `loo()`, `waic()`,
  and `redact_fit()` cover prior movement, prior implications, model comparison,
  and sharing a fit without row-level data.
* Plots: network, forest, rankogram, deviance, leverage, prior-versus-posterior,
  integration error, MCMC, and survival curves, plus three specific to cpaic (the
  population-dependent rank curve, the estimability map, and edge influence).

## Documented limitations

These are stated in the manual pages, not only here.

* **The bridging assumption is untestable.** There is by construction no
  cross-gap evidence against which to test that component effects are constant
  across sub-networks.
* **Only IPD edges are adjusted by the two-stage routes.** `cstc()` and `cmaic()`
  leave every aggregate-only edge in its own study population. `cmlnmr()` is the
  coherent single-target synthesis.
* **Marginal effects do not add.** `cmaic()` targets a marginal effect, and on a
  non-collapsible scale the additive component model is false; the resulting bias
  survives perfect matching and infinite data.
* **`cmlnmr()` reports a conditional contrast at a covariate profile**, not a
  population-standardized marginal effect.
* **One `Gamma` serves both roles.** It multiplies individual covariates and
  aggregate study means alike, so an interaction supported only by aggregate arms
  is an ecological association read as effect modification.
* **`cmlnmr()` survival status coding** (0 right, 1 event, 2 left, 3 interval) is
  not the coding `cstc()` and `cmaic()` pass to `survival::Surv()`.
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

The two agree up to Monte Carlo error but do not share a random number stream,
so the same `seed` gives different draws on each. Convergence diagnostics are
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
