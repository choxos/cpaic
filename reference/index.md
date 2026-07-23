# Package index

## Network setup

Build a (possibly disconnected) component network and code components.

- [`cpaic_network()`](https://choxos.github.io/cpaic/reference/cpaic_network.md)
  : Set up a (possibly disconnected) component network for cpaic
- [`build_C_matrix()`](https://choxos.github.io/cpaic/reference/build_C_matrix.md)
  : Build a component-coded treatment-by-component matrix

## Connectivity and estimability

Detect disconnection, and check which relative effects the component
design can actually identify. Reconnecting a network does not guarantee
that the effects you want are estimable, and both engines will otherwise
return a confident-looking number for a contrast that carries no
information.

- [`cpaic_connectivity()`](https://choxos.github.io/cpaic/reference/cpaic_connectivity.md)
  : Assess connectivity and component-bridge identifiability of a
  network
- [`estimable_effects()`](https://choxos.github.io/cpaic/reference/estimable_effects.md)
  : Which relative effects of a component network are uniquely
  estimable?
- [`estimable_effects_at()`](https://choxos.github.io/cpaic/reference/estimable_effects_at.md)
  : Which population-adjusted contrasts are estimable at a target
  population?

## Connection layer

Reconnect a disconnected network through shared components.

- [`cnma_bridge()`](https://choxos.github.io/cpaic/reference/cnma_bridge.md)
  : Reconnect a network through its additive component structure
- [`additivity_test()`](https://choxos.github.io/cpaic/reference/additivity_test.md)
  : Fit statistics for the additive component model

## Population adjustment

Anchored, population-adjusted indirect comparison across the network.

- [`cmaic()`](https://choxos.github.io/cpaic/reference/cmaic.md) :
  Component matching-adjusted indirect comparison (cMAIC)
- [`cstc()`](https://choxos.github.io/cpaic/reference/cstc.md) :
  Component simulated treatment comparison (cSTC)
- [`cmlnmr()`](https://choxos.github.io/cpaic/reference/cmlnmr.md) :
  Component-additive multilevel network meta-regression (ML-NMR)
- [`effective_sample_size()`](https://choxos.github.io/cpaic/reference/effective_sample_size.md)
  : Effective sample sizes from a cMAIC fit
- [`weight_diagnostics()`](https://choxos.github.io/cpaic/reference/weight_diagnostics.md)
  : Weight-quality diagnostics for a cMAIC fit
- [`edge_influence()`](https://choxos.github.io/cpaic/reference/edge_influence.md)
  : Does the individual patient data actually inform this contrast?

## Hierarchies

Rank treatments or components IN A TARGET POPULATION. Because component
effects are population-specific under population adjustment, so are the
rankings: a component can lead in one population and trail in another.

- [`cpaic_ranks()`](https://choxos.github.io/cpaic/reference/cpaic_ranks.md)
  : Population-adjusted treatment and component hierarchies
- [`rank_curve()`](https://choxos.github.io/cpaic/reference/rank_curve.md)
  : How a hierarchy changes across target populations
- [`rank_probs()`](https://choxos.github.io/cpaic/reference/rank_probs.md)
  : Posterior rank probabilities in a target population

## Reporting

- [`relative_effects()`](https://choxos.github.io/cpaic/reference/relative_effects.md)
  : Relative treatment effects from a cpaic fit
- [`league_table()`](https://choxos.github.io/cpaic/reference/league_table.md)
  : League table of all pairwise relative effects
- [`component_effects()`](https://choxos.github.io/cpaic/reference/component_effects.md)
  : Component effects from a cpaic fit
- [`bridge_fragility()`](https://choxos.github.io/cpaic/reference/bridge_fragility.md)
  : Bridge fragility: how much cross-sub-network drift would change a
  conclusion

## Plots

Every plot returns a ggplot object, so it can be modified with the usual
ggplot2 verbs. The network, forest, rankogram, deviance,
prior-posterior, integration-error, MCMC, and survival plots are ported
from multinma (Phillippo et al. 2020). The rank curve, the estimability
map, and the edge-influence plot are specific to cpaic: under population
adjustment the hierarchy, and the estimable set itself, are functions of
the target population.

- [`plot(`*`<cpaic_network>`*`)`](https://choxos.github.io/cpaic/reference/plot.cpaic_network.md)
  : Plot the component network
- [`forest()`](https://choxos.github.io/cpaic/reference/forest.md)
  [`plot(`*`<cpaic_effects>`*`)`](https://choxos.github.io/cpaic/reference/forest.md)
  [`plot(`*`<cpaic_bridge>`*`)`](https://choxos.github.io/cpaic/reference/forest.md)
  [`plot(`*`<cpaic_fit>`*`)`](https://choxos.github.io/cpaic/reference/forest.md)
  : Forest plot of relative or component effects
- [`plot_rank_curve()`](https://choxos.github.io/cpaic/reference/plot_rank_curve.md)
  : How the hierarchy changes across target populations
- [`plot_estimability()`](https://choxos.github.io/cpaic/reference/plot_estimability.md)
  : Map which contrasts are estimable, and on what evidence, across
  populations
- [`plot_edge_influence()`](https://choxos.github.io/cpaic/reference/plot_edge_influence.md)
  : Plot how much each edge informs a chosen contrast
- [`plot(`*`<cpaic_rank_probs>`*`)`](https://choxos.github.io/cpaic/reference/plot.cpaic_rank_probs.md)
  : Rankogram and cumulative rank plot
- [`plot(`*`<cpaic_ranks>`*`)`](https://choxos.github.io/cpaic/reference/plot.cpaic_ranks.md)
  : Plot a population-adjusted hierarchy
- [`plot(`*`<cpaic_dic>`*`)`](https://choxos.github.io/cpaic/reference/plot.cpaic_dic.md)
  : Deviance and dev-dev plots
- [`plot_leverage()`](https://choxos.github.io/cpaic/reference/plot_leverage.md)
  : Leverage plot
- [`plot_prior_posterior()`](https://choxos.github.io/cpaic/reference/plot_prior_posterior.md)
  : Prior versus posterior
- [`plot_integration_error()`](https://choxos.github.io/cpaic/reference/plot_integration_error.md)
  : Numerical integration error against the number of integration points
- [`plot(`*`<cpaic_mlnmr>`*`)`](https://choxos.github.io/cpaic/reference/plot.cpaic_mlnmr.md)
  : MCMC diagnostics for a cML-NMR fit
- [`plot_survival()`](https://choxos.github.io/cpaic/reference/plot_survival.md)
  : Fitted survival curves from a cML-NMR fit
- [`geom_km()`](https://choxos.github.io/cpaic/reference/geom_km.md) :
  Kaplan-Meier curves from the survival data behind a cML-NMR fit

## Bayesian diagnostics

- [`dic()`](https://choxos.github.io/cpaic/reference/dic.md) : Deviance
  information criterion
- [`loo(`*`<cpaic_mlnmr>`*`)`](https://choxos.github.io/cpaic/reference/loo.cpaic_mlnmr.md)
  : Pareto-smoothed importance sampling leave-one-out cross-validation
- [`waic(`*`<cpaic_mlnmr>`*`)`](https://choxos.github.io/cpaic/reference/waic.cpaic_mlnmr.md)
  : Widely applicable information criterion
- [`prior_sensitivity()`](https://choxos.github.io/cpaic/reference/prior_sensitivity.md)
  : Refit cML-NMR under tighter and looser priors
- [`prior_predictive_check()`](https://choxos.github.io/cpaic/reference/prior_predictive_check.md)
  : Summarize a prior-predictive cML-NMR fit
- [`redact_fit()`](https://choxos.github.io/cpaic/reference/redact_fit.md)
  : Strip raw individual patient data from a fitted cML-NMR object

## Data

- [`cpaic_bin_agd`](https://choxos.github.io/cpaic/reference/cpaic_bin_agd.md)
  : Example disconnected component network: aggregate contrasts
- [`cpaic_bin_ipd`](https://choxos.github.io/cpaic/reference/cpaic_bin_ipd.md)
  : Example disconnected component network: individual patient data
