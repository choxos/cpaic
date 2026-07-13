# cpaic

**Component-Based Population-Adjusted Indirect Comparison**

`cpaic` extends component network meta-analysis (cNMA) to
population-adjusted indirect comparison (PAIC), so that a
**disconnected** treatment network can be *reconnected* through shared
treatment components and, at the same time, *adjusted* for between-study
differences in effect modifiers.

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

`cpaic` composes them: the component structure supplies the bridge, and
STC / MAIC / ML-NMR supply the population adjustment.

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
  (component-additive ML-NMR), fitted with `cmdstanr`. Binary,
  continuous, count, and survival outcomes.

This is research software under active development. Read the Limitations
section before using it for a decision.

## Estimability comes first

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

`cpaic` builds on `netmeta` (cNMA), `maicplus` (MAIC weights), and
`multinma` (ML-NMR and numerical integration), all on CRAN.
[`cmlnmr()`](https://choxos.github.io/cpaic/reference/cmlnmr.md)
additionally needs `cmdstanr`.

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

# 2. Connect and population-adjust via anchored STC
cstc(net, target = c(x1 = 0), effect_modifiers = "x1")

# 3. Connect and population-adjust via anchored MAIC
cmaic(net, target = c(x1 = 0), effect_modifiers = "x1")
```

For the Bayesian model, relative effects are **population-specific**:
the model contains component by effect-modifier interactions, so
`theta_t(x) = C_t' (beta + Gamma x)`. You must name the target
population.

``` r

fit <- cmlnmr(ipd, agd, effect_modifiers = "x1", inactive = "Placebo")
relative_effects(fit, newdata = data.frame(x1 = 0.3))   # effects at x1 = 0.3
component_effects(fit, newdata = data.frame(x1 = 0.3))
```

[`relative_effects()`](https://choxos.github.io/cpaic/reference/relative_effects.md),
[`league_table()`](https://choxos.github.io/cpaic/reference/league_table.md),
[`component_effects()`](https://choxos.github.io/cpaic/reference/component_effects.md),
[`estimable_effects()`](https://choxos.github.io/cpaic/reference/estimable_effects.md),
[`additivity_test()`](https://choxos.github.io/cpaic/reference/additivity_test.md),
[`effective_sample_size()`](https://choxos.github.io/cpaic/reference/effective_sample_size.md),
[`forest()`](https://choxos.github.io/cpaic/reference/forest.md) and
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) summarize and
visualize a fit.

## Limitations

Read these before trusting a result.

- **The bridging assumption is untestable.** Reconnecting a disconnected
  network requires component effects (and, under population adjustment,
  component by effect-modifier interactions) to be *constant across
  sub-networks*. There is by construction no cross-gap evidence to test
  this against.
  [`additivity_test()`](https://choxos.github.io/cpaic/reference/additivity_test.md)
  reports the fit of the additive model *within* the observed evidence;
  a large p-value is not a licence to bridge.
- **Only IPD edges are adjusted.**
  [`cmaic()`](https://choxos.github.io/cpaic/reference/cmaic.md) and
  [`cstc()`](https://choxos.github.io/cpaic/reference/cstc.md) adjust
  the edges that have IPD; aggregate contrasts enter the bridge as
  published, in their own study populations. The pooled result is
  therefore fully adjusted to a single target only if the retained
  aggregate effects are themselves transportable.
- **Estimands differ across the methods.**
  [`cstc()`](https://choxos.github.io/cpaic/reference/cstc.md) returns a
  conditional effect at the target covariate values;
  [`cmaic()`](https://choxos.github.io/cpaic/reference/cmaic.md) returns
  a marginal weighted effect in the target population. On
  non-collapsible scales (odds ratios, hazard ratios) these differ even
  when both are correct.
- **Survival is approximate.** The aggregate survival likelihood
  approximates expected events by person-time times mean hazard, which
  is biased upward because higher-hazard individuals leave the risk set
  earlier. The `"mspline"` baseline is a step baseline whose interval
  heights are smoothed by an M-spline; it is *not* the continuous-time
  integrated M-spline likelihood of `multinma`. Right-censoring only.
- **[`cmlnmr()`](https://choxos.github.io/cpaic/reference/cmlnmr.md) is
  a common-effect model.** No between-study heterogeneity.

## Documentation

Mathematical foundations, the methods manual, and the validation study
live in the local `documentation/` folder. User-facing vignettes ship
with the package.

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
