# cpaic

**Component-Based Population-Adjusted Indirect Comparison**

`cpaic` extends component network meta-analysis (cNMA) to
population-adjusted indirect comparison (PAIC) methods, so that a
**disconnected** treatment network can be *reconnected* through shared
treatment components and, at the same time, *adjusted* for between-study
differences in effect modifiers.

Standard network meta-analysis needs a connected network. When the
network is disconnected (no common comparator linking two sub-networks),
it cannot be analyzed directly. Two ideas each solve part of the
problem:

- **Component NMA** (Rücker et al. 2020) decomposes multi-component
  treatments into additive component effects; if sub-networks *share
  components*, the component effects bridge the gap — but on aggregate
  data only, ignoring effect-modifier imbalance.
- **PAIC** (STC, MAIC, ML-NMR) adjusts for effect-modifier imbalance
  using individual patient data (IPD), but assumes a connected/anchored
  network.

`cpaic` combines them: the component structure provides the bridge, and
anchored STC / MAIC / ML-NMR provide the population adjustment.

## Status

- **Frequentist core — available.**
  [`cnma_bridge()`](https://choxos.github.io/cpaic/reference/cnma_bridge.md)
  (reconnect via components),
  [`cmaic()`](https://choxos.github.io/cpaic/reference/cmaic.md)
  (component MAIC),
  [`cstc()`](https://choxos.github.io/cpaic/reference/cstc.md)
  (component anchored STC), with binary, continuous, count, and survival
  outcomes.
- **Bayesian flagship
  [`cmlnmr()`](https://choxos.github.io/cpaic/reference/cmlnmr.md)
  (component-additive ML-NMR) — available.** Binary, continuous, count,
  and survival outcomes; fitted with `cmdstanr`. Survival uses a
  piecewise-exponential proportional-hazards model with a flexible
  (step-function) baseline set by `cut_points`.

## Installation

``` r

# install.packages("remotes")
remotes::install_github("choxos/cpaic")
```

`cpaic` builds on `netmeta` (cNMA), `maicplus` (MAIC weights), and
`multinma` (ML-NMR / numerical integration), all on CRAN.

## Quick start

``` r

library(cpaic)

# A disconnected network: sub-network {Placebo, A, B} and isolated
# sub-network {A+B, A+B+C, A+B+D}, bridged by the shared components A and B.
net <- cpaic_network(cpaic_bin_agd, ipd = cpaic_bin_ipd, sm = "OR",
                     family = "binomial", ipd_covariates = "x1",
                     inactive = "Placebo")

# Is it bridgeable? (disconnected, but components make it identifiable)
cpaic_connectivity(net)

# 1. Connect only (aggregate component NMA)
cnma_bridge(net)

# 2. Connect + population-adjust via anchored STC
cstc(net, target = c(x1 = 0), effect_modifiers = "x1")

# 3. Connect + population-adjust via anchored MAIC
cmaic(net, target = c(x1 = 0), effect_modifiers = "x1")
```

[`relative_effects()`](https://choxos.github.io/cpaic/reference/relative_effects.md),
[`league_table()`](https://choxos.github.io/cpaic/reference/league_table.md),
[`component_effects()`](https://choxos.github.io/cpaic/reference/component_effects.md),
[`additivity_test()`](https://choxos.github.io/cpaic/reference/additivity_test.md),
[`effective_sample_size()`](https://choxos.github.io/cpaic/reference/effective_sample_size.md),
[`forest()`](https://choxos.github.io/cpaic/reference/forest.md) and
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) summarize and
visualize a fit.

## Documentation

Technical documentation (mathematical foundations, methods manual,
validation study) lives in the local `documentation/` folder.
User-facing vignettes are installed with the package.

## References

- Rücker G, Petropoulou M, Schwarzer G (2020). Network meta-analysis of
  multicomponent interventions. *Biometrical Journal* 62(3):808–821.
- Phillippo DM, et al. (2020). Multilevel network meta-regression for
  population-adjusted treatment comparisons. *JRSS A* 183(3):1189–1210.

## License

GPL-3
