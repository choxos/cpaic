# cpaic <img src="man/figures/logo.png" align="right" height="120" alt="" />

**Component-Based Population-Adjusted Indirect Comparison**

`cpaic` extends component network meta-analysis (cNMA) to
population-adjusted indirect comparison (PAIC) methods, so that a
**disconnected** treatment network can be *reconnected* through shared
treatment components and, at the same time, *adjusted* for between-study
differences in effect modifiers.

Standard network meta-analysis needs a connected network. When the network
is disconnected (no common comparator linking two sub-networks), it cannot
be analyzed directly. Two ideas each solve part of the problem:

- **Component NMA** (Rücker et al. 2020) decomposes multi-component
  treatments into additive component effects; if sub-networks *share
  components*, the component effects bridge the gap — but on aggregate
  data only, ignoring effect-modifier imbalance.
- **PAIC** (STC, MAIC, ML-NMR) adjusts for effect-modifier imbalance using
  individual patient data (IPD), but assumes a connected/anchored network.

`cpaic` combines them: the component structure provides the bridge, and
anchored STC / MAIC / ML-NMR provide the population adjustment.

## Status

- **Phase 1 (frequentist core) — available now.**
  `cnma_bridge()` (reconnect via components), `cmaic()` (component MAIC),
  `cstc()` (component anchored STC), with binary, continuous, count, and
  survival outcomes.
- **Phase 2 (Bayesian flagship `cmlnmr()`, component-additive ML-NMR) — in
  development.**

## Installation

```r
# install.packages("remotes")
remotes::install_github("choxos/cpaic")
```

`cpaic` builds on `netmeta` (cNMA), `maicplus` (MAIC weights), and
`multinma` (ML-NMR / numerical integration), all on CRAN.

## Quick start

```r
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

`relative_effects()`, `league_table()`, `component_effects()`,
`additivity_test()`, `effective_sample_size()`, `forest()` and
`plot()` summarize and visualize a fit.

## Documentation

Technical documentation (mathematical foundations, methods manual,
validation study) lives in the local `documentation/` folder. User-facing
vignettes are installed with the package.

## References

- Rücker G, Petropoulou M, Schwarzer G (2020). Network meta-analysis of
  multicomponent interventions. *Biometrical Journal* 62(3):808–821.
- Phillippo DM, et al. (2020). Multilevel network meta-regression for
  population-adjusted treatment comparisons. *JRSS A* 183(3):1189–1210.

## License

GPL-3
