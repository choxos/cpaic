# CLAUDE.md — cpaic developer reference

Guidance for future Claude (and human) sessions working on **cpaic**:
Component-Based Population-Adjusted Indirect Comparison.

## What this package is

cpaic extends **component network meta-analysis (cNMA)** to
**population-adjusted indirect comparison (PAIC)** methods (STC, MAIC,
ML-NMR). The goal is to take a treatment network that is *disconnected*
(two or more sub-networks with no common comparator) and:

1.  **Reconnect** it using the additive component structure of cNMA —
    when sub-networks share treatment components, component effects
    bridge the gap (Rücker et al. 2020/2021).
2.  **Adjust** the bridged comparisons for between-study imbalance in
    effect modifiers using anchored STC, MAIC, or ML-NMR, where
    individual patient data (IPD) are available.

The result is an indirect comparison that is *both* connected and
population-adjusted. No existing package does this; it is novel
methodology. See `documentation/refs/` for the literature it builds on
and `documentation/PROGRESS.md` for the running build log.

## Two-layer architecture

                IPD trials            AgD trials
                    |                      |
       [adjustment layer: PAIC]           |
       cstc() / cmaic() / cmlnmr()        |
       -> effect-modifier-adjusted        |
          contrasts per edge              |
                    \                    /
                     v                  v
            [connection layer: component NMA]
            cnma_bridge()  (netmeta::discomb)
            component effects beta; theta = C %*% beta
                             |
                             v
            connected + population-adjusted network
            relative_effects(), marginal_effects(), predict()

- **Connection layer** = additive cNMA. Design decomposition
  `delta = X beta = B C beta`, where `B` is the edge-incidence
  (contrast) matrix and `C` is the treatment-by-component matrix.
  [`netmeta::discomb()`](https://rdrr.io/pkg/netmeta/man/discomb.html)
  estimates component effects across a disconnected network by weighted
  least squares `beta_hat = (X' W X)^+ X' W d`, with
  `Cov(theta_hat) = C (X' W X)^+ C'`. Identifiable iff
  `rank(X) = n_components`.
- **Adjustment layer** = PAIC. Each IPD-bearing edge is corrected for
  effect-modifier imbalance, producing a population-adjusted contrast
  `d` and variance that feed the connection layer.

## Method roadmap (phased)

- **Phase 1 (frequentist core)** —
  [`cnma_bridge()`](https://choxos.github.io/cpaic/reference/cnma_bridge.md),
  [`cmaic()`](https://choxos.github.io/cpaic/reference/cmaic.md),
  [`cstc()`](https://choxos.github.io/cpaic/reference/cstc.md).
- **Phase 2 (Bayesian flagship)** —
  [`cmlnmr()`](https://choxos.github.io/cpaic/reference/cmlnmr.md):
  component-additive ML-NMR (treatment effect = `C %*% beta_component`),
  integrating the individual model over each study’s covariate
  distribution. Disconnected sub-networks share component parameters and
  are therefore connected by construction. **Status: all four families
  implemented and validated** via `cmdstanr` —
  `inst/stan/cpaic_binomial.stan` (logit), `cpaic_normal.stan`
  (identity; exact at covariate means), `cpaic_poisson.stan` (log; also
  used for **survival** via the exponential / person-time formulation).
  Flexible parametric and spline survival baselines are future work
  (proportional-hazards survival is covered frequentially by
  [`cstc()`](https://choxos.github.io/cpaic/reference/cstc.md)/[`cmaic()`](https://choxos.github.io/cpaic/reference/cmaic.md)).
  Caveat: components with effect modifiers informed only by aggregate
  data are weakly identified (`prior_reg_sd` regularizes).
- **Phase 3** — validation (simulation + reproduction of published
  examples) and the documentation suite.

## Reused engines (do NOT reimplement)

| Need | Engine | Key functions |
|----|----|----|
| cNMA (connected + disconnected) | `netmeta` | `netcomb()`, `discomb()`, `netcomplex()` |
| Anchored MAIC weights | `maicplus` | `estimate_weights()`, `center_ipd()`, `bucher()`, `maic_anchored()` |
| Network data model + ML-NMR + integration | `multinma` | `set_ipd()`, `set_agd_arm()`, `set_agd_contrast()`, `combine_network()`, `add_integration()`, `nma()` |
| STC delta-method / G-computation (reference) | `mlumr` (Suggests) | `stc.R` internals — **adapted, not called** |

### Important: STC is written natively here

[`mlumr::stc()`](https://choxos.github.io/mlumr/reference/stc.html) is
**unanchored** (single index-vs-comparator pair, no common comparator).
cpaic targets **anchored/connected** comparisons, so
[`cstc()`](https://choxos.github.io/cpaic/reference/cstc.md) is
implemented in-package: fit an outcome regression on IPD with
treatment + covariate main effects + treatment×effect-modifier
interactions, predict the relative effect in the comparator (AgD)
population, then anchor (Bucher) against the AgD trial’s own
within-trial contrast. The delta-method and link machinery are adapted
from mlumr’s `stc.R`/`link_functions.R` with attribution. This is why
`mlumr` is in `Suggests`, not `Imports`.

## Source layout

    R/
      cpaic-package.R    # _PACKAGE doc, namespace imports
      data_setup.R       # cpaic_network(): wrap multinma data model + component coding
      connectivity.R     # cpaic_connectivity(): disconnection + identifiability (rank X = c)
      bridge.R           # cnma_bridge(): discomb-based reconstruction
      cmaic.R            # cmaic(): per-edge maicplus weights -> bridge
      cstc.R             # cstc(): native anchored STC -> bridge
      cmlnmr.R           # Phase 2: component-additive ML-NMR
      effects.R          # relative_effects/component_effects/marginal_effects/predict + S3
      diagnostics.R      # additivity Q-test, ESS, integration checks
      plot.R             # component network plots, forest plots
      utils.R, zzz.R
    inst/stan/           # Phase 2 component-additive ML-NMR models (compiled via cmdstanr at runtime)
    tests/testthat/
    data/, data-raw/     # example datasets (incl. constructed disconnected+IPD example)
    vignettes/           # user-facing .Rmd
    documentation/       # GITIGNORED. Technical docs, refs, validation, progress log.

## Outcome families

All four: **binary** (OR/RR/RD), **continuous/normal** (MD), **count/
Poisson** (RR), **survival/TTE** (HR, RMST). Note `maicplus` natively
covers only binary + TTE, so
[`cmaic()`](https://choxos.github.io/cpaic/reference/cmaic.md) for
continuous/count uses the weighted-GLM path directly rather than
`maic_anchored()`.

## Dev commands

``` r

devtools::load_all()                 # load package
devtools::document()                 # roxygen -> NAMESPACE + man/
devtools::test()                     # run testthat
devtools::check()                    # R CMD check
pkgdown::build_site()                # docs site
# Phase 2 Stan models compile lazily with cmdstanr; cache under tools::R_user_dir.
```

Quarto math doc:
`quarto render documentation/mathematical-foundations/cpaic-math.qmd`.

## Conventions

- roxygen2 (markdown), testthat edition 3, S3 classes (`cpaic_network`,
  `cpaic_fit`, `cpaic_bridge`, …) mirroring multinma/mlumr style.
- Distribution: CRAN-eligible (all four engine deps are on CRAN).
- American English spelling; no dash punctuation as a connector (see
  user’s global style rules). GPL-3.

## Key references (in documentation/refs/)

- **Rücker 2020/2021** — additive cNMA; cNMA vs matching in a
  disconnected network (the conceptual seed; code in appendices).
- **Wigle & Béliveau 2022** — unanchored additive cNMA (data-driven
  anchor; no anchor-misspecification bias).
- **Efthimiou 2022** — Bayesian AD+IPD cNMA with component×covariate
  interactions (closest prior art for the adjustment layer).
- **Veroniki 2026** — guidance/decision tree for cNMA on disconnected
  networks (the additivity-across-subnetworks assumption).
- **Phillippo et al. 2020** — ML-NMR (the basis for
  [`cmlnmr()`](https://choxos.github.io/cpaic/reference/cmlnmr.md)).
