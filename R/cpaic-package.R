#' @keywords internal
"_PACKAGE"

## cpaic: Component-Based Population-Adjusted Indirect Comparison
##
## The package is organized in two layers (see CLAUDE.md and
## documentation/PROGRESS.md for the full architecture):
##
##   1. Connection layer (component NMA). The additive component model of
##      Ruecker et al. (2020) decomposes multi-component treatments into
##      additive component effects. When disconnected sub-networks share
##      components, the component effects bridge the gap. Implemented in
##      `cnma_bridge()` on top of `netmeta::discomb()`.
##
##   2. Adjustment layer (population-adjusted indirect comparison). Where
##      individual patient data (IPD) are available, each evidence edge is
##      adjusted for effect-modifier imbalance with anchored STC
##      (`cstc()`), anchored MAIC (`cmaic()`, via `maicplus`), or
##      component-additive ML-NMR (`cmlnmr()`, Phase 2). The adjusted
##      contrasts are then combined through the component model.
##
## Anchored STC here is implemented natively (not via `mlumr::stc()`,
## which targets the unanchored two-trial case); the delta-method and
## G-computation machinery is adapted from 'mlumr' with attribution.

## usethis namespace: start
#' @importFrom stats coef vcov glm pnorm qnorm as.formula relevel sd setNames binomial gaussian poisson
## usethis namespace: end
NULL
