#' @keywords internal
"_PACKAGE"

## cpaic: Component-Based Population-Adjusted Indirect Comparison
##
## This package is experimental research software. Its methodology and
## implementation have not been validated for clinical, regulatory,
## reimbursement, or other decision use.
##
## The package is organized in two layers:
##
##   1. Connection layer (component NMA). The additive component model of
##      Ruecker et al. (2020) decomposes multi-component treatments into
##      additive component effects. When disconnected sub-networks share
##      components, the component effects bridge the gap. Implemented in
##      `cnma_bridge()` on top of `netmeta::discomb()`.
##
##   2. Adjustment layer (population-adjusted indirect comparison). Where
##      individual patient data (IPD) are available, an evidence edge can be
##      adjusted for effect-modifier imbalance with anchored STC
##      (`cstc()`), anchored MAIC (`cmaic()`, via `maicplus`), or modeled
##      jointly with component-additive ML-NMR (`cmlnmr()`).
##
## The two-stage `cstc()` and `cmaic()` paths are not automatically a coherent
## whole-network population adjustment. If aggregate-data edges remain, their
## published contrasts refer to their own study populations, while adjusted
## IPD edges refer to the supplied target. Further, marginal contrasts from
## `cmaic()` are not generally component-additive on nonlinear link scales.
## Affected bridges stop unless the caller explicitly sets
## `allow_experimental_bridge = TRUE`.
##
## Anchored STC is implemented natively. It fits an outcome regression and
## evaluates the conditional link-scale contrast at supplied target covariate
## means. The cML-NMR reporting methods likewise return average conditional
## link-scale effects evaluated at supplied covariate means. They do not
## marginally standardize odds, rate, or hazard ratios over a target covariate
## distribution.
##
## The Stan models in inst/stan are compiled at install time through the
## rstantools scaffolding (configure -> rstantools::rstan_config()), so
## cmlnmr() works with backend = "rstan" out of the box. backend =
## "cmdstanr" fits the same models with CmdStan when it is available.
## Everything downstream of a fit goes through the accessors in
## R/backend.R rather than either engine's API directly.

## usethis namespace: start
#' @importFrom stats coef vcov glm pnorm qnorm as.formula relevel sd setNames binomial gaussian poisson
#' @import methods
#' @import Rcpp
#' @importFrom rstan sampling
#' @importFrom rstantools rstan_config
#' @importFrom RcppParallel RcppParallelLibs
#' @useDynLib cpaic, .registration = TRUE
## usethis namespace: end
NULL
