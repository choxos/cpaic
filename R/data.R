# Bundled example data ---------------------------------------------------------

#' Example disconnected component network: aggregate contrasts
#'
#' A constructed binary-outcome network used throughout the examples and
#' tests. It is **disconnected**: sub-network 1 (`Placebo`, `A`, `B`) is
#' anchored on placebo, while sub-network 2 (`A+B`, `A+B+C`, `A+B+D`) shares
#' no treatment with sub-network 1. The shared components `A` and `B` bridge
#' the two, and components `C` and `D` are identified within sub-network 2.
#'
#' Studies `S3` and `S4` also have individual patient data
#' ([cpaic_bin_ipd]); their rows here are the *unadjusted* contrasts, which
#' [cmaic()] / [cstc()] replace with population-adjusted versions.
#'
#' @format A data frame with 5 rows and 5 columns:
#' \describe{
#'   \item{studlab}{Study label.}
#'   \item{treat1, treat2}{Treatments compared (components joined by `"+"`).}
#'   \item{TE}{Log odds ratio of `treat1` versus `treat2`.}
#'   \item{seTE}{Standard error of `TE`.}
#' }
#' The attribute `"truth"` holds the data-generating component log-odds
#' ratios.
#' @seealso [cpaic_bin_ipd], [cpaic_network()]
"cpaic_bin_agd"

#' Example disconnected component network: individual patient data
#'
#' Individual patient data for studies `S3` (`A+B+C` vs `A+B`) and `S4`
#' (`A+B+D` vs `A+B`) of the [cpaic_bin_agd] network. A single effect
#' modifier `x1` is imbalanced relative to the target population
#' (`x1 = 0`), so population adjustment changes the `C` and `D` component
#' effects.
#'
#' @format A data frame with 3200 rows and 4 columns:
#' \describe{
#'   \item{.study}{Study label (`S3` or `S4`).}
#'   \item{.trt}{Treatment arm.}
#'   \item{.y}{Binary outcome (0/1).}
#'   \item{x1}{Continuous effect modifier.}
#' }
#' The attribute `"truth"` holds the data-generating parameters.
#' @seealso [cpaic_bin_agd], [cmaic()], [cstc()]
"cpaic_bin_ipd"
