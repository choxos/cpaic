# Component simulated treatment comparison (cSTC)

Anchored STC generalized to a (possibly disconnected) component network.
For each IPD study an outcome regression is fitted with treatment main
effects, prognostic main effects, and treatment-by-effect-modifier
interactions. The effect modifiers are centered at common target means,
so the treatment coefficient is the anchored average conditional
link-scale contrast at those means. These adjusted contrasts replace the
corresponding unadjusted aggregate contrasts and
[`cnma_bridge()`](https://choxos.github.io/cpaic/reference/cnma_bridge.md)
combines them through the additive component model.

## Usage

``` r
cstc(
  network,
  target,
  effect_modifiers = NULL,
  prognostics = NULL,
  common = FALSE,
  random = TRUE,
  reference = NULL,
  allow_experimental_bridge = FALSE,
  allow_ipd_only_studies = FALSE
)
```

## Arguments

- network:

  A
  [`cpaic_network()`](https://choxos.github.io/cpaic/reference/cpaic_network.md)
  object that includes IPD.

- target:

  Named numeric vector (or list / one-row data frame) of target means
  for the effect modifiers.

- effect_modifiers:

  Covariates that interact with treatment (centered at `target`).
  Defaults to all IPD covariates.

- prognostics:

  Covariates included as main effects only. Defaults to the effect
  modifiers (so each enters as main effect + interaction).

- common, random:

  Passed to
  [`cnma_bridge()`](https://choxos.github.io/cpaic/reference/cnma_bridge.md).

- reference:

  Optional anchor (comparator) arm to use in every IPD study in which it
  appears, instead of inferring it from the aggregate row order.

- allow_experimental_bridge:

  Logical. The default `FALSE` stops when aggregate-only edges would be
  combined with target-adjusted IPD edges. Set `TRUE` only for
  explicitly exploratory sensitivity work; the fit records the retained
  edges and the reason the bridge is approximate.

- allow_ipd_only_studies:

  Logical. The default `FALSE` requires every IPD study to match exactly
  one aggregate two-arm edge. Set `TRUE` to append an IPD-derived edge
  that has no aggregate row. Such additions are recorded in the returned
  fit.

## Value

An object of class `cpaic_stc` (and `cpaic_fit`).

## Details

Unlike [`cmaic()`](https://choxos.github.io/cpaic/reference/cmaic.md)
(reweighting) this is the regression-adjustment route. The reported
treatment coefficient is the *conditional* effect at the target
effect-modifier means. Equivalently, under the fitted linear interaction
model this is the average conditional link-scale effect at the supplied
target means, not a marginal standardization. It is implemented natively
here because the `stc()` function in the mlumr package targets the
*unanchored* two-trial case; the link and standard-error machinery is
adapted from that package. (Written without the double-colon form on
purpose: mlumr is not a dependency of cpaic, and the documentation site
resolves a qualified package-and-function reference by loading that
package.)

## What the two-stage bridge does and does not adjust

Only the edges carrying individual patient data are population-adjusted
to the target means. Every aggregate-only edge keeps its published
study-specific contrast, and the additive bridge then combines all edges
as if they estimated the same component effects. Under effect
modification they do not: an aggregate edge estimates its contrast in
*its own* trial population, while the adjusted IPD edge estimates it at
the target. The two agree only when the aggregate populations resemble
the target, or when the components on those edges are not
effect-modified. Treat a cross-network contrast that leans on
aggregate-only edges as adjusted for the IPD part alone. Prefer
[`cmlnmr()`](https://choxos.github.io/cpaic/reference/cmlnmr.md) for a
joint model whose average conditional link-scale outputs are explicitly
evaluated at common target effect-modifier means.

## See also

[`cmaic()`](https://choxos.github.io/cpaic/reference/cmaic.md),
[`cnma_bridge()`](https://choxos.github.io/cpaic/reference/cnma_bridge.md)

## Examples

``` r
net <- cpaic_network(cpaic_bin_agd, ipd = cpaic_bin_ipd, sm = "OR",
                     family = "binomial", ipd_covariates = "x1",
                     inactive = "Placebo")
fit <- cstc(net, target = c(x1 = 0), effect_modifiers = "x1",
             allow_experimental_bridge = TRUE)
#> Warning: cstc() cannot form a decision-grade component bridge:
#>   - retained aggregate-only edge(s) remain in their own study populations: S1: A vs Placebo; S2: B vs Placebo; S5: A+B+C vs A+B+D
#> Use cmlnmr() for a joint model, restrict the analysis to a design in which every edge is adjusted and the estimand is additive, or set `allow_experimental_bridge = TRUE` only for explicitly exploratory sensitivity work.
relative_effects(fit)
#> Relative effects (OR, natural scale)
#>  treatment comparator estimate estimate_link se_link lower  upper   scale     z
#>          A    Placebo    1.649         0.500   0.256 0.998  2.725 natural 1.951
#>        A+B    Placebo    2.460         0.900   0.363 1.209  5.005 natural 2.483
#>      A+B+C    Placebo    4.014         1.390   0.435 1.711  9.416 natural 3.194
#>      A+B+D    Placebo    4.669         1.541   0.430 2.009 10.850 natural 3.582
#>          B    Placebo    1.492         0.400   0.256 0.903  2.466 natural 1.560
#>      p
#>  0.051
#>  0.013
#>  0.001
#>  0.000
#>  0.119
#>   `se_link` is on the log-ratio scale; the interval is back-transformed.
additivity_test(fit)
#> Additive component model: fit statistics
#>   Total lack of fit (Q.additive): Q = 2.669, df = 1, p = 0.102
#>   Additivity restrictions (Q.diff): not available; no standard NMA
#>     is estimable on a disconnected network.
#>   Note: neither statistic tests whether component effects are constant
#>   ACROSS sub-networks, which is the assumption that bridges the gap.
#>   That assumption is untestable from the data and must be justified
#>   clinically.
```
