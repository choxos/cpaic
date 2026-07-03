# Component simulated treatment comparison (cSTC)

Anchored STC generalized to a (possibly disconnected) component network.
For each IPD study an outcome regression is fitted with treatment main
effects, prognostic main effects, and treatment-by-effect-modifier
interactions; the effect modifiers are centered at a common `target`
population so the treatment coefficient is the anchored,
population-adjusted contrast in that population. These adjusted
contrasts replace the corresponding unadjusted aggregate contrasts and
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
  random = TRUE
)
```

## Arguments

- network:

  A
  [`cpaic_network()`](https://choxos.github.io/cpaic/reference/cpaic_network.md)
  object that includes IPD.

- target:

  Named numeric vector (or list / one-row data frame) of
  target-population means for the effect modifiers.

- effect_modifiers:

  Covariates that interact with treatment (centered at `target`).
  Defaults to all IPD covariates.

- prognostics:

  Covariates included as main effects only. Defaults to the effect
  modifiers (so each enters as main effect + interaction).

- common, random:

  Passed to
  [`cnma_bridge()`](https://choxos.github.io/cpaic/reference/cnma_bridge.md).

## Value

An object of class `cpaic_stc` (and `cpaic_fit`).

## Details

Unlike [`cmaic()`](https://choxos.github.io/cpaic/reference/cmaic.md)
(reweighting) this is the regression-adjustment route. The reported
treatment coefficient is the *conditional* effect at the target
effect-modifier means (not a marginal standardization); for collapsible
measures the two coincide. It is implemented natively here because
[`mlumr::stc()`](https://choxos.github.io/mlumr/reference/stc.html)
targets the *unanchored* two-trial case; the link and standard-error
machinery is adapted from that package.

## See also

[`cmaic()`](https://choxos.github.io/cpaic/reference/cmaic.md),
[`cnma_bridge()`](https://choxos.github.io/cpaic/reference/cnma_bridge.md)

## Examples

``` r
net <- cpaic_network(cpaic_bin_agd, ipd = cpaic_bin_ipd, sm = "OR",
                     family = "binomial", ipd_covariates = "x1",
                     inactive = "Placebo")
fit <- cstc(net, target = c(x1 = 0), effect_modifiers = "x1")
relative_effects(fit)
#> Relative effects (OR, back-transformed)
#>  treatment comparator estimate    se lower  upper     z     p
#>          A    Placebo    1.649 0.256 0.998  2.725 1.951 0.051
#>        A+B    Placebo    2.460 0.363 1.209  5.005 2.483 0.013
#>      A+B+C    Placebo    4.014 0.435 1.711  9.416 3.194 0.001
#>      A+B+D    Placebo    4.669 0.430 2.009 10.850 3.582 0.000
#>          B    Placebo    1.492 0.256 0.903  2.466 1.560 0.119
additivity_test(fit)
#> Additivity (Cochran Q) test of the component model
#>   Q = 2.669, df = 1, p = 0.102
```
