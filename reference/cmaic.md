# Component matching-adjusted indirect comparison (cMAIC)

Anchored MAIC generalized to a (possibly disconnected) component
network. Each IPD study is reweighted with
[`maicplus::estimate_weights()`](https://hta-pharma.github.io/maicplus/main/reference/estimate_weights.html)
so that its effect-modifier distribution matches a common `target`
population; the resulting population-adjusted within-study contrasts
(with bootstrap standard errors that propagate the weighting
uncertainty) then replace the corresponding unadjusted aggregate
contrasts. Finally
[`cnma_bridge()`](https://choxos.github.io/cpaic/reference/cnma_bridge.md)
combines all contrasts through the additive component model, yielding
relative effects that are both connected across sub-networks and
adjusted to the target population.

## Usage

``` r
cmaic(
  network,
  target,
  effect_modifiers = NULL,
  target_sd = NULL,
  n_boot = 500,
  seed = NULL,
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

  Named numeric vector (or one-row data frame / list) giving the
  target-population means of the effect modifiers.

- effect_modifiers:

  Character vector of covariates to match on (defaults to all IPD
  covariates). Matching only on effect modifiers is the anchored-MAIC
  convention.

- target_sd:

  Optional named numeric vector of target standard deviations; when
  supplied, second moments are matched as well.

- n_boot:

  Number of bootstrap resamples for the adjusted-contrast standard
  errors. Default `500`.

- seed:

  Optional RNG seed for reproducible bootstrap.

- common, random:

  Passed to
  [`cnma_bridge()`](https://choxos.github.io/cpaic/reference/cnma_bridge.md).

## Value

An object of class `cpaic_maic` (also inheriting `cpaic_bridge`
structure via `$bridge`), with the bridged fit, per-study effective
sample sizes, and the target population.

## Non-collapsibility and the additive model

cMAIC returns a **marginal** effect in the target population, and the
additive component model assumes effects add. On a non-collapsible scale
(the odds ratio, the hazard ratio) **marginal effects do not add**, even
when every conditional effect does. In one simulated target population
the marginal log-odds ratios satisfied
`marginal(A) + marginal(B) = 0.6615` while `marginal(A+B) = 0.6411`; the
additive model is simply false on that scale. cMAIC therefore carries a
small **irreducible bias** (about +0.02 log-OR there) that survives
perfect matching and infinite sample size. It is small relative to a
typical standard error (about 0.25) but it does not vanish with more
data.

Marginal component effects are not *generally* additive; they add
exactly when the standardized treatment effects remain affine in the
component design. Additivity is therefore a property of the conditional
link scale that the marginal scale inherits only approximately, and the
error does not vanish with sample size. Where it is material,
[`cstc()`](https://choxos.github.io/cpaic/reference/cstc.md) or
[`cmlnmr()`](https://choxos.github.io/cpaic/reference/cmlnmr.md), which
target a conditional effect and inherit additivity exactly, are
preferable. Note also that the two-stage route combines a conditional
adjusted edge with aggregate edges reported on a marginal scale, so it
should be regarded as approximate. See
`documentation/validation/VALIDATION.md`.

## See also

[`cstc()`](https://choxos.github.io/cpaic/reference/cstc.md),
[`cnma_bridge()`](https://choxos.github.io/cpaic/reference/cnma_bridge.md)

## Examples

``` r
net <- cpaic_network(cpaic_bin_agd, ipd = cpaic_bin_ipd, sm = "OR",
                     family = "binomial", ipd_covariates = "x1",
                     inactive = "Placebo")
# \donttest{
fit <- cmaic(net, target = c(x1 = 0), effect_modifiers = "x1",
             n_boot = 100, seed = 1)
relative_effects(fit)
#> Relative effects (OR, back-transformed)
#>  treatment comparator estimate    se lower  upper     z     p
#>          A    Placebo    1.649 0.401 0.752  3.615 1.248 0.212
#>        A+B    Placebo    2.460 0.567 0.810  7.466 1.589 0.112
#>      A+B+C    Placebo    4.941 0.672 1.323 18.448 2.377 0.017
#>      A+B+D    Placebo    5.324 0.666 1.443 19.647 2.510 0.012
#>          B    Placebo    1.492 0.401 0.680  3.271 0.999 0.318
effective_sample_size(fit)
#>       S3       S4 
#> 207.4202 358.1461 
# }
```
