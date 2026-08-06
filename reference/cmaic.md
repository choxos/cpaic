# Component matching-adjusted indirect comparison (cMAIC)

Anchored MAIC generalized to a (possibly disconnected) component
network. Each IPD study is reweighted with
[`maicplus::estimate_weights()`](https://hta-pharma.github.io/maicplus/main/reference/estimate_weights.html)
so that its requested effect-modifier moments match a common `target`;
the resulting target-matched within-study contrasts (with bootstrap
standard errors that propagate the weighting uncertainty) then replace
the corresponding unadjusted aggregate contrasts. Finally
[`cnma_bridge()`](https://choxos.github.io/cpaic/reference/cnma_bridge.md)
combines all contrasts through the additive component model. The bridge
is gated because retained aggregate edges and nonlinear marginal effects
can make that synthesis incoherent.

## Usage

``` r
cmaic(
  network,
  target,
  effect_modifiers = NULL,
  target_sd = NULL,
  n_boot = 500,
  min_boot_success = 0.8,
  seed = NULL,
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

  Named numeric vector (or one-row data frame / list) giving target
  means of the effect modifiers.

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

- min_boot_success:

  Minimum fraction of bootstrap resamples that must succeed for a
  contrast. The enforced count is
  `max(ceiling(min_boot_success * n_boot), min(20, n_boot))`, so a run
  with fewer than 20 requested resamples requires every resample to
  succeed. Below this threshold the edge is rejected rather than given a
  fragile standard error from a selected subset. Default `0.8`.

- seed:

  Optional RNG seed for reproducible bootstrap. The caller's global RNG
  state is restored on exit, so calling `cmaic()` does not perturb a
  downstream random stream.

- common, random:

  Passed to
  [`cnma_bridge()`](https://choxos.github.io/cpaic/reference/cnma_bridge.md).

- reference:

  Optional anchor (comparator) arm to use in every IPD study in which it
  appears, instead of inferring it from the aggregate row order.

- allow_experimental_bridge:

  Logical. The default `FALSE` stops when aggregate-only edges would be
  combined with target-matched IPD edges, or when a non-Gaussian cMAIC
  contrast would be forced through an additive component model. Set
  `TRUE` only for explicitly exploratory sensitivity work; the fit
  records the exact approximation reasons.

- allow_ipd_only_studies:

  Logical. The default `FALSE` requires every IPD study to match exactly
  one aggregate two-arm edge. Set `TRUE` to append an IPD-derived edge
  that has no aggregate row. Such additions are recorded in the returned
  fit.

## Value

An object of class `cpaic_maic` (also inheriting `cpaic_bridge`
structure via `$bridge`), with the bridged fit, per-study effective
sample sizes, and the target moments. Bootstrap diagnostic fields
include `$bootstrap_draws`, `$bootstrap_summary`, `$bootstrap_failures`,
`$bootstrap_failure_table`, `$bootstrap_mcse_method`, and
`$bootstrap_success_rule`. A threshold failure raises a
`cpaic_bootstrap_error` condition carrying the same diagnostic
information.

## What the two-stage bridge does and does not adjust

Only the edges carrying individual patient data are population-adjusted
to the target moments. Every aggregate-only edge keeps its published
study-specific contrast, and the additive bridge then combines all edges
as if they estimated the same component effects. Under effect
modification they do not: an aggregate edge estimates its contrast in
*its own* trial population, while the reweighted IPD edge estimates it
at the target. The two agree only when the aggregate populations
resemble the target, or when the components on those edges are not
effect-modified. Treat a cross-network contrast that leans on
aggregate-only edges as adjusted for the IPD part alone. Prefer
[`cmlnmr()`](https://choxos.github.io/cpaic/reference/cmlnmr.md) for a
joint model whose average conditional link-scale outputs are explicitly
evaluated at common target effect-modifier means.

## Non-collapsibility and the additive model

cMAIC returns a **marginal** effect in the reweighted IPD sample, and
the additive component model assumes effects add. On a non-collapsible
scale (the odds ratio, the hazard ratio) **marginal effects do not
add**, even when every conditional effect does. In one simulated target
population the marginal log-odds ratios satisfied
`marginal(A) + marginal(B) = 0.6615` while `marginal(A+B) = 0.6411`; the
additive model is simply false on that scale. cMAIC therefore carries an
**irreducible approximation error** that survives perfect matching and
infinite sample size. Its size is problem-specific and cannot be assumed
negligible.

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
should be regarded as approximate.

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
             n_boot = 100, seed = 1,
             allow_experimental_bridge = TRUE)
#> Warning: cmaic() cannot form a decision-grade component bridge:
#>   - retained aggregate-only edge(s) remain in their own study populations: S1: A vs Placebo; S2: B vs Placebo; S5: A+B+C vs A+B+D
#>   - cMAIC estimates marginal binomial contrasts, which are not generally additive in the component design on a nonlinear link scale
#> Use cmlnmr() for a joint model, restrict the analysis to a design in which every edge is adjusted and the estimand is additive, or set `allow_experimental_bridge = TRUE` only for explicitly exploratory sensitivity work.
relative_effects(fit)
#> Relative effects (OR, back-transformed)
#>  treatment comparator estimate estimate_link se_link lower  upper   scale     z
#>          A    Placebo    1.649         0.500   0.401 0.752  3.615 natural 1.248
#>        A+B    Placebo    2.460         0.900   0.567 0.810  7.466 natural 1.589
#>      A+B+C    Placebo    4.941         1.597   0.672 1.323 18.448 natural 2.377
#>      A+B+D    Placebo    5.324         1.672   0.666 1.443 19.647 natural 2.510
#>          B    Placebo    1.492         0.400   0.401 0.680  3.271 natural 0.999
#>      p
#>  0.212
#>  0.112
#>  0.017
#>  0.012
#>  0.318
#>   `se_link` is on the link (log) scale; the interval is back-transformed.
effective_sample_size(fit)
#>       S3       S4 
#> 207.4202 358.1461 
# }
```
