# Relative treatment effects from a cpaic fit

Tidies the (random- or common-effects) relative effects from the
component-bridged model: every treatment versus a chosen reference, or
all pairwise comparisons. Effects are reported on the natural scale of
the summary measure (e.g. odds ratios) unless `backtransf = FALSE`.

## Usage

``` r
relative_effects(
  object,
  reference = NULL,
  all_contrasts = FALSE,
  backtransf = TRUE,
  level = 0.95,
  ...
)
```

## Arguments

- object:

  A fitted cpaic object (`cpaic_bridge`, `cpaic_maic`, `cpaic_stc`, or
  `cpaic_mlnmr`).

- reference:

  Reference treatment. Defaults to the network reference.

- all_contrasts:

  If `TRUE`, return all pairwise comparisons instead of versus the
  reference.

- backtransf:

  If `TRUE` (default) back-transform log-scale measures (OR/RR/HR/...)
  by exponentiating.

- level:

  Confidence level for the intervals. Default `0.95`.

- ...:

  Unused.

## Value

A data frame with columns `treatment`, `comparator`, `estimate`, `se`
(link scale), `lower`, `upper`, and `z`/`p` for frequentist fits. For
[`cmlnmr()`](https://choxos.github.io/cpaic/reference/cmlnmr.md)
(Bayesian) fits the intervals are credible intervals and the final
column is `pr_gt0`, the posterior probability that the effect (on the
link scale) exceeds zero, instead of `z`/`p`.
