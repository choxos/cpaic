# Relative treatment effects from a cpaic fit

Tidies the relative effects of the fitted model: every treatment versus
a chosen reference, or all pairwise comparisons. Effects are reported on
the natural scale of the summary measure (e.g. odds ratios) unless
`backtransf = FALSE`.

## Usage

``` r
relative_effects(
  object,
  reference = NULL,
  all_contrasts = FALSE,
  backtransf = TRUE,
  level = 0.95,
  newdata = NULL,
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

- newdata:

  For [`cmlnmr()`](https://choxos.github.io/cpaic/reference/cmlnmr.md)
  fits: a one-row data frame giving the target population's
  effect-modifier values. Required when the model has effect modifiers.

- ...:

  Unused.

## Value

A data frame with columns `treatment`, `comparator`, `estimate`, `se`
(link scale), `lower`, `upper`, and `z`/`p` for frequentist fits. For
[`cmlnmr()`](https://choxos.github.io/cpaic/reference/cmlnmr.md)
(Bayesian) fits the intervals are credible intervals and the final
column is `pr_gt0`, the posterior probability that the effect (on the
link scale) exceeds zero, instead of `z`/`p`.

## Details

Relative effects that the component design cannot uniquely identify
(their contrast vector lies outside the row space of `X = B C`) are
returned as `NA` rather than as pseudoinverse or prior-driven artefacts.
See
[`estimable_effects()`](https://choxos.github.io/cpaic/reference/estimable_effects.md).

For [`cmlnmr()`](https://choxos.github.io/cpaic/reference/cmlnmr.md)
fits the model contains component x effect-modifier interactions, so
relative effects are **population-specific**:
`theta_t(x) = C_t' (beta + gamma x)`. You must name the target
population through `newdata`; there is no population-free relative
effect.
