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
  estimand = NULL,
  target = NULL,
  weights = NULL,
  measure = NULL,
  baseline_study = NULL,
  times = NULL,
  random_effect = "population",
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
  fits: a one-row data frame giving target effect-modifier means.
  Required when the model has effect modifiers.

- estimand:

  For [`cmlnmr()`](https://choxos.github.io/cpaic/reference/cmlnmr.md)
  fits, either `"average_conditional_link"` or `"marginal"`. The latter
  delegates to
  [`marginal_effects()`](https://choxos.github.io/cpaic/reference/marginal_effects.md).

- target, weights, measure, baseline_study, times, random_effect:

  Arguments used only for `estimand = "marginal"`; see
  [`marginal_effects()`](https://choxos.github.io/cpaic/reference/marginal_effects.md).

- ...:

  Unused.

## Value

For frequentist and average-conditional outputs, a data frame with
columns `treatment`, `comparator`, `estimate`, `estimate_link`,
`se_link`, `lower`, `upper`, `scale`, and `z`/`p` or Bayesian `pr_gt0`.
`estimate_link` and `se_link` are on the model link scale.
`estimand = "marginal"` returns the schema documented by
[`marginal_effects()`](https://choxos.github.io/cpaic/reference/marginal_effects.md),
with canonical `estimate_contrast`, `se_contrast`, and `contrast_scale`
columns because natural differences are not model-link quantities.

## Details

Relative effects that the component design cannot uniquely identify
(their contrast vector lies outside the row space of `X = B C`) are
returned as `NA` rather than as pseudoinverse or prior-driven artefacts.
See
[`estimable_effects()`](https://choxos.github.io/cpaic/reference/estimable_effects.md).

For [`cmlnmr()`](https://choxos.github.io/cpaic/reference/cmlnmr.md)
fits, `estimand = "average_conditional_link"` reports component-additive
link-scale contrasts at supplied covariate means:
`theta_t(x) = C_t' (beta + gamma x)`. Because this expression is linear
in `x`, evaluating it at `E[X]` gives the average conditional link-scale
effect. Set `estimand = "marginal"` to delegate to
[`marginal_effects()`](https://choxos.github.io/cpaic/reference/marginal_effects.md)
and standardize treatment-specific outcomes over an explicit target
covariate distribution before forming contrasts.

The average-conditional path is narrower than it may look. A one-row
`newdata` profile is not a target covariate distribution, and
exponentiating a link-scale contrast does not create a marginal
standardized effect. The marginal path therefore requires `target`
instead. In either path, interactions informed only by aggregate arms
can reflect ecological rather than within-study effect modification; see
[`estimable_effects_at()`](https://choxos.github.io/cpaic/reference/estimable_effects_at.md).
