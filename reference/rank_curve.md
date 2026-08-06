# How a hierarchy changes across target effect-modifier means

Recomputes
[`cpaic_ranks()`](https://choxos.github.io/cpaic/reference/cpaic_ranks.md)
over a grid of target means. This exposes how the hierarchy of average
conditional link-scale effects changes with the chosen mean. It does not
standardize effects over a sequence of target distributions.

## Usage

``` r
rank_curve(
  object,
  em,
  values,
  at = NULL,
  what = c("treatment", "component"),
  lower_is_better = FALSE,
  include_screen_only = FALSE,
  ...
)
```

## Arguments

- object:

  A [`cmlnmr()`](https://choxos.github.io/cpaic/reference/cmlnmr.md)
  fit.

- em:

  Name of the effect modifier to vary.

- values:

  Numeric vector of target values for `em`.

- at:

  Optional named vector fixing the other effect modifiers. Defaults to 0
  for each.

- what, lower_is_better, include_screen_only:

  See
  [`cpaic_ranks()`](https://choxos.github.io/cpaic/reference/cpaic_ranks.md).

- ...:

  Unused.

## Value

A data frame with one row per (element, target value), giving `sucra`,
`mean_rank`, `p_best`, and `estimate`. Failed target values are retained
as one `status = "failed"` row with `NA` metrics and an explanatory
`error`.

## See also

[`cpaic_ranks()`](https://choxos.github.io/cpaic/reference/cpaic_ranks.md)

## Examples

``` r
if (FALSE) {
rank_curve(fit, em = "x1", values = seq(-1, 1, by = 0.25), what = "component")
}
```
