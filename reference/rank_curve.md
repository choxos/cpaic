# How a hierarchy changes across target populations

Recomputes
[`cpaic_ranks()`](https://choxos.github.io/cpaic/reference/cpaic_ranks.md)
over a grid of target populations, so that the population dependence of
the hierarchy is visible. Under population adjustment a component's rank
is a function of the target, and this is the object that shows it.

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
`mean_rank`, `p_best` and `estimate`, plus `estimable`.

## See also

[`cpaic_ranks()`](https://choxos.github.io/cpaic/reference/cpaic_ranks.md)

## Examples

``` r
if (FALSE) {
rank_curve(fit, em = "x1", values = seq(-1, 1, by = 0.25), what = "component")
}
```
