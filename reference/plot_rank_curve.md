# How the hierarchy changes across target effect-modifier means

The component effects are `beta + Gamma x`, so a component's rank is a
function of the target means `x` and components can cross. A single
hierarchy quoted without its target means is therefore incomplete. This
plot shows the family of average conditional link-scale hierarchies
across a mean grid. It does not standardize over different population
distributions.

## Usage

``` r
plot_rank_curve(
  x,
  em = NULL,
  values = NULL,
  at = NULL,
  what = c("treatment", "component"),
  lower_is_better = FALSE,
  metric = c("sucra", "mean_rank", "p_best"),
  ...
)
```

## Arguments

- x:

  A [`cmlnmr()`](https://choxos.github.io/cpaic/reference/cmlnmr.md)
  fit, or the data frame returned by
  [`rank_curve()`](https://choxos.github.io/cpaic/reference/rank_curve.md).

- em:

  Name of the effect modifier to vary. Required when `x` is a fit.

- values:

  Numeric vector of target values for `em`. Required when `x` is a fit.

- at:

  Optional named vector fixing the other effect modifiers.

- what, lower_is_better:

  See
  [`cpaic_ranks()`](https://choxos.github.io/cpaic/reference/cpaic_ranks.md).

- metric:

  Which ranking metric to trace: `"sucra"` (default), `"mean_rank"`, or
  `"p_best"`.

- ...:

  Passed to
  [`rank_curve()`](https://choxos.github.io/cpaic/reference/rank_curve.md)
  when `x` is a fit.

## Value

A `ggplot` object.

## Details

There is no counterpart in multinma, which ranks in one population at a
time.

## References

Wigle A, Beliveau A, Nikolakopoulou A, Lin L (2026). Creating Treatment
and Component Hierarchies in Component Network Meta-Analysis.

## See also

[`rank_curve()`](https://choxos.github.io/cpaic/reference/rank_curve.md),
[`cpaic_ranks()`](https://choxos.github.io/cpaic/reference/cpaic_ranks.md),
[`plot_estimability()`](https://choxos.github.io/cpaic/reference/plot_estimability.md)

## Examples

``` r
if (FALSE) {
plot_rank_curve(fit, em = "x1", values = seq(-1, 1, by = 0.25),
                what = "component")
}
```
