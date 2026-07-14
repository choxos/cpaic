# Plot a population-adjusted hierarchy

Plots the ranking metrics of a
[`cpaic_ranks()`](https://choxos.github.io/cpaic/reference/cpaic_ranks.md)
hierarchy. Ranking metrics depend on the set being ranked, so read them
alongside the relative effects, never instead of them.

## Usage

``` r
# S3 method for class 'cpaic_ranks'
plot(x, y, ..., metric = c("sucra", "mean_rank", "median_rank", "p_best"))
```

## Arguments

- x:

  A `cpaic_ranks` object from
  [`cpaic_ranks()`](https://choxos.github.io/cpaic/reference/cpaic_ranks.md).

- y:

  Unused, for compatibility with the
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) generic.

- ...:

  Unused.

- metric:

  Which metric to plot: `"sucra"` (default), `"mean_rank"`,
  `"median_rank"`, or `"p_best"`.

## Value

A `ggplot` object.

## See also

[`cpaic_ranks()`](https://choxos.github.io/cpaic/reference/cpaic_ranks.md),
[`plot_rank_curve()`](https://choxos.github.io/cpaic/reference/plot_rank_curve.md)

## Examples

``` r
if (FALSE) {
plot(cpaic_ranks(fit, newdata = data.frame(x1 = 0)))
}
```
