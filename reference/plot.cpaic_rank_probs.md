# Rankogram and cumulative rank plot

Ported from
[`multinma::plot.nma_rank_probs()`](https://dmphillippo.github.io/multinma/reference/plot.nma_summary.html)
(Phillippo et al. 2020). The rankogram gives the posterior probability
of each rank; the cumulative version gives the probability of being
ranked among the best `k`, whose normalized area is SUCRA.

## Usage

``` r
# S3 method for class 'cpaic_rank_probs'
plot(x, y, ...)
```

## Arguments

- x:

  A `cpaic_rank_probs` object from
  [`rank_probs()`](https://choxos.github.io/cpaic/reference/rank_probs.md).

- y:

  Unused, for compatibility with the
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) generic.

- ...:

  Unused.

## Value

A `ggplot` object.

## Details

Both are computed **in a named target population**, because a
population-adjusted hierarchy is not population-free.

## See also

[`rank_probs()`](https://choxos.github.io/cpaic/reference/rank_probs.md),
[`plot_rank_curve()`](https://choxos.github.io/cpaic/reference/plot_rank_curve.md)

## Examples

``` r
if (FALSE) {
plot(rank_probs(fit, newdata = data.frame(x1 = 0), what = "component"))
}
```
