# Leverage plot

Ported from
[`multinma::plot.nma_dic()`](https://dmphillippo.github.io/multinma/reference/plot.nma_dic.html)
with `type = "leverage"` (Phillippo et al. 2020). Each data point's
leverage (its contribution to the effective number of parameters) is
plotted against its signed square root residual deviance, with contours
of constant DIC contribution. Points outside the `DIC = 3` contour are
usually the ones spoiling the fit.

## Usage

``` r
plot_leverage(object, ..., dic_contours = 1:4)
```

## Arguments

- object:

  A [`cmlnmr()`](https://choxos.github.io/cpaic/reference/cmlnmr.md)
  fit.

- ...:

  Unused.

- dic_contours:

  Numeric vector of DIC contours to draw. Default `1:4`.

## Value

A `ggplot` object.

## Details

Leverage is the pointwise `pV` penalty, half the posterior variance of
the point's deviance, matching the penalty used by
[`dic()`](https://choxos.github.io/cpaic/reference/dic.md). Because
deviances covary across points, the pointwise leverages need not sum
exactly to the model's total `pV`.

Leverage plots need a saturated model and so are **not available for
survival outcomes**, where the censored contributions have no saturated
reference. multinma declines them for the same reason.

## See also

[`dic()`](https://choxos.github.io/cpaic/reference/dic.md),
[`plot.cpaic_dic()`](https://choxos.github.io/cpaic/reference/plot.cpaic_dic.md)

## Examples

``` r
if (FALSE) {
plot_leverage(fit)
}
```
