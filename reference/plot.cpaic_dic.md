# Deviance and dev-dev plots

Ported from
[`multinma::plot.nma_dic()`](https://dmphillippo.github.io/multinma/reference/plot.nma_dic.html)
(Phillippo et al. 2020). With a single
[`dic()`](https://choxos.github.io/cpaic/reference/dic.md) object the
plot shows each data point's contribution to the posterior mean
deviance; points contributing much more than the rest are fitted poorly.
With two [`dic()`](https://choxos.github.io/cpaic/reference/dic.md)
objects it draws the **dev-dev plot**: points below the line of equality
fit better under the second model, points above it fit better under the
first.

## Usage

``` r
# S3 method for class 'cpaic_dic'
plot(x, y = NULL, ..., labels = c("Model 1", "Model 2"))
```

## Arguments

- x:

  A `cpaic_dic` object from
  [`dic()`](https://choxos.github.io/cpaic/reference/dic.md).

- y:

  An optional second `cpaic_dic` object, for a dev-dev plot.

- ...:

  Unused.

- labels:

  Model names used on the axes of a dev-dev plot.

## Value

A `ggplot` object.

## Details

[`dic()`](https://choxos.github.io/cpaic/reference/dic.md) stores the
posterior mean deviance per data point, not the residual deviance. The
two differ by the saturated log likelihood, which for the binomial and
Poisson families is a function of the data alone: it is the same under
both models, so it shifts both axes equally and the line of equality
keeps its meaning. For a Gaussian likelihood the saturated term also
involves the model's own `sigma`, so two Gaussian models with very
different residual variance are shifted by different amounts; read that
comparison with care.

For posterior uncertainty and for the leverage plot, which need the
saturated model explicitly, call
[`plot_leverage()`](https://choxos.github.io/cpaic/reference/plot_leverage.md)
on the fitted model itself.

## See also

[`dic()`](https://choxos.github.io/cpaic/reference/dic.md),
[`plot_leverage()`](https://choxos.github.io/cpaic/reference/plot_leverage.md)

## Examples

``` r
if (FALSE) {
plot(dic(fit_fe), dic(fit_re), labels = c("Fixed", "Random"))
}
```
