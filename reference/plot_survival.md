# Fitted survival curves from a cML-NMR fit

Ported from
[`multinma::plot.surv_nma_summary()`](https://dmphillippo.github.io/multinma/reference/plot.nma_summary.html)
(Phillippo et al. 2020). Draws the model-implied survival function for
each study arm, averaged over that arm's own covariate distribution (the
IPD covariates for an IPD arm, the integration points for an aggregate
arm), with posterior credible bands. Add
[`geom_km()`](https://choxos.github.io/cpaic/reference/geom_km.md) to
overlay the observed Kaplan-Meier curves; systematic departures indicate
that the baseline hazard is too rigid.

## Usage

``` r
plot_survival(object, ..., times = NULL, ndraws = 200L, level = 0.95)
```

## Arguments

- object:

  A [`cmlnmr()`](https://choxos.github.io/cpaic/reference/cmlnmr.md) fit
  with `family = "survival"`.

- ...:

  Unused.

- times:

  Time grid. Defaults to 100 points spanning the observed times.

- ndraws:

  Number of posterior draws to summarize over. Default `200`.

- level:

  Credible level for the band. Default `0.95`.

## Value

A `ggplot` object.

## Details

The baseline hazard is whatever the model was fitted with: a piecewise
exponential step function, or a continuous cubic M-spline. Its posterior
enters the curve through the integrated basis, so the bands include
baseline uncertainty.

## See also

[`geom_km()`](https://choxos.github.io/cpaic/reference/geom_km.md),
[`cmlnmr()`](https://choxos.github.io/cpaic/reference/cmlnmr.md)

## Examples

``` r
if (FALSE) {
plot_survival(fit) + geom_km(fit)
}
```
