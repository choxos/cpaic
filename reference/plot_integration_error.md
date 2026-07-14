# Numerical integration error against the number of integration points

Ported from
[`multinma::plot_integration_error()`](https://dmphillippo.github.io/multinma/reference/plot_integration_error.html)
(Phillippo et al. 2020). Aggregate arms are fitted by integrating the
individual-level model over the study's covariate distribution with
Sobol' quasi-Monte-Carlo points. The integration error at `N` points is
the estimate using the first `N` points minus the estimate using all of
them; the typical convergence rate of QMC integration, `1/N`, is drawn
for reference. If the error has not settled well inside the `1/N`
envelope by `n_int`, refit with more integration points.

## Usage

``` r
plot_integration_error(
  x,
  ...,
  int_thin = NULL,
  ndraws = 200L,
  show_expected_rate = TRUE
)
```

## Arguments

- x:

  A [`cmlnmr()`](https://choxos.github.io/cpaic/reference/cmlnmr.md)
  fit.

- ...:

  Unused.

- int_thin:

  Report the error every `int_thin` points. Default is `n_int / 8`,
  rounded up.

- ndraws:

  Number of posterior draws to summarize over. Default `200`.

- show_expected_rate:

  Draw the `1/N` convergence envelope? Default `TRUE`.

## Value

A `ggplot` object.

## Details

cpaic does not save cumulative integration points inside Stan, so the
integrated aggregate-arm quantity is reconstructed here from the
posterior draws and the (deterministic) Sobol' point set. This is exact,
not an approximation, but it is not free: subsample the posterior with
`ndraws` on a large fit.

Not available for `family = "survival"`, where the aggregate
contribution is a `log_sum_exp` over integration points of the
likelihood rather than an integrated mean outcome; multinma declines
this plot for survival models too. Nor for a Gaussian model fitted with
all-normal margins, which is *exact* at the covariate means and uses a
single integration point.

## See also

[`cmlnmr()`](https://choxos.github.io/cpaic/reference/cmlnmr.md)

## Examples

``` r
if (FALSE) {
plot_integration_error(fit)
}
```
