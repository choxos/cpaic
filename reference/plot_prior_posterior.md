# Prior versus posterior

Ported from
[`multinma::plot_prior_posterior()`](https://dmphillippo.github.io/multinma/reference/plot_prior_posterior.html)
(Phillippo et al. 2020). Posteriors are drawn as histograms, priors as
lines. Where a posterior simply reproduces its prior, the data carry no
information about that parameter, and any quantity that leans on it is
prior-driven rather than estimated. This is the visual counterpart of
[`prior_sensitivity()`](https://choxos.github.io/cpaic/reference/prior_sensitivity.md).

## Usage

``` r
plot_prior_posterior(x, ..., prior = NULL, bins = 40)
```

## Arguments

- x:

  A [`cmlnmr()`](https://choxos.github.io/cpaic/reference/cmlnmr.md)
  fit.

- ...:

  Unused.

- prior:

  Which priors to show. Any of `"intercept"` (`mu`), `"beta"` (component
  effects), `"regression"` (`breg`), `"gamma"` (component by
  effect-modifier interactions), and `"tau"` (heterogeneity, random
  effects only). Defaults to all that the model used.

- bins:

  Number of histogram bins for the posterior. Default `40`.

## Value

A `ggplot` object.

## Details

It matters most for the component by effect-modifier interactions
`gamma`: interactions informed only by aggregate arms are weakly
identified, and `prior_gamma_scale` then does real regularization.

## See also

[`prior_sensitivity()`](https://choxos.github.io/cpaic/reference/prior_sensitivity.md),
[`prior_predictive_check()`](https://choxos.github.io/cpaic/reference/prior_predictive_check.md)

## Examples

``` r
if (FALSE) {
plot_prior_posterior(fit, prior = "gamma")
}
```
