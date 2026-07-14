# MCMC diagnostics for a cML-NMR fit

Ported from
[`multinma::plot.stan_nma()`](https://dmphillippo.github.io/multinma/reference/summary.stan_nma.html)
(Phillippo et al. 2020), which hands the posterior draws to `bayesplot`.
Diverging transitions, a high maximum `Rhat`, or a low effective sample
size all mean the posterior has not been explored, and nothing
downstream of it is trustworthy.

## Usage

``` r
# S3 method for class 'cpaic_mlnmr'
plot(
  x,
  y,
  ...,
  type = c("trace", "density", "hist", "pairs", "rhat", "neff"),
  pars = NULL
)
```

## Arguments

- x:

  A [`cmlnmr()`](https://choxos.github.io/cpaic/reference/cmlnmr.md)
  fit.

- y:

  Unused, for compatibility with the
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) generic.

- ...:

  Passed to the underlying `bayesplot` function.

- type:

  `"trace"` (default), `"density"`, `"hist"`, `"pairs"`, `"rhat"`, or
  `"neff"`.

- pars:

  Character vector of parameter names to show. Defaults to the component
  effects `beta`, the interactions `gamma`, and, for a random-effects
  model, the heterogeneity `tau`.

## Value

A `ggplot` object (`"pairs"` returns a `bayesplot` grid).

## See also

[`cmlnmr()`](https://choxos.github.io/cpaic/reference/cmlnmr.md),
[`plot_prior_posterior()`](https://choxos.github.io/cpaic/reference/plot_prior_posterior.md)

## Examples

``` r
if (FALSE) {
plot(fit, type = "trace")
plot(fit, type = "pairs", pars = c("beta[1]", "tau[1]"))
}
```
