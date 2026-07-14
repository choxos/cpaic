# Deviance information criterion

Computes DIC with the variance penalty `pV`, following the
survival-model implementation in multinma. The pointwise deviance is
`-2 * log_lik`; the effective parameter count is half the posterior
variance of total deviance.

## Usage

``` r
dic(x, ...)
```

## Arguments

- x:

  A fitted model.

- ...:

  Unused.

## Value

For [`cmlnmr()`](https://choxos.github.io/cpaic/reference/cmlnmr.md)
fits, a `cpaic_dic` object with DIC, mean deviance, the `pV` penalty,
and pointwise mean deviance.
