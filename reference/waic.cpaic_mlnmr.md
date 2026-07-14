# Widely applicable information criterion

Widely applicable information criterion

## Usage

``` r
# S3 method for class 'cpaic_mlnmr'
waic(x, ...)
```

## Arguments

- x:

  A [`cmlnmr()`](https://choxos.github.io/cpaic/reference/cmlnmr.md)
  fit.

- ...:

  Passed to
  [`loo::waic.matrix()`](https://mc-stan.org/loo/reference/waic.html).

## Value

A `waic` object from the `loo` package.
