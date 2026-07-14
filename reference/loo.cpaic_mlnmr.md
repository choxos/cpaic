# Pareto-smoothed importance sampling leave-one-out cross-validation

Pareto-smoothed importance sampling leave-one-out cross-validation

## Usage

``` r
# S3 method for class 'cpaic_mlnmr'
loo(x, ...)
```

## Arguments

- x:

  A [`cmlnmr()`](https://choxos.github.io/cpaic/reference/cmlnmr.md)
  fit.

- ...:

  Passed to
  [`loo::loo.matrix()`](https://mc-stan.org/loo/reference/loo.html).

## Value

A `psis_loo` object from the `loo` package.
