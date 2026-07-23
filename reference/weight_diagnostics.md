# Weight-quality diagnostics for a cMAIC fit

Per IPD study, the effective sample size, weight-entropy efficiency,
coefficient of variation, largest normalized weight, mass in the top 5%
of weights, and the largest residual effect-modifier imbalance after
weighting. A high maximum weight or low entropy efficiency signals a few
dominant individuals, which the effective sample size alone can hide.

## Usage

``` r
weight_diagnostics(object)
```

## Arguments

- object:

  A [`cmaic()`](https://choxos.github.io/cpaic/reference/cmaic.md) fit.

## Value

A data frame, one row per IPD study.

## See also

[`cmaic()`](https://choxos.github.io/cpaic/reference/cmaic.md),
[`effective_sample_size()`](https://choxos.github.io/cpaic/reference/effective_sample_size.md)
