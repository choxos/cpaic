# Pareto-smoothed importance sampling leave-one-out cross-validation

This is **observation-level** LOO: it measures within-study
interpolation (leaving out one IPD patient or one reconstructed
pseudo-observation), which is not the scientific question in a
disconnected network. It does **not** validate cross-gap prediction (a
new study, a held-out treatment contrast, or a held-out sub-network); a
good pointwise LOO can coexist with a wrong cross-gap extrapolation.
Grouped leave-one-study-out is not yet implemented, so `unit` accepts
only `"observation"`.

## Usage

``` r
# S3 method for class 'cpaic_mlnmr'
loo(x, unit = "observation", ...)
```

## Arguments

- x:

  A [`cmlnmr()`](https://choxos.github.io/cpaic/reference/cmlnmr.md)
  fit.

- unit:

  Predictive unit; only `"observation"` is supported.

- ...:

  Passed to
  [`loo::loo.matrix()`](https://mc-stan.org/loo/reference/loo.html).

## Value

A `psis_loo` object from the `loo` package.
