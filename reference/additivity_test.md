# Test the additivity assumption of the component model

Returns the Cochran Q statistic for the additive component network
meta-analysis (from
[`netmeta::discomb()`](https://rdrr.io/pkg/netmeta/man/discomb.html)). A
small p-value indicates lack of fit of the additive model, i.e. evidence
of component interactions; consider adding interaction terms.

## Usage

``` r
additivity_test(object, ...)
```

## Arguments

- object:

  A `cpaic_bridge` or `cpaic_fit` object.

- ...:

  Unused.

## Value

A one-row data frame with `Q`, `df`, and `pval`.
