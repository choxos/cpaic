# Fit statistics for the additive component model

Returns the Cochran Q statistics from
[`netmeta::discomb()`](https://rdrr.io/pkg/netmeta/man/discomb.html).

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

A one-row data frame with `Q`, `df`, `pval` (total additive-model lack
of fit) and `Q.diff`, `df.diff`, `pval.diff` (the nested additivity
test, `NA` on a disconnected network).

## Details

Two different statistics are reported, and only one of them tests
additivity:

- `Q` (`Q.additive`) is the **total lack of fit** of the additive
  component model, pooling ordinary heterogeneity/inconsistency with any
  failure of additivity. It is not a test of additivity.

- `Q.diff = Q.additive - Q.standard` is the **nested test of the
  additivity restrictions** themselves, and is the statistic to read. It
  exists only when a standard (non-additive) NMA is also estimable, i.e.
  on a connected network; on a disconnected network it is `NA`.

**Neither statistic can test the assumption that actually bridges a
disconnected network**, namely that component effects (and, under
population adjustment, component x effect-modifier interactions) are
*constant across sub-networks*. There is by construction no cross-gap
evidence against which to test it. A large p-value here is therefore not
a licence to bridge; the assumption must be defended on clinical grounds
(Veroniki et al. 2026).
