# Example disconnected component network: aggregate contrasts

A constructed binary-outcome network used throughout the examples and
tests. It is **disconnected**: sub-network 1 (`Placebo`, `A`, `B`) is
anchored on placebo, while sub-network 2 (`A+B`, `A+B+C`, `A+B+D`)
shares no treatment with sub-network 1. The shared components `A` and
`B` bridge the two, and components `C` and `D` are identified within
sub-network 2.

## Usage

``` r
cpaic_bin_agd
```

## Format

A data frame with 5 rows and 5 columns:

- studlab:

  Study label.

- treat1, treat2:

  Treatments compared (components joined by `"+"`).

- TE:

  Log odds ratio of `treat1` versus `treat2`.

- seTE:

  Standard error of `TE`.

The attribute `"truth"` holds the data-generating component log-odds
ratios.

## Details

Studies `S3` and `S4` also have individual patient data
([cpaic_bin_ipd](https://choxos.github.io/cpaic/reference/cpaic_bin_ipd.md));
their rows here are the *unadjusted* contrasts, which
[`cmaic()`](https://choxos.github.io/cpaic/reference/cmaic.md) /
[`cstc()`](https://choxos.github.io/cpaic/reference/cstc.md) replace
with population-adjusted versions.

## See also

[cpaic_bin_ipd](https://choxos.github.io/cpaic/reference/cpaic_bin_ipd.md),
[`cpaic_network()`](https://choxos.github.io/cpaic/reference/cpaic_network.md)
