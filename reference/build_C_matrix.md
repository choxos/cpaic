# Build a component-coded treatment-by-component matrix

Splits multi-component treatment labels (e.g. `"A + B"`) on `sep.comps`
and returns the binary treatment-by-component design matrix `C`, where
`C[t, c] = 1` if treatment `t` contains component `c`. The `inactive`
component (e.g. placebo) is represented by an all-zero row, matching the
convention of
[`netmeta::netcomb()`](https://rdrr.io/pkg/netmeta/man/netcomb.html).

## Usage

``` r
build_C_matrix(treatments, sep.comps = "+", inactive = NULL)
```

## Arguments

- treatments:

  Character vector of (unique) treatment labels.

- sep.comps:

  Single character separating components in a treatment label. Default
  `"+"`.

- inactive:

  Optional name of the inactive/reference treatment or component (mapped
  to a zero row / dropped as a column).

## Value

A binary matrix with one row per treatment and one column per component
(treatments as row names, components as column names).

## Examples

``` r
build_C_matrix(c("A", "B", "A + B", "placebo"), inactive = "placebo")
#>         A B
#> A       1 0
#> B       0 1
#> A + B   1 1
#> placebo 0 0
```
