# Example disconnected component network: individual patient data

Individual patient data for studies `S3` (`A+B+C` vs `A+B`) and `S4`
(`A+B+D` vs `A+B`) of the
[cpaic_bin_agd](https://choxos.github.io/cpaic/reference/cpaic_bin_agd.md)
network. A single effect modifier `x1` is imbalanced relative to the
target population (`x1 = 0`), so population adjustment changes the `C`
and `D` component effects.

## Usage

``` r
cpaic_bin_ipd
```

## Format

A data frame with 3200 rows and 4 columns:

- .study:

  Study label (`S3` or `S4`).

- .trt:

  Treatment arm.

- .y:

  Binary outcome (0/1).

- x1:

  Continuous effect modifier.

The attribute `"truth"` holds the data-generating parameters.

## See also

[cpaic_bin_agd](https://choxos.github.io/cpaic/reference/cpaic_bin_agd.md),
[`cmaic()`](https://choxos.github.io/cpaic/reference/cmaic.md),
[`cstc()`](https://choxos.github.io/cpaic/reference/cstc.md)
