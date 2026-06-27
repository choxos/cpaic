# League table of all pairwise relative effects

League table of all pairwise relative effects

## Usage

``` r
league_table(object, backtransf = TRUE, level = 0.95, digits = 2)
```

## Arguments

- object:

  A `cpaic_bridge` / `cpaic_fit` object.

- backtransf, level:

  See
  [`relative_effects()`](https://choxos.github.io/cpaic/reference/relative_effects.md).

- digits:

  Rounding for the printed cells.

## Value

A character matrix (treatments x treatments); cell `[i, j]` is the
effect of the row treatment versus the column treatment with its
confidence interval.
