# Forest plot of relative or component effects

A ggplot2 forest plot of the estimates in a cpaic fit. Ported from
[`multinma::plot.nma_summary()`](https://dmphillippo.github.io/multinma/reference/plot.nma_summary.html)
(Phillippo et al. 2020) and re-implemented on ggplot2 alone.

## Usage

``` r
forest(
  x,
  ...,
  what = c("relative", "component"),
  order = c("estimate", "alphabetical", "none"),
  ref_line = NULL,
  point_size = 1.2,
  show_na = TRUE
)

# S3 method for class 'cpaic_effects'
plot(x, y, ...)

# S3 method for class 'cpaic_bridge'
plot(x, y, ...)

# S3 method for class 'cpaic_fit'
plot(x, y, ...)
```

## Arguments

- x:

  A `cpaic_effects` data frame (from
  [`relative_effects()`](https://choxos.github.io/cpaic/reference/relative_effects.md)),
  a fitted cpaic object (`cpaic_bridge`, `cpaic_maic`, `cpaic_stc`,
  `cpaic_mlnmr`), or a component-effect data frame from
  [`component_effects()`](https://choxos.github.io/cpaic/reference/component_effects.md).

- ...:

  Passed to
  [`relative_effects()`](https://choxos.github.io/cpaic/reference/relative_effects.md)
  /
  [`component_effects()`](https://choxos.github.io/cpaic/reference/component_effects.md)
  when `x` is a fit (for example `newdata` for a
  [`cmlnmr()`](https://choxos.github.io/cpaic/reference/cmlnmr.md) fit).

- what:

  `"relative"` (default) for relative effects, or `"component"` for the
  incremental effect of each component.

- order:

  Row ordering: `"estimate"` (default, most to least favorable),
  `"alphabetical"`, or `"none"` (the order in the input).

- ref_line:

  Position of the vertical reference line. Defaults to the null value of
  the summary measure (`1` on a back-transformed ratio scale, `0`
  otherwise); `NA` draws none.

- point_size:

  Size of the point-estimate marker. Default `1.2`.

- show_na:

  Show non-estimable contrasts as labelled empty rows? Default `TRUE`.
  Setting this to `FALSE` hides evidence that the network cannot answer
  part of your question, so leave it on unless you have a reason.

- y:

  Unused, for compatibility with the
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) generic.

## Value

A `ggplot` object.

## Details

Contrasts that the component design cannot identify are **shown**,
labelled `not estimable`, rather than silently dropped. Dropping them
would leave the reader with a plot that looks complete when it is not;
see
[`estimable_effects()`](https://choxos.github.io/cpaic/reference/estimable_effects.md)
and
[`estimable_effects_at()`](https://choxos.github.io/cpaic/reference/estimable_effects_at.md).

A table with several comparators (from
`relative_effects(all_contrasts = TRUE)`) is faceted, one panel per
comparator. Ratio measures are drawn on a log axis, on which they are
symmetric. Pass `level` through `...` to change the interval width.

## See also

[`relative_effects()`](https://choxos.github.io/cpaic/reference/relative_effects.md),
[`component_effects()`](https://choxos.github.io/cpaic/reference/component_effects.md),
[`league_table()`](https://choxos.github.io/cpaic/reference/league_table.md)

## Examples

``` r
net <- cpaic_network(cpaic_bin_agd, sm = "OR", inactive = "Placebo")
br <- cnma_bridge(net)
forest(br)

forest(br, what = "component")
```
