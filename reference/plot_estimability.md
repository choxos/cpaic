# Map which contrasts are estimable, and on what evidence, across populations

Estimability under population adjustment is itself a function of the
target population: a contrast identified at the covariate origin need
not be identified in a population where the relevant component by
effect-modifier interactions are not pinned down. This plot evaluates
[`estimable_effects_at()`](https://choxos.github.io/cpaic/reference/estimable_effects_at.md)
over a grid of target populations and tiles the result, separating
contrasts identified by **IPD** (a within-study interaction, which
randomization protects) from those identified only **ecologically**,
from between-study differences in aggregate covariate means (which
randomization does not protect; Berlin et al. 2002).

## Usage

``` r
plot_estimability(object, em, values, at = NULL, reference = NULL, ...)
```

## Arguments

- object:

  A [`cmlnmr()`](https://choxos.github.io/cpaic/reference/cmlnmr.md)
  fit.

- em:

  Name of the effect modifier to vary across the grid.

- values:

  Numeric vector of target values for `em`.

- at:

  Optional named vector fixing the other effect modifiers. Defaults to 0
  for each.

- reference:

  Reference treatment. Defaults to the fit's reference.

- ...:

  Unused.

## Value

A `ggplot` object.

## Details

There is no counterpart in multinma.

## References

Berlin JA, Santanna J, Schmid CH, Szczech LA, Feldman HI (2002).
Individual patient- versus group-level data meta-regressions for the
investigation of treatment effect modifiers. *Statistics in Medicine*,
21(3), 371–387.

## See also

[`estimable_effects_at()`](https://choxos.github.io/cpaic/reference/estimable_effects_at.md),
[`plot_rank_curve()`](https://choxos.github.io/cpaic/reference/plot_rank_curve.md)

## Examples

``` r
if (FALSE) {
plot_estimability(fit, em = "x1", values = seq(-1, 1, by = 0.5))
}
```
