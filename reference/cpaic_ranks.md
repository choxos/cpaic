# Population-adjusted treatment and component hierarchies

Ranks treatments or components *in a named target population*, following
the workflow of Wigle et al. (2026) but with every quantity evaluated at
the target's effect-modifier values. Because the component effects are
population-specific under population adjustment, so is the hierarchy: a
component may rank first in one population and last in another.

## Usage

``` r
cpaic_ranks(
  object,
  newdata = NULL,
  what = c("treatment", "component"),
  set = NULL,
  lower_is_better = FALSE,
  ...
)
```

## Arguments

- object:

  A [`cmlnmr()`](https://choxos.github.io/cpaic/reference/cmlnmr.md)
  fit.

- newdata:

  A one-row data frame giving the target population's effect-modifier
  values. Required when the model has effect modifiers.

- what:

  `"treatment"` (default) or `"component"`. Ranking components by their
  incremental effect is only meaningful in an additive model.

- set:

  Optional character vector restricting the elements to rank (the set
  `S` of Wigle et al.). Defaults to all treatments (or all components).

- lower_is_better:

  If `TRUE`, a smaller effect is preferred (e.g. mortality). Default
  `FALSE` (a larger effect is preferred).

- ...:

  Unused.

## Value

A data frame, ordered from most to least preferred, with columns
`element`, `estimate` (posterior mean of the relative effect versus the
reference, on the link scale), `p_best`, `median_rank`, `mean_rank` and
`sucra`. The `dropped` attribute lists elements excluded as not
estimable in this target population.

## Details

Elements whose relative effect is not estimable at that target
population are **dropped from the ranking set** rather than ranked from
a prior-driven posterior, and are reported in the `dropped` attribute.
This is Step 3 of the Wigle et al. workflow, and it matters more here
than in the aggregate-data case, because the estimable set depends on
the target (see
[`estimable_effects_at()`](https://choxos.github.io/cpaic/reference/estimable_effects_at.md)).

Ranking metrics depend on the set being ranked, so they are not
comparable across different sets. Report them alongside the relative
effects, never instead of them.

## References

Wigle A, Beliveau A, Nikolakopoulou A, Lin L (2026). Creating Treatment
and Component Hierarchies in Component Network Meta-Analysis.

## See also

[`estimable_effects_at()`](https://choxos.github.io/cpaic/reference/estimable_effects_at.md),
[`relative_effects()`](https://choxos.github.io/cpaic/reference/relative_effects.md)

## Examples

``` r
if (FALSE) {
# Which component is best for a patient population with x1 = 0.5?
cpaic_ranks(fit, newdata = data.frame(x1 = 0.5), what = "component")
}
```
