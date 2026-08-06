# Treatment and component hierarchies at target covariate means

Ranks treatments or components using average conditional link-scale
effects evaluated at supplied effect-modifier means. Because the
component contrast is linear in those means, the hierarchy can change
with them. This is not a hierarchy of marginal effects standardized over
a covariate distribution.

## Usage

``` r
cpaic_ranks(
  object,
  newdata = NULL,
  what = c("treatment", "component"),
  set = NULL,
  lower_is_better = FALSE,
  include_screen_only = FALSE,
  estimand = "average_conditional_link",
  ...
)
```

## Arguments

- object:

  A [`cmlnmr()`](https://choxos.github.io/cpaic/reference/cmlnmr.md)
  fit.

- newdata:

  A one-row data frame giving target effect-modifier means. Required
  when the model has effect modifiers.

- what:

  `"treatment"` (default) or `"component"`. Ranking components by their
  incremental effect is only meaningful in an additive model.

- set:

  Optional character vector restricting the elements to rank (the set
  `S` of Wigle et al.). Defaults to all treatments (or all components).

- lower_is_better:

  If `TRUE`, a smaller effect is preferred (e.g. mortality). Default
  `FALSE` (a larger effect is preferred).

- include_screen_only:

  If `FALSE` (default), elements whose relative effect is identified
  only by aggregate arms (a first-order screen that can be optimistic
  under a nonlinear link) are excluded from the hierarchy and reported
  in the `dropped_screen` attribute. Set `TRUE` to rank them as an
  explicitly exploratory hierarchy.

  The test here is whether the individual patient data identify the
  element, which is not identical to `basis == "exact"` in
  [`estimable_effects_at()`](https://choxos.github.io/cpaic/reference/estimable_effects_at.md).
  That column additionally excludes **survival** from `"exact"`, because
  a flexible baseline hazard adds support-dependent nuisance parameters
  the covariate-support argument does not see. Survival elements
  identified by IPD are therefore still ranked by default; dropping
  every survival element from every survival hierarchy would leave
  nothing to rank. Read a survival hierarchy alongside
  [`estimable_effects_at()`](https://choxos.github.io/cpaic/reference/estimable_effects_at.md)
  rather than on its own.

- estimand:

  The only implemented value is `"average_conditional_link"`. Marginal
  standardized rankings are not yet implemented and are rejected
  explicitly.

- ...:

  Unused.

## Value

A data frame, ordered from most to least preferred, with columns
`element`, `estimate` (posterior mean of the relative effect versus the
reference, on the link scale), `p_best`, `median_rank`, `mean_rank` and
`sucra`. The `dropped` attribute lists elements excluded as not
estimable at these target means.

## Details

Elements whose relative effect is not estimable at those target means
are **dropped from the ranking set** rather than ranked from a
prior-driven posterior, and are reported in the `dropped` attribute.
This is Step 3 of the Wigle et al. workflow, and it matters more here
than in the aggregate-data case, because the estimable set depends on
the target means (see
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
# Which component ranks best when the target mean of x1 is 0.5?
cpaic_ranks(fit, newdata = data.frame(x1 = 0.5), what = "component")
}
```
