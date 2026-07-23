# Posterior rank probabilities in a target population

The full rank distribution behind
[`cpaic_ranks()`](https://choxos.github.io/cpaic/reference/cpaic_ranks.md):
the posterior probability that each treatment (or component) takes each
rank, **in a named target population**. Ported from
[`multinma::posterior_rank_probs()`](https://dmphillippo.github.io/multinma/reference/posterior_ranks.html)
(Phillippo et al. 2020) and extended, because under population
adjustment the hierarchy is a function of the target: the component
effects are `beta + Gamma x`, so the ranks move with `x`.

## Usage

``` r
rank_probs(
  object,
  newdata = NULL,
  what = c("treatment", "component"),
  lower_is_better = FALSE,
  cumulative = FALSE,
  include_screen_only = FALSE,
  ...
)
```

## Arguments

- object:

  A [`cmlnmr()`](https://choxos.github.io/cpaic/reference/cmlnmr.md)
  fit.

- newdata:

  A one-row data frame giving the target population's effect-modifier
  values.

- what:

  `"treatment"` (default) or `"component"`.

- lower_is_better:

  If `TRUE`, a smaller effect is preferred.

- cumulative:

  Return cumulative rank probabilities (the quantity SUCRA summarizes)
  instead of the rankogram? Default `FALSE`.

- include_screen_only:

  If `FALSE` (default), elements identified only by aggregate arms (a
  first-order screen) are excluded, as in
  [`cpaic_ranks()`](https://choxos.github.io/cpaic/reference/cpaic_ranks.md).

- ...:

  Unused.

## Value

A data frame of class `cpaic_rank_probs` with one row per (element,
rank) and columns `element`, `rank_position`, and `probability`.

## Details

Elements that are not estimable at the target population are dropped
from the ranking set rather than ranked from the prior, exactly as in
[`cpaic_ranks()`](https://choxos.github.io/cpaic/reference/cpaic_ranks.md)
(Step 3 of Wigle et al. 2026); they are listed in the `dropped`
attribute.

## References

Wigle A, Beliveau A, Nikolakopoulou A, Lin L (2026). Creating Treatment
and Component Hierarchies in Component Network Meta-Analysis.

## See also

[`cpaic_ranks()`](https://choxos.github.io/cpaic/reference/cpaic_ranks.md),
[`rank_curve()`](https://choxos.github.io/cpaic/reference/rank_curve.md),
[`plot.cpaic_rank_probs()`](https://choxos.github.io/cpaic/reference/plot.cpaic_rank_probs.md)

## Examples

``` r
if (FALSE) {
rp <- rank_probs(fit, newdata = data.frame(x1 = 0.5), what = "component")
plot(rp)
}
```
