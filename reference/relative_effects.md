# Relative treatment effects from a cpaic fit

Tidies the relative effects of the fitted model: every treatment versus
a chosen reference, or all pairwise comparisons. Effects are reported on
the natural scale of the summary measure (e.g. odds ratios) unless
`backtransf = FALSE`.

## Usage

``` r
relative_effects(
  object,
  reference = NULL,
  all_contrasts = FALSE,
  backtransf = TRUE,
  level = 0.95,
  newdata = NULL,
  estimand = NULL,
  ...
)
```

## Arguments

- object:

  A fitted cpaic object (`cpaic_bridge`, `cpaic_maic`, `cpaic_stc`, or
  `cpaic_mlnmr`).

- reference:

  Reference treatment. Defaults to the network reference.

- all_contrasts:

  If `TRUE`, return all pairwise comparisons instead of versus the
  reference.

- backtransf:

  If `TRUE` (default) back-transform log-scale measures (OR/RR/HR/...)
  by exponentiating.

- level:

  Confidence level for the intervals. Default `0.95`.

- newdata:

  For [`cmlnmr()`](https://choxos.github.io/cpaic/reference/cmlnmr.md)
  fits: a one-row data frame giving target effect-modifier means.
  Required when the model has effect modifiers.

- estimand:

  For [`cmlnmr()`](https://choxos.github.io/cpaic/reference/cmlnmr.md)
  fits, the only implemented value is `"average_conditional_link"`.
  Requests for `"marginal"` fail explicitly.

- ...:

  Unused.

## Value

A data frame with columns `treatment`, `comparator`, `estimate`,
`estimate_link`, `se_link`, `lower`, `upper`, `scale`, and `z`/`p` for
frequentist fits. `estimate_link` and `se_link` are always on the
model's link scale. `estimate`, `lower`, and `upper` use the scale named
in `scale`. For
[`cmlnmr()`](https://choxos.github.io/cpaic/reference/cmlnmr.md)
(Bayesian) fits the intervals are credible intervals and the final
column is `pr_gt0`, the posterior probability that the effect (on the
link scale) exceeds zero, instead of `z`/`p`.

## Details

Relative effects that the component design cannot uniquely identify
(their contrast vector lies outside the row space of `X = B C`) are
returned as `NA` rather than as pseudoinverse or prior-driven artefacts.
See
[`estimable_effects()`](https://choxos.github.io/cpaic/reference/estimable_effects.md).

For [`cmlnmr()`](https://choxos.github.io/cpaic/reference/cmlnmr.md)
fits the model contains component x effect-modifier interactions, so
relative effects depend on the supplied covariate means:
`theta_t(x) = C_t' (beta + gamma x)`. Because this expression is linear
in `x`, evaluating it at `E[X]` gives the average conditional link-scale
effect. It is not a marginal standardized OR, RR, or HR.

Two things this returns are narrower than they may look, both documented
in full under
[`cmlnmr()`](https://choxos.github.io/cpaic/reference/cmlnmr.md). The
value is the **average conditional link-scale** contrast at target means
`x`, not the marginal effect standardized over a population with a
distribution of covariates; on a non-collapsible scale (odds ratio,
hazard ratio) those differ, so `newdata = <a study's covariate means>`
does not give that study's population-average effect. And where the
interactions are informed only by aggregate arms, `x` is being applied
to an ecological gradient rather than to within-study effect
modification; see
[`estimable_effects_at()`](https://choxos.github.io/cpaic/reference/estimable_effects_at.md),
whose `identified_by` column separates the two.
