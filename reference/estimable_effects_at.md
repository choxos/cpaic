# Which population-adjusted contrasts are estimable at a target population?

Extends the row-space criterion of Wigle et al. (2026) from the
component main effects to the population-adjusted estimand
`theta_t(x) = C_t' (beta + Gamma x)`. A relative effect is identified by
the first-order information if and only if its augmented contrast vector
`(1, x) %x% (C_t - C_u)` lies in the row space of the information design
(see the file header for how that design is built from the IPD and
aggregate evidence).

## Usage

``` r
estimable_effects_at(object, newdata = NULL, reference = NULL, ...)
```

## Arguments

- object:

  A `cpaic_mlnmr` fit.

- newdata:

  A one-row data frame giving the target population's effect-modifier
  values. Defaults to the covariate origin.

- reference:

  Reference treatment. Defaults to the fit's reference.

- ...:

  Unused.

## Value

A data frame with `treatment`, `comparator`, `estimable`,
`identified_by` (`"IPD"`, `"aggregate"`, or `"none"`) and `basis`
(`"exact"`, `"first-order screen"`, or `"not identified"`); see the
section below.

## Details

Because the criterion depends on `x`, **the estimable set can depend on
the target population**: a contrast estimable at the covariate origin
need not be estimable in a target population where the component by
effect-modifier interactions are not identified.

## Strength of the guarantee

The `basis` column states how much the criterion actually proves for
each contrast, which is not the same for every row.

- `"exact"`:

  Either the contrast is identified by IPD, or the link is the identity.
  IPD identification is exact under **any** link: the IPD likelihood is
  an ordinary regression in arm and covariates, so the within-study
  arm-by-covariate variation pins down `m'beta` and `m'Gamma` directly.
  Under an identity link an aggregate arm's mean is linear in the
  parameters, so aggregate identification is exact too.

- `"first-order screen"`:

  The contrast is identified only through aggregate arms, under a
  nonlinear link (logit, log). The aggregate likelihood is then an
  integral over the covariate distribution, and a study pins the
  contrast down at a variance-weighted mean rather than at its raw
  covariate mean. The criterion has the right rank structure, so it
  finds under-determined contrasts correctly, but the anchor point is
  shifted and it can be **optimistic**. With a log link, one aggregate
  study and a symmetric covariate `P(x = -1) = P(x = +1) = 1/2`, the arm
  means are `exp(mu)` and `exp(mu + beta) cosh(gamma)`, so the data
  identify only `beta + log cosh(gamma)`, not `beta` itself. Treat such
  a contrast as reported under an additional smoothness assumption, and
  check it with
  [`prior_sensitivity()`](https://choxos.github.io/cpaic/reference/prior_sensitivity.md).

- `"not identified"`:

  Not in the row space of the first-order information. Any number
  reported here would be the prior, not the data. These contrasts are
  returned as `NA` by
  [`relative_effects()`](https://choxos.github.io/cpaic/reference/relative_effects.md)
  and dropped by
  [`cpaic_ranks()`](https://choxos.github.io/cpaic/reference/cpaic_ranks.md).

## References

Wigle A, Beliveau A, Nikolakopoulou A, Lin L (2026). Creating Treatment
and Component Hierarchies in Component Network Meta-Analysis.

## See also

[`prior_sensitivity()`](https://choxos.github.io/cpaic/reference/prior_sensitivity.md),
[`relative_effects()`](https://choxos.github.io/cpaic/reference/relative_effects.md),
[`cpaic_ranks()`](https://choxos.github.io/cpaic/reference/cpaic_ranks.md)
