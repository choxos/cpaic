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

A data frame with `treatment`, `comparator`, `estimable`, and
`identified_by` (`"IPD"`, `"aggregate"`, or `"none"`).

## Details

Because the criterion depends on `x`, **the estimable set can depend on
the target population**: a contrast estimable at the covariate origin
need not be estimable in a target population where the component by
effect-modifier interactions are not identified.

## References

Wigle A, Beliveau A, Nikolakopoulou A, Lin L (2026). Creating Treatment
and Component Hierarchies in Component Network Meta-Analysis.
