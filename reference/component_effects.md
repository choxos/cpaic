# Component effects from a cpaic fit

Component effects from a cpaic fit

## Usage

``` r
component_effects(object, newdata = NULL, ...)

# S3 method for class 'cpaic_mlnmr'
component_effects(object, newdata = NULL, level = 0.95, ...)
```

## Arguments

- object:

  A fitted cpaic object (`cpaic_bridge`, `cpaic_maic`, `cpaic_stc`, or
  `cpaic_mlnmr`).

- newdata:

  For [`cmlnmr()`](https://choxos.github.io/cpaic/reference/cmlnmr.md)
  fits: a one-row data frame giving the target population's
  effect-modifier values, at which the component effects
  `beta + Gamma x` are reported. With `newdata = NULL` (default) the
  component *main* effects `beta` are returned; these are the effects at
  the covariate origin and are not population-adjusted.

- ...:

  Passed to methods.

- level:

  Credible level for the component-effect intervals (default `0.95`),
  for [`cmlnmr()`](https://choxos.github.io/cpaic/reference/cmlnmr.md)
  fits.

## Value

A data frame of component effects (estimate, se, CI, p-value).
Components that the design cannot identify are returned as `NA`.
