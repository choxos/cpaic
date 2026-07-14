# Refit cML-NMR under tighter and looser priors

Prior movement is an empirical identification diagnostic. Contrasts that
move substantially when a weakly identified prior is changed should not
be interpreted as data-driven. This helper reuses the principle in
`documentation/validation/estimability_gamma.R`.

## Usage

``` r
prior_sensitivity(
  object,
  newdata,
  prior = c("gamma", "beta", "all"),
  tighter = 0.5,
  looser = 2,
  reference = NULL,
  ...
)
```

## Arguments

- object:

  A [`cmlnmr()`](https://choxos.github.io/cpaic/reference/cmlnmr.md)
  fit.

- newdata:

  One target population, as for
  [`relative_effects()`](https://choxos.github.io/cpaic/reference/relative_effects.md).

- prior:

  Which scales to vary: the interaction prior, component-effect prior,
  or all configurable scale priors.

- tighter, looser:

  Positive multipliers for the fitted prior scales.

- reference:

  Reference treatment. Defaults to the fit reference.

- ...:

  Named arguments overriding the stored refit call, such as fewer
  sampling iterations for a screening run.

## Value

A `cpaic_prior_sensitivity` object containing the movement table and the
tighter and looser fits.
