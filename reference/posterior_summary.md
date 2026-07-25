# Posterior summary of a component ML-NMR fit

Summarizes the posterior draws of a
[`cmlnmr()`](https://choxos.github.io/cpaic/reference/cmlnmr.md) fit:
the component effects `beta`, the component by effect-modifier
interactions `gamma`, the study baselines `mu`, the heterogeneity `tau`,
and any other sampled block, with the usual convergence quantities
alongside.

Use this rather than reaching into `fit$fit`. That object is whatever
the sampler backend returned: an S4 `stanfit` under `backend = "rstan"`
and an R6 object under `backend = "cmdstanr"`. The two share no
accessors, so `fit$fit$summary(...)` works on one and fails on the
other. This function works on both and returns the same columns either
way.

## Usage

``` r
posterior_summary(x, variables = NULL, ...)
```

## Arguments

- x:

  A `cpaic_mlnmr` fit from
  [`cmlnmr()`](https://choxos.github.io/cpaic/reference/cmlnmr.md).

- variables:

  Character vector of Stan variable names to summarize, for example
  `"tau"` or `c("beta", "gamma")`. Naming a block returns one row per
  element of it. The default summarizes every sampled block the fit has.

- ...:

  Further summary functions passed to
  [`posterior::summarise_draws()`](https://mc-stan.org/posterior/reference/draws_summary.html),
  for example `"quantile2"` or a function. With none given, the default
  set is returned.

## Value

A data frame with one row per scalar parameter. With the default
summaries the columns are `variable`, `mean`, `median`, `sd`, `mad`,
`q5`, `q95`, `rhat`, `ess_bulk`, and `ess_tail`. Because both backends
are summarized through the same code, `rhat`, `ess_bulk`, and `ess_tail`
are the same quantities whichever engine produced the fit.

## See also

[`cmlnmr()`](https://choxos.github.io/cpaic/reference/cmlnmr.md) for the
fit,
[`relative_effects()`](https://choxos.github.io/cpaic/reference/relative_effects.md)
and
[`component_effects()`](https://choxos.github.io/cpaic/reference/component_effects.md)
for effects on the outcome scale rather than the parameters themselves,
and
[`redact_fit()`](https://choxos.github.io/cpaic/reference/redact_fit.md),
which strips the draws.

## Examples

``` r
if (FALSE) { # \dontrun{
fit <- cmlnmr(ipd, agd, effect_modifiers = "x1", inactive = "Placebo")
posterior_summary(fit, "tau")
min(posterior_summary(fit)$ess_bulk)
} # }
```
