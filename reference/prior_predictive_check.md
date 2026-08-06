# Summarize a prior-predictive cML-NMR fit

`cmlnmr(prior_predictive = TRUE)` samples from the prior without adding
the observed likelihood. This helper compares a simple statistic of the
observed outcomes with the corresponding replicated outcomes. Survival
replications are event-by-observed-time indicators because the censoring
process is not modeled. For Poisson models, the check stops if any
generated count exceeded Stan's safe random-number range. Affected
values are explicit sentinels and are never summarized as data.

## Usage

``` r
prior_predictive_check(object, statistic = c("mean", "sd"), level = 0.95)
```

## Arguments

- object:

  A [`cmlnmr()`](https://choxos.github.io/cpaic/reference/cmlnmr.md) fit
  created with `prior_predictive = TRUE`.

- statistic:

  Either `"mean"` or `"sd"`.

- level:

  Central prior-predictive interval level.

## Value

A data frame with observed and replicated summaries for IPD and AgD.
