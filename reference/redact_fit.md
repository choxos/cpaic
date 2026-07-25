# Strip raw individual patient data from a fitted cML-NMR object

A serialized
[`cmlnmr()`](https://choxos.github.io/cpaic/reference/cmlnmr.md) fit
retains the individual patient data in its refit arguments and
observed-outcome slots. `redact_fit()` removes those, so a saved or
shared object carries no row-level data; the posterior draws, component
design, diagnostics, and estimability information are preserved.

## Usage

``` r
redact_fit(object)
```

## Arguments

- object:

  A [`cmlnmr()`](https://choxos.github.io/cpaic/reference/cmlnmr.md)
  fit.

## Value

The fit with raw individual patient data removed, marked redacted.

## Details

After redaction the object can no longer be refitted, so
[`prior_sensitivity()`](https://choxos.github.io/cpaic/reference/prior_sensitivity.md)
will not run on it. Under either backend the sampler object in `fit$fit`
may still hold the model data it was sampled with, so for a fully
data-free artifact save only what you need from the posterior, such as
the output of
[`posterior_summary()`](https://choxos.github.io/cpaic/reference/posterior_summary.md),
rather than the fit itself.
