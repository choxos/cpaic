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
will not run on it. The underlying `cmdstanr` fit may still hold the
model data it was sampled with; for a fully data-free artifact, save
only the posterior draws (for example `fit$fit$draws()`).
