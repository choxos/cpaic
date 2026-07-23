# Bridge fragility: how much cross-sub-network drift would change a conclusion

In a disconnected network the cross-gap contrast exists only because the
component effects are assumed constant *across* sub-networks. That
assumption cannot be tested from the data, because there is no cross-gap
evidence. `bridge_fragility()` quantifies how sensitive a requested
contrast is to a violation of it.

## Usage

``` r
bridge_fragility(
  object,
  treatment,
  comparator = NULL,
  newdata = NULL,
  threshold = 0,
  plausible_drift = NULL,
  ...
)
```

## Arguments

- object:

  A [`cmlnmr()`](https://choxos.github.io/cpaic/reference/cmlnmr.md)
  fit.

- treatment, comparator:

  The contrast to assess. `comparator` defaults to the fit reference.

- newdata:

  A one-row data frame giving the target population's effect-modifier
  values (required when the model has effect modifiers).

- threshold:

  Decision boundary on the link scale. Default `0` (no effect).

- plausible_drift:

  Optional per-component drift bound (link scale) at which to report the
  posterior probability that the conclusion is robust.

- ...:

  Unused.

## Value

An object of class `cpaic_fragility`: the contrast, the L1 drift
loading, the posterior of the bridge fragility threshold, and (if
`plausible_drift` is given) the probability the conclusion survives it.

## Details

On the linear-predictor scale the contrast is \\D = m'(\beta + \Gamma
x)\\ with \\m = C_t - C_u\\. A cross-sub-network drift \\\Delta\\ in the
component effects shifts it to \\D + m'\Delta\\. Bounding each
component's drift by \\\|\Delta_c\| \le d\\, the worst-case shift is \\d
\sum_c \|m_c\|\\, so the smallest per-component drift that moves the
contrast to a decision threshold \\\tau\\ (default 0, on the link scale)
is the **bridge fragility threshold** \$\$\mathrm{BFT} = \|D - \tau\| /
\textstyle\sum_c \|m_c\|,\$\$ reported per posterior draw. A small BFT
means a clinically trivial amount of un-testable drift would overturn
the conclusion. This is a conservative worst-case over the component
*main-effect* drift; interaction drift \\\Lambda\\ is not included, so
the true fragility is no larger than reported.

## See also

[`relative_effects()`](https://choxos.github.io/cpaic/reference/relative_effects.md),
[`estimable_effects_at()`](https://choxos.github.io/cpaic/reference/estimable_effects_at.md)

## Examples

``` r
if (FALSE) {
bridge_fragility(fit, treatment = "A+B", newdata = data.frame(x1 = 0))
}
```
