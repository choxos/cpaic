# Kaplan-Meier curves from the survival data behind a cML-NMR fit

Ported from
[`multinma::geom_km()`](https://dmphillippo.github.io/multinma/reference/geom_km.html)
(Phillippo et al. 2020). Returns a list of ggplot2 layers, so it can be
added to an existing plot; adding it to
[`plot_survival()`](https://choxos.github.io/cpaic/reference/plot_survival.md)
overlays the observed data on the fitted survival curves.

## Usage

``` r
geom_km(object, ..., curve_args = list(), cens_args = list())
```

## Arguments

- object:

  A [`cmlnmr()`](https://choxos.github.io/cpaic/reference/cmlnmr.md) fit
  with `family = "survival"`.

- ...:

  Passed to
  [`survival::survfit()`](https://rdrr.io/pkg/survival/man/survfit.html).

- curve_args, cens_args:

  Optional lists of arguments customizing the curves
  ([`ggplot2::geom_step()`](https://ggplot2.tidyverse.org/reference/geom_path.html))
  and the censoring marks
  ([`ggplot2::geom_point()`](https://ggplot2.tidyverse.org/reference/geom_point.html)).

## Value

A list of ggplot2 layers.

## Details

Only status `1` counts as an event; statuses `0`, `2`, and `3` (right,
left, and interval censoring) are treated as censored for the empirical
curve.

## See also

[`plot_survival()`](https://choxos.github.io/cpaic/reference/plot_survival.md)

## Examples

``` r
if (FALSE) {
plot_survival(fit) + geom_km(fit)
}
```
