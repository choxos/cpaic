# Effective sample sizes from a cMAIC fit

The effective sample size summarizes the precision lost to reweighting.
It is **not** a validity diagnostic: a healthy ESS says the weights are
well behaved, not that the reweighted edge is relevant to your estimand.
Use
[`edge_influence()`](https://choxos.github.io/cpaic/reference/edge_influence.md)
to ask whether the IPD informs the contrast at all.

## Usage

``` r
effective_sample_size(object, ...)
```

## Arguments

- object:

  A `cpaic_maic` object.

- ...:

  Unused.

## Value

A named numeric vector of effective sample sizes per IPD study.

## See also

[`edge_influence()`](https://choxos.github.io/cpaic/reference/edge_influence.md)
