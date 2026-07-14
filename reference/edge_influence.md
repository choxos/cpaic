# Does the individual patient data actually inform this contrast?

Population adjustment only helps if the adjusted edges actually *carry*
the contrast you care about. They need not. In a component bridge the
estimate of a contrast `m' beta` is a weighted combination of the
observed edges,

## Usage

``` r
edge_influence(object, treatment, comparator = NULL, tol = 1e-08, ...)
```

## Arguments

- object:

  A `cpaic_bridge`, `cpaic_maic` or `cpaic_stc` object.

- treatment, comparator:

  The contrast of interest. `comparator` defaults to the network
  reference.

- tol:

  Influence weights below this (relative to the largest) are treated as
  zero.

- ...:

  Unused.

## Value

A data frame with one row per edge: `studlab`, `treat1`, `treat2`,
`has_ipd`, and `influence` (the weight `w_j`). Edges are ordered by
absolute influence. A warning is issued if any IPD edge has no influence
on the requested contrast.

## Details

\$\$m'\hat\beta = \underbrace{m' (X'WX)^{+} X'W}\_{w} \\ d ,\$\$

so edge `j` influences the answer only through its weight `w_j`. An IPD
edge with `w_j` of zero contributes nothing to that contrast, and
adjusting it changes nothing.

This matters because the usual diagnostic cannot detect the problem. In
simulation, putting the IPD on an edge that does not bridge the gap left
cMAIC numerically identical to the unadjusted analysis (bias +0.374,
coverage 0.676) while
[`effective_sample_size()`](https://choxos.github.io/cpaic/reference/effective_sample_size.md)
happily reported an ESS of 999 out of 1000. A healthy ESS says the
*weights* are well behaved; it says nothing about whether the reweighted
edge is relevant to your estimand.

## See also

[`effective_sample_size()`](https://choxos.github.io/cpaic/reference/effective_sample_size.md),
[`estimable_effects()`](https://choxos.github.io/cpaic/reference/estimable_effects.md)

## Examples

``` r
net <- cpaic_network(cpaic_bin_agd, sm = "OR", inactive = "Placebo")
br <- cnma_bridge(net)
edge_influence(br, treatment = "A+B+C")
#>   studlab treat1  treat2 has_ipd influence
#> 1      S1      A Placebo   FALSE 1.0000000
#> 2      S2      B Placebo   FALSE 1.0000000
#> 3      S3  A+B+C     A+B   FALSE 0.7065349
#> 4      S4  A+B+D     A+B   FALSE 0.2934651
#> 5      S5  A+B+C   A+B+D   FALSE 0.2934651
```
