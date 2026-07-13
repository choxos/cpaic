# Which relative effects of a component network are uniquely estimable?

The additive component model identifies a relative effect
`theta_i - theta_j = (C_i - C_j)' beta` only when the contrast vector
`C_i - C_j` lies in the **row space** of the design matrix `X = B C`
(Wigle et al. 2026). Full column rank of `X` (rank equal to the number
of components) is *sufficient* for every contrast to be estimable, but
it is **not necessary**: a disconnected, rank-deficient component
network can still identify many cross-sub-network treatment contrasts.

## Usage

``` r
estimable_effects(object, reference = NULL, ...)
```

## Arguments

- object:

  A
  [`cpaic_network()`](https://choxos.github.io/cpaic/reference/cpaic_network.md),
  [`cpaic_connectivity()`](https://choxos.github.io/cpaic/reference/cpaic_connectivity.md),
  `cpaic_bridge` or `cpaic_mlnmr` object.

- reference:

  Reference treatment. Defaults to the network reference.

- ...:

  Unused.

## Value

A data frame with one row per treatment, giving the `treatment`, the
`comparator` (the reference), and `estimable` (logical).

## Details

Checking this matters because both engines otherwise return a
finite-looking answer for a contrast that carries no information: the
frequentist weighted least squares through the Moore-Penrose
pseudoinverse, and the Bayesian model through the prior.

## References

Wigle A, Beliveau A, Nikolakopoulou A, Lin L (2026). Creating Treatment
and Component Hierarchies in Component Network Meta-Analysis.

## Examples

``` r
net <- cpaic_network(cpaic_bin_agd, sm = "OR", inactive = "Placebo")
estimable_effects(net)
#>   treatment comparator estimable
#> 1         A    Placebo      TRUE
#> 2       A+B    Placebo      TRUE
#> 3     A+B+C    Placebo      TRUE
#> 4     A+B+D    Placebo      TRUE
#> 5         B    Placebo      TRUE
```
