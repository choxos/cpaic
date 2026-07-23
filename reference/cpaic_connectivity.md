# Assess connectivity and component-bridge identifiability of a network

Reports whether the treatment network is connected and, for a
disconnected network, which relative effects the additive component
structure makes estimable.

## Usage

``` r
cpaic_connectivity(network, tol = 1e-08)
```

## Arguments

- network:

  A
  [`cpaic_network()`](https://choxos.github.io/cpaic/reference/cpaic_network.md)
  object.

- tol:

  Numerical tolerance for the rank and null-space computations.

## Value

An object of class `cpaic_connectivity`: a list with `connected`
(logical), `n_subnetworks`, `subnetworks` (list of treatment-label
vectors), `bridging_components` (components shared across sub-networks),
`rank` and `n_components`, `identifiable` (logical:
`rank == n_components`), `null_space`, `estimable_components`,
`estimable` (a data frame of estimable relative effects versus the
reference), and the `B`/`C`/`X` matrices.

## Details

Two distinct questions are answered, and they are not the same (Wigle et
al. 2026):

- **Are all component effects identified?** Yes if and only if the
  component design matrix `X = B C` has full column rank (`rank(X)`
  equal to `n_components`). Reported as `identifiable`.

- **Is a particular relative effect estimable?** Yes if and only if its
  contrast vector lies in the row space of `X`. Full column rank is
  sufficient but **not necessary**, so a rank-deficient network can
  still identify useful cross-sub-network contrasts. Reported per
  treatment in `estimable` (see
  [`estimable_effects()`](https://choxos.github.io/cpaic/reference/estimable_effects.md)).

## References

Wigle A, Beliveau A, Nikolakopoulou A, Lin L (2026). Creating Treatment
and Component Hierarchies in Component Network Meta-Analysis.

## See also

[`estimable_effects()`](https://choxos.github.io/cpaic/reference/estimable_effects.md)

## Examples

``` r
net <- cpaic_network(cpaic_bin_agd, sm = "OR", inactive = "Placebo")
cpaic_connectivity(net)
#> cpaic connectivity
#>   Connected network: FALSE
#>   Sub-networks:      2
#>     [1] 3 treatments
#>     [2] 3 treatments
#>   Bridging components: A, B
#>     (components that OCCUR in more than one sub-network; occurrence is
#>      not identifiability, homogeneity, or influence for any contrast.)
#>   Component design:  rank(X) = 4 / 4 components -> all component effects identified
#>   Estimable effects: 5 / 5 vs Placebo
```
