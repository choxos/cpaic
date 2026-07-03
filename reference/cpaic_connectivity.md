# Assess connectivity and component-bridge identifiability of a network

Reports whether the treatment network is connected, and – crucially for
a disconnected network – whether the additive component structure makes
all component effects estimable. A disconnected network is *bridgeable*
only when sub-networks share enough components that the component design
matrix `X = B C` has full column rank
(`rank(X) = number of components`). Otherwise component effects cannot
be uniquely identified and the bridge would produce arbitrary estimates.

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

  Numerical tolerance for the rank computation.

## Value

An object of class `cpaic_connectivity`: a list with `connected`
(logical), `n_subnetworks`, `subnetworks` (list of treatment-label
vectors), `bridging_components` (components shared across sub-networks),
`rank` and `n_components`, `identifiable` (logical:
`rank == n_components`), and the `B`/`C`/`X` matrices.

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
#>   Component identifiability: rank(X) = 4 / 4 components -> IDENTIFIABLE
```
