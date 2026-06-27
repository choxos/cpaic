# Reconnect a network through its additive component structure

Fits the additive component network meta-analysis (cNMA) model of Rücker
et al. (2020) to the aggregate contrast data, using
[`netmeta::discomb()`](https://rdrr.io/pkg/netmeta/man/discomb.html).
When the network is disconnected but its sub-networks share components,
the additive model estimates all component effects and so reconstructs
the (otherwise unavailable) relative effects *across* sub-networks. This
is the "connect first" step; population adjustment is layered on by
[`cmaic()`](https://choxos.github.io/cpaic/reference/cmaic.md) /
[`cstc()`](https://choxos.github.io/cpaic/reference/cstc.md), which
replace unadjusted contrasts with adjusted ones before calling this
function.

## Usage

``` r
cnma_bridge(network, common = FALSE, random = TRUE, ...)
```

## Arguments

- network:

  A
  [`cpaic_network()`](https://choxos.github.io/cpaic/reference/cpaic_network.md)
  object, or an aggregate contrast data frame (then `sm` and the column
  names are taken from `...`).

- common, random:

  Fit common- and/or random-effects models.

- ...:

  Additional arguments passed to
  [`netmeta::discomb()`](https://rdrr.io/pkg/netmeta/man/discomb.html)
  (e.g. `tau.preset`).

## Value

An object of class `cpaic_bridge` wrapping the
[`netmeta::discomb()`](https://rdrr.io/pkg/netmeta/man/discomb.html)
fit, with tidied component and treatment effects.

## References

Rücker G, Petropoulou M, Schwarzer G (2020). Network meta-analysis of
multicomponent interventions. *Biometrical Journal*, 62(3), 808–821.

## Examples

``` r
net <- cpaic_network(cpaic_bin_agd, sm = "OR", inactive = "Placebo")
br <- cnma_bridge(net)
component_effects(br)
#>   component  estimate        se     lower    upper statistic      pval
#> 1         A 0.5000000 1.1922140 -1.836697 2.836697 0.4193878 0.6749328
#> 2         B 0.4000000 1.1922140 -1.936697 2.736697 0.3355102 0.7372402
#> 3         C 0.7170248 0.9734562 -1.190914 2.624964 0.7365763 0.4613800
#> 4         D 0.3250136 0.9728622 -1.581761 2.231788 0.3340798 0.7383193
relative_effects(br)
#> Relative effects (OR, back-transformed)
#>  treatment comparator estimate    se lower   upper     z     p
#>          A    Placebo    1.649 1.192 0.159  17.059 0.419 0.675
#>        A+B    Placebo    2.460 1.686 0.090  66.993 0.534 0.593
#>      A+B+C    Placebo    5.038 1.947 0.111 228.801 0.831 0.406
#>      A+B+D    Placebo    3.404 1.947 0.075 154.510 0.629 0.529
#>          B    Placebo    1.492 1.192 0.144  15.436 0.336 0.737
```
