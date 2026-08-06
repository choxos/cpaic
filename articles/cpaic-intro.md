# Getting started with cpaic

``` r

library(cpaic)
```

> **Research use only.** cpaic is experimental. Its methodology and
> implementation have not been validated for clinical, regulatory,
> reimbursement, or other decision use.

## The problem

Network meta-analysis needs a *connected* network. When the evidence
splits into two or more sub-networks with no common comparator, the
network is **disconnected** and the treatments in different sub-networks
cannot be compared directly.

`cpaic` can reconnect such a network through the *additive component*
structure of the treatments. It also implements experimental adjustment
for between-study differences in effect modifiers using anchored STC,
MAIC, and ML-NMR. Algebraic reconnection does not by itself make the
result transportable to one population.

## A disconnected example

The bundled data describe a binary-outcome network in two pieces:

- sub-network 1, anchored on placebo: `Placebo`, `A`, `B`;
- sub-network 2, isolated: `A+B`, `A+B+C`, `A+B+D`.

No treatment is shared between the two pieces, so the network is
disconnected. The shared components `A` and `B` bridge it.

``` r

net <- cpaic_network(cpaic_bin_agd, ipd = cpaic_bin_ipd, sm = "OR",
                     family = "binomial", ipd_covariates = "x1",
                     inactive = "Placebo")
net
#> cpaic component network
#>   Summary measure:   OR
#>   Treatments:        6
#>   Components:        4 (A, B, C, D)
#>   AgD comparisons:   5
#>   Reference:         Placebo
#>   Inactive:          Placebo
#>   IPD studies:       2 (binomial; 3200 patients)
#>   Connected:         FALSE | components bridgeable: TRUE
```

## Is the network bridgeable?

A disconnected network can be bridged only if the shared components make
all component effects identifiable, that is
`rank(X) = number of components`.

``` r

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

The report confirms two sub-networks, identifies `A` and `B` as the
bridging components, and shows the component effects are identifiable.

``` r

plot(net)
```

![Component network coloured by
sub-network](cpaic-intro_files/figure-html/unnamed-chunk-4-1.png)

## Step 1: connect with component NMA

[`cnma_bridge()`](https://choxos.github.io/cpaic/reference/cnma_bridge.md)
fits the additive component model and reconstructs the relative effects
across the gap.

``` r

br <- cnma_bridge(net)
component_effects(br)
#>   component  estimate        se     lower    upper statistic      pval
#> 1         A 0.5000000 1.1922140 -1.836697 2.836697 0.4193878 0.6749328
#> 2         B 0.4000000 1.1922140 -1.936697 2.736697 0.3355102 0.7372402
#> 3         C 0.7170248 0.9734562 -1.190914 2.624964 0.7365763 0.4613800
#> 4         D 0.3250136 0.9728622 -1.581761 2.231788 0.3340798 0.7383193
```

## Step 2: adjust for effect modifiers

Components `C` and `D` come from the IPD studies, whose effect modifier
`x1` is imbalanced relative to the target mean (`x1 = 0`). Anchored STC
fits an outcome regression with treatment-by-`x1` interactions and
evaluates the conditional link-scale treatment effect at that target
mean.

``` r

fit_stc <- cstc(net, target = c(x1 = 0), effect_modifiers = "x1",
                allow_experimental_bridge = TRUE)
#> Warning: cstc() cannot form a decision-grade component bridge:
#>   - retained aggregate-only edge(s) remain in their own study populations: S1: A vs Placebo; S2: B vs Placebo; S5: A+B+C vs A+B+D
#> Use cmlnmr() for a joint model, restrict the analysis to a design in which every edge is adjusted and the estimand is additive, or set `allow_experimental_bridge = TRUE` only for explicitly exploratory sensitivity work.
component_effects(fit_stc)
#>   component  estimate        se        lower     upper statistic        pval
#> 1         A 0.5000000 0.2563324 -0.002402322 1.0024023  1.950592 0.051105590
#> 2         B 0.4000000 0.2563324 -0.102402322 0.9024023  1.560474 0.118647988
#> 3         C 0.4896667 0.2406290  0.018042458 0.9612910  2.034944 0.041856471
#> 4         D 0.6408956 0.2317142  0.186744196 1.0950470  2.765889 0.005676788
```

Anchored MAIC instead reweights each IPD study to the target moments.
Its marginal weighted contrast is not generally component-additive on a
nonlinear scale such as the log odds scale.

``` r

fit_maic <- cmaic(net, target = c(x1 = 0), effect_modifiers = "x1",
                  n_boot = 100, seed = 1,
                  allow_experimental_bridge = TRUE)
#> Warning: cmaic() cannot form a decision-grade component bridge:
#>   - retained aggregate-only edge(s) remain in their own study populations: S1: A vs Placebo; S2: B vs Placebo; S5: A+B+C vs A+B+D
#>   - cMAIC estimates marginal binomial contrasts, which are not generally additive in the component design on a nonlinear link scale
#> Use cmlnmr() for a joint model, restrict the analysis to a design in which every edge is adjusted and the estimand is additive, or set `allow_experimental_bridge = TRUE` only for explicitly exploratory sensitivity work.
effective_sample_size(fit_maic)
#>       S3       S4 
#> 207.4202 358.1461
```

Both calls require explicit opt-in because this example retains
unadjusted AgD edges from their original study populations. The fitted
objects record the override. It is not a claim that the resulting
network has one coherent target estimand.

Adjusting the IPD edges moves the `C` and `D` estimates relative to the
unadjusted bridge, while the estimates driven by retained aggregate
edges are unchanged. This comparison is descriptive because the opt-in
bridge mixes study-population and target-specific contrasts:

``` r

data.frame(
  component = component_effects(br)$component,
  naive     = round(component_effects(br)$estimate, 3),
  cSTC      = round(component_effects(fit_stc)$estimate, 3),
  cMAIC     = round(component_effects(fit_maic)$estimate, 3)
)
#>   component naive  cSTC cMAIC
#> 1         A 0.500 0.500 0.500
#> 2         B 0.400 0.400 0.400
#> 3         C 0.717 0.490 0.697
#> 4         D 0.325 0.641 0.772
```

## Reporting

``` r

relative_effects(fit_stc)
#> Relative effects (OR, back-transformed)
#>  treatment comparator estimate estimate_link se_link lower  upper   scale     z
#>          A    Placebo    1.649         0.500   0.256 0.998  2.725 natural 1.951
#>        A+B    Placebo    2.460         0.900   0.363 1.209  5.005 natural 2.483
#>      A+B+C    Placebo    4.014         1.390   0.435 1.711  9.416 natural 3.194
#>      A+B+D    Placebo    4.669         1.541   0.430 2.009 10.850 natural 3.582
#>          B    Placebo    1.492         0.400   0.256 0.903  2.466 natural 1.560
#>      p
#>  0.051
#>  0.013
#>  0.001
#>  0.000
#>  0.119
#>   `se_link` is on the link (log) scale; the interval is back-transformed.
additivity_test(fit_stc)
#> Additive component model: fit statistics
#>   Total lack of fit (Q.additive): Q = 2.669, df = 1, p = 0.102
#>   Additivity restrictions (Q.diff): not available; no standard NMA
#>     is estimable on a disconnected network.
#>   Note: neither statistic tests whether component effects are constant
#>   ACROSS sub-networks, which is the assumption that bridges the gap.
#>   That assumption is untestable from the data and must be justified
#>   clinically.
```

``` r

forest(fit_stc)
```

![Forest plot of relative effects versus
placebo](cpaic-intro_files/figure-html/unnamed-chunk-10-1.png)

## Where next

- [`vignette("cpaic-methods")`](https://choxos.github.io/cpaic/articles/cpaic-methods.md)
  covers the statistical framework in depth.
- The full mathematical foundations and a validation study are provided
  with the development sources.
