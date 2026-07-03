# The cpaic statistical framework

This vignette describes the statistical framework behind cpaic. The
companion mathematical-foundations document (shipped with the
development sources) gives the full derivations; here we summarize the
model and show how each piece maps to a function.

## Two layers

cpaic targets networks that are *disconnected*: the treatments split
into two or more sub-networks with no common comparator, so standard
network meta-analysis cannot compare across the gap. Two layers solve
the two distinct problems.

1.  **Connection layer (component network meta-analysis).**
    Multi-component treatments are decomposed into additive component
    effects. When sub-networks share components, the component effects
    bridge the gap.
2.  **Adjustment layer (population adjustment).** Where individual
    patient data (IPD) are available, each evidence edge is corrected
    for between-study imbalance in effect modifiers, using anchored STC,
    MAIC, or ML-NMR.

The output is an indirect comparison that is both connected and
population-adjusted.

## The additive component model

Let $`\delta`$ be the vector of observed relative effects (one per
comparison), $`B`$ the edge-incidence (contrast) matrix mapping
comparisons to treatments, and $`C`$ the treatment-by-component matrix
with $`C_{tc} = 1`$ if treatment $`t`$ contains component $`c`$.
Treatment effects are additive in the component effects $`\beta`$,
``` math
 \theta = C\beta, \qquad \delta = B\theta = BC\beta = X\beta, 
```
with $`X = BC`$ the component design matrix. With inverse-variance
weights $`W`$, the component effects are estimated by weighted least
squares
``` math
 \hat\beta = (X^\top W X)^{+} X^\top W d, \qquad
   \mathrm{Cov}(\hat\theta) = C (X^\top W X)^{+} C^\top, 
```
where $`(\cdot)^{+}`$ is the Moore-Penrose inverse and $`d`$ the data
vector (Rücker et al. 2020). The additivity assumption is checked with a
Cochran $`Q`$ statistic. This is implemented in
[`cnma_bridge()`](https://choxos.github.io/cpaic/reference/cnma_bridge.md),
a wrapper around
[`netmeta::discomb()`](https://rdrr.io/pkg/netmeta/man/discomb.html).

## Connecting a disconnected network

A disconnected network can be bridged only if the shared components make
the component effects identifiable, that is $`\mathrm{rank}(X) = K`$,
the number of components.
[`cpaic_connectivity()`](https://choxos.github.io/cpaic/reference/cpaic_connectivity.md)
detects the sub-networks, lists the bridging components, and reports
identifiability.

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

When identifiable,
[`cnma_bridge()`](https://choxos.github.io/cpaic/reference/cnma_bridge.md)
reconstructs the relative effects across the gap from the component
effects.

``` r

component_effects(cnma_bridge(net))
#>   component  estimate        se     lower    upper statistic      pval
#> 1         A 0.5000000 1.1922140 -1.836697 2.836697 0.4193878 0.6749328
#> 2         B 0.4000000 1.1922140 -1.936697 2.736697 0.3355102 0.7372402
#> 3         C 0.7170248 0.9734562 -1.190914 2.624964 0.7365763 0.4613800
#> 4         D 0.3250136 0.9728622 -1.581761 2.231788 0.3340798 0.7383193
```

## Anchored simulated treatment comparison (cSTC)

For each IPD study,
[`cstc()`](https://choxos.github.io/cpaic/reference/cstc.md) fits an
outcome regression with treatment main effects, prognostic main effects,
and treatment-by-effect-modifier interactions, with the effect modifiers
centered at a target population. On the link scale,
``` math
 g\{E(y \mid \text{arm } t, x)\}
   = \mu + \beta_t + \gamma_t^\top (x - \bar x_{\text{target}}) + \dots 
```
so the treatment coefficient $`\beta_t`$ is the population-adjusted
contrast at the target (the interaction term vanishes at
$`x = \bar x_{\text{target}}`$). This is the anchored generalization of
regression-based standardization; it is implemented natively because
[`mlumr::stc()`](https://choxos.github.io/mlumr/reference/stc.html)
targets the unanchored two-trial case.

``` r

net_ipd <- cpaic_network(cpaic_bin_agd, ipd = cpaic_bin_ipd, sm = "OR",
                         family = "binomial", ipd_covariates = "x1",
                         inactive = "Placebo")
component_effects(cstc(net_ipd, target = c(x1 = 0), effect_modifiers = "x1"))
#>   component  estimate        se        lower     upper statistic        pval
#> 1         A 0.5000000 0.2563324 -0.002402322 1.0024023  1.950592 0.051105590
#> 2         B 0.4000000 0.2563324 -0.102402322 0.9024023  1.560474 0.118647988
#> 3         C 0.4896667 0.2406290  0.018042458 0.9612910  2.034944 0.041856471
#> 4         D 0.6408956 0.2317142  0.186744196 1.0950470  2.765889 0.005676788
```

## Anchored matching-adjusted indirect comparison (cMAIC)

[`cmaic()`](https://choxos.github.io/cpaic/reference/cmaic.md) reweights
each IPD study so that its effect-modifier distribution matches the
target population, using entropy-balancing weights
$`w_i = \exp(\tilde x_i^\top \alpha)`$ with $`\tilde x_i`$ the centered
effect modifiers (Phillippo et al. 2020). The effective sample size is
$`\mathrm{ESS} = (\sum_i w_i)^2 / \sum_i w_i^2`$. The weighted
within-study contrasts, with bootstrap standard errors that propagate
the weighting uncertainty, are then combined through the component
model.

``` r

fit_maic <- cmaic(net_ipd, target = c(x1 = 0), effect_modifiers = "x1",
                  n_boot = 100, seed = 1)
effective_sample_size(fit_maic)
#>       S3       S4 
#> 207.4202 358.1461
```

## The unification

For a single IPD edge anchored on a common comparator, combining the
adjusted contrast with the aggregate comparator contrast through the
component model is exactly a Bucher indirect comparison. The component
model generalizes this to a network: many adjusted and unadjusted
contrasts are combined simultaneously through
$`\hat\beta = (X^\top W X)^{+} X^\top W d`$, and the bridge supplies the
contrasts that are otherwise unavailable across the disconnection. The
key assumption is that component effects (and their effect-modifier
interactions) are constant across sub-networks (Rücker et al. 2021).

## Component-additive ML-NMR

[`cmlnmr()`](https://choxos.github.io/cpaic/reference/cmlnmr.md) places
the component structure inside multilevel network meta-regression: the
relative effect of an arm is $`C\beta`$ rather than a free per-treatment
parameter, and aggregate arms are fitted by integrating the individual
model over each study’s covariate distribution. Disconnected
sub-networks share the component parameters, so the network is connected
by construction. All four outcome families are supported; survival uses
a proportional-hazards model with a flexible baseline
(piecewise-exponential by default, or a smooth M-spline), set through
`cut_points` and `baseline`. Aggregate covariates are integrated with a
Gaussian copula whose correlation is estimated from the individual
patient data. A component whose effect-modifier interaction is informed
only by aggregate data is weakly identified and relies on the regression
prior.

``` r

# Requires cmdstanr; see ?cmlnmr.
fit <- cmlnmr(ipd, agd, effect_modifiers = "x1", inactive = "Placebo",
              family = "binomial")
component_effects(fit)
```

## Assumptions and caveats

- **Additivity** of component effects (test with
  [`additivity_test()`](https://choxos.github.io/cpaic/reference/additivity_test.md);
  add interaction terms if violated).
- **Identifiability**: a disconnected network is bridgeable only when
  [`cpaic_connectivity()`](https://choxos.github.io/cpaic/reference/cpaic_connectivity.md)
  reports `identifiable = TRUE`.
- **Cross-population transportability** of effect modifiers (the
  standard population-adjustment assumption), here extended to constancy
  of component effects across sub-networks.
- **Non-collapsibility**: MAIC targets a marginal effect and STC a
  conditional effect, so for noncollapsible measures they answer
  slightly different questions.

## References

Phillippo, David M., Sofia Dias, A. E. Ades, et al. 2020. “Multilevel
Network Meta-Regression for Population-Adjusted Treatment Comparisons.”
*Journal of the Royal Statistical Society: Series A* 183 (3): 1189–210.
<https://doi.org/10.1111/rssa.12579>.

Rücker, Gerta, Maria Petropoulou, and Guido Schwarzer. 2020. “Network
Meta-Analysis of Multicomponent Interventions.” *Biometrical Journal* 62
(3): 808–21. <https://doi.org/10.1002/bimj.201800167>.

Rücker, Gerta, Susanne Schmitz, and Guido Schwarzer. 2021. “Component
Network Meta-Analysis Compared to a Matching Method in a Disconnected
Network.” *Biometrical Journal* 63 (2): 447–61.
