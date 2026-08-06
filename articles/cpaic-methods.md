# The cpaic statistical framework

> **Research use only.** cpaic is experimental. Its methodology and
> implementation have not been validated for clinical, regulatory,
> reimbursement, or other decision use.

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

The component model can identify a cross-gap contrast, but scientific
transportability is an additional assumption. A two-stage result is not
a coherent single-target population adjustment when target-specific IPD
contrasts are mixed with retained AgD contrasts from their original
study populations.

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

A full set of component effects is identifiable when
$`\mathrm{rank}(X) = K`$, the number of components. Even in a
rank-deficient design, a particular contrast can be identifiable when
its contrast vector lies in the row space of $`X`$.
[`cpaic_connectivity()`](https://choxos.github.io/cpaic/reference/cpaic_connectivity.md)
detects the sub-networks and reports the component-design diagnostics;
[`estimable_effects()`](https://choxos.github.io/cpaic/reference/estimable_effects.md)
checks the requested contrasts.

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
centered at supplied target covariate means. On the link scale,
``` math
 g\{E(y \mid \text{arm } t, x)\}
   = \mu + \beta_t + \gamma_t^\top (x - \bar x_{\text{target}}) + \dots 
```
so the treatment coefficient $`\beta_t`$ is the conditional link-scale
contrast at the supplied target covariate means (the interaction term
vanishes at $`x = \bar x_{\text{target}}`$). Because this contrast is
linear in the effect modifiers on the link scale, it is also the average
conditional link-scale contrast at those means. It is not an
outcome-scale marginal standardization. The model is implemented
natively.

``` r

net_ipd <- cpaic_network(cpaic_bin_agd, ipd = cpaic_bin_ipd, sm = "OR",
                         family = "binomial", ipd_covariates = "x1",
                         inactive = "Placebo")
component_effects(cstc(net_ipd, target = c(x1 = 0),
                       effect_modifiers = "x1",
                       allow_experimental_bridge = TRUE))
#> Warning: cstc() cannot form a decision-grade component bridge:
#>   - retained aggregate-only edge(s) remain in their own study populations: S1: A vs Placebo; S2: B vs Placebo; S5: A+B+C vs A+B+D
#> Use cmlnmr() for a joint model, restrict the analysis to a design in which every edge is adjusted and the estimand is additive, or set `allow_experimental_bridge = TRUE` only for explicitly exploratory sensitivity work.
#>   component  estimate        se        lower     upper statistic        pval
#> 1         A 0.5000000 0.2563324 -0.002402322 1.0024023  1.950592 0.051105590
#> 2         B 0.4000000 0.2563324 -0.102402322 0.9024023  1.560474 0.118647988
#> 3         C 0.4896667 0.2406290  0.018042458 0.9612910  2.034944 0.041856471
#> 4         D 0.6408956 0.2317142  0.186744196 1.0950470  2.765889 0.005676788
```

## Anchored matching-adjusted indirect comparison (cMAIC)

[`cmaic()`](https://choxos.github.io/cpaic/reference/cmaic.md) reweights
each IPD study so that its effect-modifier moments match the supplied
target moments, using entropy-balancing weights
$`w_i = \exp(\tilde x_i^\top \alpha)`$ with $`\tilde x_i`$ the centered
effect modifiers (Phillippo et al. 2020). The effective sample size is
$`\mathrm{ESS} = (\sum_i w_i)^2 / \sum_i w_i^2`$. The weighted
within-study contrasts, with bootstrap standard errors that propagate
the weighting uncertainty, are then combined through the component
model.

That last step is an additional methodology constraint. A marginal
weighted contrast is not generally additive across treatment components
on a nonlinear link scale. Also, if AgD-only edges remain, those
contrasts still refer to their original study populations. cMAIC
therefore stops by default outside the Gaussian identity-link case and
whenever retained AgD edges would be mixed in. The experimental override
does not repair either estimand mismatch.

``` r

fit_maic <- cmaic(net_ipd, target = c(x1 = 0), effect_modifiers = "x1",
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

## The unification

For a single anchored comparison, the contrast algebra resembles a
Bucher indirect comparison. A network bridge combines many adjusted and
unadjusted contrasts through
$`\hat\beta = (X^\top W X)^{+} X^\top W d`$. That algebra does not
establish that every input contrast shares the same target estimand. It
also does not test the required constancy of component effects and
effect-modifier interactions across disconnected sub-networks (Rücker et
al. 2021).

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
`cut_points` and `baseline`. Its analytic row-level likelihood handles
events, right/left/interval censoring, and delayed entry. Aggregate
survival arms require reconstructed event and censoring rows; event
counts plus person-time are rejected. The finite quasi-Monte Carlo
covariate integration is approximate. Restricted mean survival time is
not implemented. Aggregate covariates are integrated with a Gaussian
copula whose correlation is estimated from the individual patient data.
A component whose effect-modifier interaction is informed only by
aggregate data is weakly identified and relies on the regression prior.

Reporting methods evaluate $`(C_t-C_u)'(\beta + \Gamma x)`$ at a one-row
vector of supplied covariate means. Since this expression is linear in
$`x`$, it is an average conditional link-scale effect at those means. It
is not a marginally standardized odds, rate, or hazard ratio. The
aggregate likelihood integrates outcomes over covariate distributions
for model fitting, but that integration does not turn the reporting
function into marginal standardization.

``` r

# Fits with rstan by default; backend = "cmdstanr" is the alternative.
fit <- cmlnmr(ipd, agd, effect_modifiers = "x1", inactive = "Placebo",
              family = "binomial")
component_effects(fit)
```

## Assumptions and caveats

- **Additivity** of component effects.
  [`additivity_test()`](https://choxos.github.io/cpaic/reference/additivity_test.md)
  can diagnose lack of fit where evidence exists, but it cannot validate
  an unobserved cross-gap bridge.
- **Identifiability**: a disconnected network is bridgeable only when
  [`cpaic_connectivity()`](https://choxos.github.io/cpaic/reference/cpaic_connectivity.md)
  reports `identifiable = TRUE`.
- **Cross-population transportability** of effect modifiers (the
  standard population-adjustment assumption), here extended to constancy
  of component effects across sub-networks.
- **Two-stage estimand mixing**: retained AgD contrasts remain tied to
  their study populations. The default gate stops this mixed bridge
  unless the caller explicitly opts into research-only behavior.
- **Nonlinearity and non-collapsibility**: MAIC targets a marginal
  weighted contrast and STC a conditional link-scale contrast. Marginal
  contrasts are not generally component-additive on nonlinear links.
- **cML-NMR reporting**: target summaries are average conditional
  link-scale effects at covariate means, not marginally standardized
  effects.

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
