# cpaic component-additive ML-NMR and Stan code review

Status: BLOCK
Recommendation: REQUEST_CHANGES
Date: 2026-06-27

## Skill perspective check

- remove-ai-slops: consulted before judging tests and production complexity. The audited tests violate this perspective in part: several cML-NMR tests are finite-output smoke tests that can pass with statistically wrong or unidentified estimates.
- programming: consulted before judging maintainability and tests. No R/Stan-specific reference is available in that skill, so I applied its general criteria: observable-behavior tests, avoid implementation-mirroring tests, avoid unvalidated boundary data, and reject false confidence from weak tests.

## Coverage reviewed

Reviewed every line of:

- R/cmlnmr.R
- inst/stan/cpaic_binomial.stan
- inst/stan/cpaic_normal.stan
- inst/stan/cpaic_poisson.stan
- man/cmlnmr.Rd
- tests/testthat/test-cmlnmr.R
- documentation/validation/simulation.R

Relevant supporting files reviewed:

- R/data_setup.R, especially build_C_matrix()
- documentation/validation/reproduce_multinma.R
- documentation/validation/VALIDATION.md
- documentation/mathematical-foundations/cpaic-math.qmd ML-NMR section
- DESCRIPTION and NAMESPACE

Verification run:

- `rtk Rscript -e 'devtools::test(filter="cmlnmr", reporter="summary")'`
- Result: completed with no reported failures. This compiles/samples the target Stan models but does not address the high-severity statistical gaps below.

Read-only note:

- Source files were not edited. The worktree already had untracked `docs/`; after the run, `.omo/` is also untracked for evidence artifacts.

## CRITICAL

No CRITICAL findings.

## HIGH

### HIGH 1: Component effects can be unidentified, but cmlnmr() samples and reports them anyway

References:

- R/cmlnmr.R:30-31 claims disconnected subnetworks are connected by construction.
- R/cmlnmr.R:101-111 builds the component design matrix for all treatments.
- R/cmlnmr.R:131-138 sends the design to Stan without a rank or estimability check.
- R/cmlnmr.R:166-179 reports posterior summaries for every beta component.
- R/data_setup.R:22-44 builds component columns but does not check whether observed contrasts identify them.
- man/cmlnmr.Rd:72-77 repeats the connected-by-construction claim.

Problem:

The model estimates one beta per component, but cmlnmr() never checks whether the observed treatment/component design has full rank after accounting for study intercepts. Shared component parameters do not by themselves identify all component effects. Some component combinations only identify sums.

Evidence:

```r
C <- build_C_matrix(c("Placebo", "A+B"), inactive = "Placebo")
#         A B
# Placebo 0 0
# A+B     1 1
qr(C["A+B", , drop = FALSE] - C["Placebo", , drop = FALSE])$rank
# 1, but there are 2 component columns
```

The same problem occurs with no inactive/reference treatment:

```r
C <- build_C_matrix(c("A", "B"), inactive = NULL)
#   A B
# A 1 0
# B 0 1
qr(C["B", , drop = FALSE] - C["A", , drop = FALSE])$rank
# 1, but there are 2 component columns
```

Impact:

`component_effects()` can return posterior means and intervals for component effects that are not identified by the likelihood. Those summaries are then prior- and parameterization-driven, not evidence-driven. This directly undermines the package's component-bridging claim.

Recommended fix:

Before fitting, construct the estimable treatment/component contrast design after accounting for study-specific baselines and check rank against `ncol(C)`. If deficient, stop with the non-estimable components or estimable linear combinations. Require an inactive/reference convention or impose an explicit identifying constraint. Update the docs to say component sharing can identify effects only when the design is full rank.

### HIGH 2: `inactive = NULL` leaves the model without an anchored reference by default

References:

- R/cmlnmr.R:81-88 defaults `inactive = NULL`.
- R/cmlnmr.R:101-111 builds C from all treatment labels.
- R/data_setup.R:30-39 only creates an all-zero row when inactive is supplied.
- man/cmlnmr.Rd:7-31 documents the same default.

Problem:

With the default `inactive = NULL`, every treatment label becomes an active component column. In an arm-based model with study intercepts, absolute component effects are not anchored; only within-study differences are identified. cmlnmr() has no `reference` argument and does not force one component/treatment to zero.

Impact:

Users can call the documented default and receive absolute component summaries that are not interpretable as relative effects. This is especially dangerous because the fit may converge numerically and the tests only check finite estimates for several families.

Recommended fix:

Require `inactive` or a new explicit `reference`/constraint argument for cmlnmr(). If the package intends to support networks without a true inactive component, implement a contrast parameterization and report relative effects or estimable component contrasts only.

### HIGH 3: Aggregate covariate integration assumes independent normal margins, not the joint AgD covariate distribution

References:

- R/cmlnmr.R:7-10 says correlated Gaussian-copula integration is future work.
- R/cmlnmr.R:18-20 maps independent Sobol dimensions through marginal qnorm calls.
- R/cmlnmr.R:121-128 only reads `x_mean` and `x_sd`, with no covariance/correlation input.
- inst/stan/cpaic_binomial.stan:66-75, inst/stan/cpaic_poisson.stan:58-67 average nonlinear inverse-link/rate values over those points.
- documentation/mathematical-foundations/cpaic-math.qmd:399-401 says the aggregate likelihood integrates over the target covariate distribution by quasi-Monte Carlo with a Gaussian copula.

Problem:

For multiple covariates under logit/log links, the integrated population-average outcome depends on the joint distribution, including covariance. The implementation supplies independent normal margins and has no way to pass AgD correlations or an empirical integration distribution. The mathematical documentation overstates this as Gaussian-copula integration.

Evidence:

With two standard normal covariates and the same marginal means/sds, changing only correlation changes the integrated outcome:

```r
# logit example, eta = -0.5 + 1.2*x1 + 1.2*x2
# independent: 0.4176001
# corr 0.8:    0.4306437
# corr -0.8:   0.3908731

# Poisson/log example, eta = -1 + 0.8*x1 + 0.8*x2
# independent: 0.6981654
# corr 0.8:    1.16468
# corr -0.8:   0.4181407
```

Impact:

For multi-effect-modifier AgD, the aggregate likelihood can target the wrong population-average probability/rate/hazard. The bias can be material under the log link.

Recommended fix:

Either restrict cmlnmr() to one aggregate effect modifier unless independence is explicitly accepted, or accept per-study covariance/correlation matrices and generate correlated integration points. The docs should stop claiming Gaussian-copula integration until correlation is actually implemented.

## MEDIUM

### MEDIUM 1: Zero-IPD and aggregate-only shapes are inconsistent between R and Stan

References:

- R/cmlnmr.R:131-136 sets `N_ipd`, `X_ipd`, and `em_idx`.
- inst/stan/cpaic_binomial.stan:43-46, inst/stan/cpaic_normal.stan:33-36, inst/stan/cpaic_poisson.stan:34-38 read `em_idx[1, q]` in transformed data.

Problem:

The Stan files declare `N_ipd` with lower bound 0, but transformed data indexes row 1 of `em_idx` whenever `Q > 0`. On the R side, `matrix(rep(seq_len(Q), each = nrow(ipd)), nrow = nrow(ipd))` becomes a `0 x 0` matrix when `nrow(ipd) == 0`, even if `Q > 0`.

Evidence:

```r
matrix(rep(seq_len(2), each = 0), nrow = 0)
# <0 x 0 matrix>
```

Impact:

Any aggregate-only or empty-IPD call fails at the Stan data/transformed-data boundary rather than with a clear front-door error. If aggregate-only ML-NMR is not intended, the R function and Stan declarations should reject it explicitly.

Recommended fix:

Either require `nrow(ipd) > 0` in R with a clear message, or make `em_idx` a length-Q vector and guard Stan transformed data so it does not index absent IPD rows.

### MEDIUM 2: Runtime diagnostics are not checked before returning an apparently valid fit

References:

- R/cmlnmr.R:160-179 samples, extracts beta draws, and returns summaries.
- tests/testthat/test-cmlnmr.R:43 checks only max Rhat for beta in one binary recovery test.
- tests/testthat/test-cmlnmr.R:70, 86, 102 only check finite component estimates for gaussian, poisson, and survival.

Problem:

cmlnmr() does not inspect divergences, treedepth saturation, E-BFMI, effective sample size, Rhat across all parameters, or CmdStan sampler warnings before returning the object. The tests do not require clean diagnostics for the family smoke tests.

Impact:

The API can return polished component summaries from a fit with serious HMC pathologies. This is particularly risky because weak identifiability and broad priors make such pathologies plausible.

Recommended fix:

After sampling, call CmdStanR diagnostic summaries and either warn or error according to severity. Tests should assert no divergences and acceptable Rhat/ESS for all core parameters in at least one representative model per family.

### MEDIUM 3: Boundary data validation is largely deferred to Stan or accidental R behavior

References:

- R/cmlnmr.R:97-118 reads required columns without explicit missing-column checks.
- R/cmlnmr.R:120-128 accepts `n_int`, means, and sds without domain checks.
- R/cmlnmr.R:141-157 coerces outcomes, counts, sample sizes, exposures, and standard errors without validating finite values or legal ranges.
- inst/stan/cpaic_poisson.stan:16-18 and 34-39 allow zero exposure but take `log(offset_ipd)`.

Problem:

Invalid inputs become `NULL`, `NA`, `NaN`, rejected Stan data, or obscure sampler errors. Examples include negative AgD SDs, binomial `r > n`, non-binary survival status, non-positive SEs, zero exposure with positive counts, and invalid `n_int`.

Impact:

Users get late, low-context failures and, in some cases, a model fit to unintended data after silent coercion.

Recommended fix:

Add R-side validation for required columns and family-specific domains before constructing Stan data. Validate `n_int >= 1`, finite means/sds, `sd >= 0`, binomial counts within sample sizes, positive SEs, valid event indicators, and exposure/time consistency.

### MEDIUM 4: Summaries and generated quantities are too narrow for ML-NMR validation

References:

- inst/stan/cpaic_binomial.stan:89-92 has an empty generated quantities block.
- inst/stan/cpaic_normal.stan and inst/stan/cpaic_poisson.stan have no generated quantities block.
- R/cmlnmr.R:166-179 summarizes only `beta`.

Problem:

The fit object contains raw CmdStan draws, but cmlnmr() only exposes component main effects. It does not compute treatment effects `C %*% beta`, component interactions, relative effects at target covariate values, integrated AgD predictions, log likelihood, or posterior predictive quantities.

Impact:

Users cannot easily validate whether the model fits AgD arms, inspect effect modification, compare target-population relative effects, or run standard Bayesian checks without manually reconstructing model internals.

Recommended fix:

Add generated quantities or R-side posterior summaries for treatment effects, target-population effects, integrated predictions, and log likelihood. Expose gamma summaries separately from main component effects.

## LOW

### LOW 1: The Stan interface advertises P and Q separation, but R always sets P = Q

References:

- inst/stan/cpaic_binomial.stan:18-19 describes P covariates and Q effect modifiers as a subset.
- R/cmlnmr.R:99-110 uses only `effect_modifiers` as X.
- R/cmlnmr.R:131-135 sets `P = Q` and `em_idx` to all columns.

Problem:

The Stan programs are shaped for prognostic covariates plus an effect-modifier subset, but cmlnmr() has no `prognostic_covariates` argument and forces every covariate into both `breg` and treatment interaction terms.

Impact:

This limits practical ML-NMR modeling and can encourage users to include prognostic-only variables as effect modifiers just to adjust baselines.

Recommended fix:

Either simplify the Stan comments/interface to the current design, or add separate prognostic and effect-modifier arguments with a real `em_idx`.

### LOW 2: Validation script is relevant to STC/MAIC, not the cmlnmr Stan implementation

References:

- documentation/validation/simulation.R:23-50 simulates and fits `cnma_bridge()`, `cstc()`, and optionally `cmaic()`.
- documentation/validation/simulation.R never calls `cmlnmr()`.

Problem:

This script supports population-adjusted component methods generally, but it does not validate the audited cmlnmr Stan likelihood, AgD integration, or CmdStan behavior.

Impact:

It should not be cited as evidence that component-additive ML-NMR is correct.

Recommended fix:

Keep the script for STC/MAIC validation, but add separate cmlnmr simulations that test rank-deficient designs, multi-covariate AgD integration, each family likelihood, and sampler diagnostics.

## Major test gaps

- No test expects cmlnmr() to reject rank-deficient component designs.
- No test covers the default `inactive = NULL` identifiability problem.
- No test covers multiple correlated effect modifiers in AgD.
- No test compares cmlnmr integrated AgD likelihood values against known analytic or Monte Carlo expectations.
- No test verifies family-specific invalid inputs are rejected before Stan.
- No test requires clean HMC diagnostics for gaussian, poisson, or survival.
- Survival tests use all events and do not cover censoring/status 0.
- Generated quantities, treatment-level effects, gamma summaries, and log likelihood are untested because they are absent.

## Blockers

1. Add identifiability/rank checking and stop or reparameterize when component effects are not estimable.
2. Require an explicit inactive/reference/constraint convention instead of allowing unanchored defaults.
3. Fix or explicitly restrict multi-covariate AgD integration so the joint covariate distribution is represented correctly.
4. Add diagnostic checks and meaningful tests that fail for the current rank/integration problems.

