# cpaic PAIC adjustment layer code review

Goal: strict read-only statistical/code audit of cpaic's PAIC adjustment layer, cMAIC and cSTC.

Status: BLOCK

Recommendation: REQUEST_CHANGES

Scope reviewed:
- `R/cmaic.R`: every line
- `R/cstc.R`: every line
- `man/cmaic.Rd`: every line
- `man/cstc.Rd`: every line
- `tests/testthat/test-paic.R`: every line
- `documentation/validation/reproduce_maicplus.R`: every line
- `documentation/validation/reproduce_multinma.R`: every line; relevant only as a cML-NMR interaction/validation comparator, not as evidence for cMAIC/cSTC
- Interaction points: `R/data_setup.R`, `R/bridge.R`, `R/effects.R`, `documentation/validation/simulation.R`, `documentation/validation/VALIDATION.md`

Skill-perspective check:
- `omo:remove-ai-slops` was loaded and applied as a read-only overfit/slop review pass. The diff/test surface violates that perspective through narrow, implementation-mirroring tests that create false confidence around OR-only two-arm behavior while documented behavior is broader.
- `omo:programming` was loaded and applied as a maintainability/test-shape perspective. The PAIC layer violates that perspective at input boundaries: summary-measure/family contracts and target values are not parsed into a valid internal contract before model fitting.

Verification run:
- `rtk Rscript -e 'devtools::test(filter = "paic", reporter = "summary")'`
  - Result: PASS, reported `paic: ......`
- Targeted probes:
  - cSTC with `family = "binomial", sm = "RR"` stored `TE = 0.538997`; expected log RR was `0.405465`; the stored value matched logistic log OR `0.538997`.
  - cMAIC with `family = "binomial", sm = "RR"` stored the same log OR instead of log RR.
  - cMAIC with `target_sd = c(x = 0.1)` and `target_sd = c(x = 10)` produced identical `TE = 0.563136993427` and `ESS = 196.023995051912`.
  - Three-arm IPD studies failed in both cMAIC and cSTC with netmeta multi-arm comparison errors when all pairwise rows were required.
  - `target = c(x = NA_real_)` in cSTC reached a cryptic downstream model error, not a target-validation error.
  - `n_boot = 0` in cMAIC produced missing `TE`/`seTE` downstream and no explicit bootstrap-argument error.

## CRITICAL

None.

## HIGH

### 1. cMAIC and cSTC ignore `sm` and can inject the wrong estimand into `netmeta`

Files/lines:
- `R/cmaic.R:28-38`
- `R/cmaic.R:64-79`
- `R/cmaic.R:202-204`
- `R/cstc.R:41-52`
- `R/cstc.R:146-148`
- `R/data_setup.R:65-73`
- `man/cpaic_network.Rd:38-49`

Problem:
`cpaic_network()` advertises summary measures such as OR, RR, MD, and HR. The PAIC helpers receive `network$sm`, but the model-fitting code does not use it. For binomial outcomes, both layers always use `stats::binomial()` with the default logit link. That returns log odds ratios even when the network summary measure is `RR` or any other binomial measure.

Evidence:
For a simple A vs C binary dataset with risks 0.30 and 0.20:
- Expected `sm = "RR"` contrast: `log(0.30 / 0.20) = 0.405465`
- Logistic OR contrast: `qlogis(0.30) - qlogis(0.20) = 0.538997`
- Both `cstc()` and `cmaic()` stored `0.538997` in the adjusted AgD contrast while the network said `sm = "RR"`.

Impact:
This mixes scales in the same `netmeta::discomb()` call. A log OR generated from IPD can be combined with aggregate log RR evidence under `sm = "RR"`, producing numerically valid-looking but statistically invalid component effects and relative effects. This is a blocker for the documented binary-outcome surface outside OR.

Recommended fix:
Parse and enforce valid `family`/`sm` contracts before fitting. Either:
- support only combinations the code actually estimates, for example binomial OR, gaussian MD, poisson IRR/rate ratio, survival HR, and error clearly for unsupported combinations; or
- choose model/link and contrast extraction to match `sm`, for example RR/log link and RD/identity link where statistically appropriate.

Add regression tests proving the adjusted `TE` scale equals the requested `sm`, not just the family default.

### 2. `target_sd` is documented as second-moment matching but is a no-op

Files/lines:
- `R/cmaic.R:49-57`
- `R/cmaic.R:70-72`
- `R/cmaic.R:93-95`
- `R/cmaic.R:135-136`
- `R/cmaic.R:174`
- `man/cmaic.Rd:28-29`

Problem:
`.cpaic_center()` creates `<em>_sq_CENTERED` columns when `target_sd` is supplied, but `em_centered_cols` is always `paste0(effect_modifiers, "_CENTERED")`. The square-centered columns are never passed to `maicplus::estimate_weights()` in either the point estimate or bootstrap paths.

Evidence:
On a Gaussian two-arm IPD example, extreme target SD inputs produced identical results:
- `target_sd = c(x = 0.1)`: `TE = 0.563136993427`, `ESS = 196.023995051912`
- `target_sd = c(x = 10)`: `TE = 0.563136993427`, `ESS = 196.023995051912`

Impact:
Users are told second moments are matched when `target_sd` is supplied, but the fitted weights ignore them. This can materially change MAIC weights, effective sample size, and adjusted effects when variance imbalance matters.

Recommended fix:
When `target_sd` is supplied for an effect modifier, include the generated square-centered column in `centered_colnames` for both the main fit and every bootstrap refit. Validate `target_sd` as finite, non-negative, scalar numeric input. Add tests that weighted first and second moments match the target and that changing `target_sd` can change weights/ESS.

### 3. Multi-arm IPD studies are not correctly adjusted or rejected

Files/lines:
- `R/cmaic.R:83-111`
- `R/cmaic.R:191-204`
- `R/cmaic.R:234-251`
- `R/cstc.R:57-65`
- `R/cstc.R:135-148`
- `R/data_setup.R:51-55`

Problem:
For each IPD study, cMAIC/cSTC choose one reference arm from the first matching AgD row and emit only non-reference-vs-reference contrasts. `netmeta` contrast-level multi-arm studies require the complete set of pairwise comparisons for a k-arm study. The replacement helper then updates only exact-orientation matches and leaves any other rows untouched or appends new rows.

Evidence:
For a three-arm IPD study with A, B, and C:
- With all three AgD pairwise rows present but `TE`/`seTE = NA`, both cMAIC and cSTC fail downstream with:
  `After removing comparisons with missing treatment effects or standard errors, study 'S1' has a wrong number of comparisons.`
- With non-missing placeholder contrasts, `netmeta` detects inconsistent multi-arm treatment estimates rather than receiving a complete adjusted multi-arm set.

Impact:
The layer is unsafe for multi-arm IPD studies. It can fail cryptically, and with real unadjusted values it risks mixing adjusted and unadjusted within-study contrasts. Variance/covariance handling for multi-arm adjusted contrasts is also absent.

Recommended fix:
Either explicitly reject IPD studies with more than two arms before fitting, with a targeted error message, or implement complete multi-arm adjustment. A complete implementation must replace the entire study's contrast set and handle multi-arm covariance consistently with `netmeta` expectations.

## MEDIUM

### 4. Target population inputs are not parsed or validated before model fitting

Files/lines:
- `R/cmaic.R:168-172`
- `R/cmaic.R:49-57`
- `R/cstc.R:117-121`
- `R/cstc.R:19-22`

Problem:
The code checks that each effect modifier has a non-null target entry, but it does not require a finite numeric scalar. `NA`, `NaN`, `Inf`, character values, or length > 1 entries can reach centering and model fitting.

Evidence:
`cstc(..., target = c(x = NA_real_))` failed with the downstream error `contrasts can be applied only to factors with 2 or more levels`, not with a clear target-population validation error.

Impact:
Bad target inputs produce misleading model errors or silent recycling. That weakens failure-mode clarity and makes statistical misuse harder to detect.

Recommended fix:
Parse `target` into a named numeric scalar vector/list at the public boundary. Require every effect modifier to have exactly one finite numeric mean. Apply the same validation to `target_sd`.

### 5. Bootstrap failure handling can produce invalid or opaque SEs

Files/lines:
- `R/cmaic.R:84-106`
- `R/cmaic.R:137-139`
- `R/cmaic.R:157-158`

Problem:
`n_boot` is not validated. Bootstrap weight/model errors are swallowed, failed replicates remain `NA`, and standard errors are computed with `na.rm = TRUE` without tracking how many usable bootstrap estimates remain per contrast.

Evidence:
`n_boot = 0` reached `netmeta` with missing estimates/SEs and failed with `No network meta-analysis possible`, rather than a direct cMAIC argument error.

Impact:
All-failed or mostly failed bootstraps can silently produce `NA` or unstable `seTE` values until a downstream failure, or worse, produce SEs from too few replicates without warning.

Recommended fix:
Require `n_boot` to be a positive integer. Count successful bootstrap estimates per contrast, expose/report the count, and fail or warn when the usable replicate count is below a defensible threshold.

### 6. cSTC documentation overstates G-computation/delta-method behavior

Files/lines:
- `R/cstc.R:19-31`
- `R/cstc.R:57-65`
- `R/cstc.R:82-85`
- `man/cstc.Rd:43-47`

Problem:
The implementation centers effect modifiers and takes the main treatment coefficient as the adjusted contrast. It does not perform explicit G-computation over a target covariate distribution, and it does not implement visible delta-method machinery beyond using the fitted coefficient covariance.

Impact:
For nonlinear or noncollapsible models, this is a conditional contrast at supplied target means, not necessarily a marginal target-population average contrast. The current wording can lead users to overinterpret the estimand.

Recommended fix:
Either implement the claimed G-computation/delta-method estimator, including a clear target distribution contract, or revise the documentation to say the estimator is the conditional anchored STC contrast at target covariate means.

## LOW

### 7. Formula construction is brittle for non-syntactic covariate names

Files/lines:
- `R/cstc.R:24-31`
- `R/cstc.R:34-50`

Problem:
cSTC builds formulas by pasting raw column names. Non-syntactic effect-modifier or prognostic names such as `age group` or `x-1` can misparse or fail.

Impact:
This is a preventable failure mode for user-provided IPD column names.

Recommended fix:
Use `reformulate()` where possible or quote/backtick variable names safely before constructing formulas.

## Test and validation gaps

Existing tests:
- `tests/testthat/test-paic.R` covers one two-arm binomial OR cMAIC reproduction, one two-arm binomial OR cSTC reproduction, and one cSTC-vs-naive movement check.
- `documentation/validation/reproduce_maicplus.R` repeats the same two-arm binomial OR cMAIC path.
- `documentation/validation/simulation.R` gives useful Gaussian cSTC/cMAIC simulation evidence, but it is not in the requested unit-test file and does not cover binary non-OR links.
- `documentation/validation/reproduce_multinma.R` validates cML-NMR, not cMAIC/cSTC.

Major missing coverage:
- `family = "binomial"` with `sm = "RR"` and `sm = "RD"` or explicit rejection of unsupported binary measures.
- Any `family`/`sm` mismatch rejection.
- `target_sd` second-moment matching.
- Multi-arm IPD studies: explicit rejection or complete adjusted contrast replacement.
- Multiple IPD studies and orientation-sensitive replacement.
- Poisson with exposure, survival HR, and Gaussian coverage in automated tests.
- Bootstrap edge cases: `n_boot = 0`, failed bootstrap replicates, valid replicate counts, and SE non-missingness.
- Target validation for missing, non-finite, non-numeric, duplicate, and length > 1 target entries.
- Tests that assert observable statistical contracts rather than only reproducing the same implementation formula.

## Blockers

Approval blockers:
- Correct or explicitly reject unsupported `sm`/`family` combinations so adjusted IPD contrasts are never combined on the wrong scale.
- Make `target_sd` real second-moment matching or remove/document it as unsupported.
- Reject or correctly implement multi-arm IPD adjustment before calling `netmeta`.
- Add tests for the above contracts. The current tests are too narrow and overfit the happy path.

