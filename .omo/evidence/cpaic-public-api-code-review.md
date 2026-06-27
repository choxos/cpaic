# cpaic Public API Code Review

Goal: strict audit of the cpaic R package public API, data model, component coding, connectivity, effects, diagnostics, and plotting layers.

Scope reviewed:
- Source: `R/data_setup.R`, `R/connectivity.R`, `R/bridge.R`, `R/effects.R`, `R/diagnostics.R`, `R/plot.R`, `R/data.R`, `R/cpaic-package.R`
- Namespace: `NAMESPACE`
- Docs: `man/build_C_matrix.Rd`, `man/cpaic_network.Rd`, `man/cpaic_connectivity.Rd`, `man/cnma_bridge.Rd`, `man/component_effects.Rd`, `man/relative_effects.Rd`, `man/league_table.Rd`, `man/additivity_test.Rd`, `man/effective_sample_size.Rd`, `man/forest.Rd`, `man/plot.cpaic_network.Rd`, `man/cpaic_bin_agd.Rd`, `man/cpaic_bin_ipd.Rd`, `man/cpaic-package.Rd`, and shallow interaction docs for `man/cmaic.Rd`, `man/cstc.Rd`, `man/cmlnmr.Rd`
- Tests: `tests/testthat/test-bridge.R`, `tests/testthat/test-package.R`

Skill-perspective check:
- `omo:remove-ai-slops`: loaded and applied as an overfit/slop review over tests and production code. No deletion-only tests or production over-extraction were found. Violations found: low-signal dependency-only package test, and missing tests around core bridge failure modes.
- `omo:programming`: loaded. It has no R-specific reference files, so its general criteria were applied: strict boundary parsing, useful tests over implementation-mirroring tests, no untyped escape hatches, and avoiding needless complexity. Violations found: invalid data model states are accepted at the public boundary, and tests miss observable failure behavior.

Verification:
- `rtk Rscript -e 'testthat::test_local(filter = "bridge|package", load_package = "source")'`: PASS, 18 passed.
- Runtime probe: non-identifiable disconnected network had `rank(X) = 2 / 4` and `identifiable = FALSE`, but `cnma_bridge(..., common = TRUE, random = FALSE)` still returned a `cpaic_bridge` object with all component estimates `NA`.
- Runtime probe: `cnma_bridge(cpaic_bin_agd, sm = "OR", inactive = "Placebo")` failed with `inherits(network, "cpaic_network") is not TRUE`, despite documentation claiming aggregate data frame input.
- Runtime probe: `cpaic_network()` accepted IPD with the `.y` outcome column removed and printed a valid-looking binomial IPD network.
- Runtime probe: `relative_effects(br, reference = "NotATreatment")` failed with raw `subscript out of bounds`.
- Runtime probe: bundled data dimensions and names matched their Rd files.
- Worktree note: `?? docs/` was already present before report writing. This review added only this report artifact under `.omo/evidence/`.

## CRITICAL

No CRITICAL findings.

## HIGH

1. `cnma_bridge()` does not enforce component identifiability before fitting.
   - Location: `R/bridge.R:32`, `R/bridge.R:37-63`, `R/connectivity.R:75-77`, `tests/testthat/test-bridge.R:50-64`
   - Problem: `cpaic_connectivity()` computes the rank-based identifiability flag, and the test suite verifies a non-bridgeable network is unidentifiable, but `cnma_bridge()` never calls that check or rejects such a network. A runtime probe using the same non-bridgeable shape returned a bridge object with `NA` component effects.
   - Why it matters: The package goal is component-based reconnection of disconnected networks only when components identify the bridge. Returning a valid-looking `cpaic_bridge` object for a non-identifiable network lets downstream `component_effects()`, `relative_effects()`, `league_table()`, and diagnostics operate on an invalid fit.
   - Recommended fix: In `cnma_bridge()`, call `cpaic_connectivity(network)` before `netmeta::discomb()`. Stop with a clear error whenever `identifiable` is false, including connected but rank-deficient component codings. Store the connectivity object on successful bridge fits. Add a regression test that the non-bridgeable fixture in `test-bridge.R:50-64` causes `cnma_bridge()` to error before fitting.

2. `cpaic_network()` accepts invalid IPD data models that contradict its public contract.
   - Location: `R/data_setup.R:61-63`, `R/data_setup.R:144-158`, `R/data_setup.R:170-172`, `man/cpaic_network.Rd:32-34`
   - Problem: The docs say IPD must contain study, treatment, outcome, and covariates. The implementation validates only study, treatment, and covariates for non-survival families. It does not validate `ipd_outcome` for binomial, gaussian, or poisson families, and does not validate a supplied `ipd_exposure` column. A runtime probe removed `.y` from `cpaic_bin_ipd`; `cpaic_network(..., family = "binomial")` still returned and printed a valid-looking network.
   - Why it matters: This admits impossible `cpaic_network` objects at the package boundary. The PAIC layers then receive malformed state and fail later with less local, harder-to-debug errors, or may accidentally operate on the wrong columns if names drift.
   - Recommended fix: Validate family-specific IPD columns in `cpaic_network()`: `ipd_outcome` for binomial, gaussian, and poisson; `ipd_time` and `ipd_status` for survival; and `ipd_exposure` whenever it is non-null. Also validate that required columns are not `NULL`, not missing from `ipd`, and have plausible storage mode for the selected family.

## MEDIUM

1. `cnma_bridge()` documentation advertises raw aggregate data frame input that the implementation rejects.
   - Location: `R/bridge.R:15-16`, `R/bridge.R:32-35`, `man/cnma_bridge.Rd:10-11`
   - Problem: The Rd and roxygen say `network` may be a `cpaic_network()` object or an aggregate contrast data frame, but line 33 immediately requires `inherits(network, "cpaic_network")`. A runtime probe with `cpaic_bin_agd` failed with `inherits(network, "cpaic_network") is not TRUE`.
   - Why it matters: This is a public API contract mismatch. Users following the help page cannot call the function as documented.
   - Recommended fix: Either implement the documented data-frame path by forwarding to `cpaic_network()` with column-name arguments from `...`, or remove that claim from roxygen/Rd and make the error message explicit.

2. `cnma_bridge()` accepts nonsensical effect-model flags and silently reports common effects.
   - Location: `R/bridge.R:17-19`, `R/bridge.R:48-49`, `R/bridge.R:62-63`
   - Problem: `common` and `random` can both be `FALSE`. The wrapper still sets `effect <- "common"` whenever `random` is false and then builds a component table from common/fixed fields. A runtime probe with `common = FALSE, random = FALSE` returned and printed common effects.
   - Why it matters: A user explicitly requested no common and no random model, but the wrapper returns a model anyway. This undermines the meaning of the exported arguments and can hide misconfigured analyses.
   - Recommended fix: Reject `!common && !random` before calling `netmeta::discomb()`. If both are true, document which table `component_effects()` will return or add an `effect` selector.

3. `relative_effects()` lacks public validation for `reference`, producing a raw subscript error.
   - Location: `R/effects.R:39-49`, `R/effects.R:53-75`, `man/relative_effects.Rd:19`
   - Problem: When a user supplies a reference treatment absent from `fit$TE.*`, `build()` indexes `TE[t1, t2]` directly and throws `subscript out of bounds`.
   - Why it matters: This is an exported API. A bad reference should be diagnosed at the boundary with available treatment names, not as an internal matrix-indexing failure.
   - Recommended fix: After deriving `trts`, validate `reference %in% trts` and stop with a clear message. Also validate `level` is a scalar numeric in `(0, 1)`.

4. Tests cover the connectivity flag but not the bridge safety contract.
   - Location: `tests/testthat/test-bridge.R:20-48`, `tests/testthat/test-bridge.R:50-64`
   - Problem: The tests prove one bridgeable fixture and one non-identifiable connectivity result, but they never assert that `cnma_bridge()` refuses non-identifiable networks or returns finite effects only for valid networks.
   - Why it matters: The current suite passes while the main bridge API can produce an invalid bridge object for a non-identifiable network.
   - Recommended fix: Add tests for `expect_error(cnma_bridge(non_bridgeable_net), "not identifiable|cannot be bridged")`, finite component tables on bridgeable fits, and `relative_effects()` behavior after a valid bridge.

## LOW

1. `component_effects()` documentation omits supported `cpaic_mlnmr` dispatch and leaves scale/back-transformation ambiguous.
   - Location: `R/bridge.R:121-125`, `NAMESPACE:5-7`, `R/cmlnmr.R:204`, `man/component_effects.Rd:10-15`
   - Problem: The generic docs list `cpaic_bridge`, `cpaic_maic`, and `cpaic_stc`, but `NAMESPACE` exports `component_effects.cpaic_mlnmr`. For cNMA bridge fits, `component_table()` sets `attr(out, "backtransf") <- isTRUE(fit$backtransf)` while returning link-scale component estimates.
   - Why it matters: Users cannot infer from the help page which classes are supported or whether component estimates are link-scale/log-scale or natural-scale.
   - Recommended fix: Update the generic docs to include `cpaic_mlnmr`, and explicitly state the scale of component effects for log measures. Avoid setting a `backtransf` attribute that suggests transformed values unless values are actually transformed.

2. `tests/testthat/test-package.R` is dependency-only smoke coverage.
   - Location: `tests/testthat/test-package.R:1-5`
   - Problem: The test name says package loads and engine dependencies are available, but it only checks `requireNamespace()` for `netmeta`, `maicplus`, and `multinma`. It does not call `library(cpaic)` or a minimal exported function.
   - Why it matters: It can pass without proving any package API behavior. It is not harmful, but it provides weak confidence relative to the stated purpose.
   - Recommended fix: Either rename it to dependency availability, or add a true package-load/export smoke test such as constructing a tiny `cpaic_network()` object after `load_package = "source"`.

## No-Finding Notes

- `build_C_matrix()` basic component splitting, inactive zero-row handling, and bundled test expectations are coherent for the documented examples.
- `cpaic_connectivity()` rank computation and returned matrix fields are internally consistent for the reviewed fixtures; the defect is that `cnma_bridge()` does not enforce its result.
- `additivity_test()` and `effective_sample_size()` are small wrappers with clear dispatch or class checks in the scoped code.
- `plot.cpaic_network()` and `forest()` are simple base plotting wrappers with no direct correctness defect found in static review, though they have no scoped tests.
- `R/data.R` documentation matches loaded object dimensions, names, and truth attributes in the current checkout.
- `NAMESPACE` is coherent for the scoped exports and S3 methods, aside from the documentation mismatch around `component_effects.cpaic_mlnmr`.

## Status

codeQualityStatus: BLOCK
recommendation: REQUEST_CHANGES
blockers:
- Enforce bridge identifiability in `cnma_bridge()` before fitting.
- Validate family-specific IPD columns in `cpaic_network()` before constructing a `cpaic_network` object.

