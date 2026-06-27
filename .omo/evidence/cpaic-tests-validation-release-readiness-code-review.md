# cpaic Tests, Validation, Documentation, and Release Audit

Date: 2026-06-27
Worktree: `/Users/choxos/Documents/GitHub/cpaic`
Status: BLOCK
Recommendation: REQUEST_CHANGES

## Skill Perspective Check

- `code-review` skill loaded. CodeRabbit was not run because this audit is local/read-only and CodeRabbit would transmit code externally.
- `remove-ai-slops` skill loaded and applied as a review lens. Violations found: tautological dependency tests, implementation-mirroring tests, finite-only smoke tests, and narrative validation without fail conditions.
- `programming` skill loaded and applied as a testing/maintainability lens. Its language-specific rules do not target R, but its test-shape criteria do apply. Violations found: brittle implementation-mirroring tests and missing executable boundary checks for major claims.

## Verdict

The package can build and mostly checks in a relaxed local setup, but it is not release-ready. The blockers are not cosmetic: the public URLs in package metadata are 404, the generated pkgdown site publishes `CLAUDE.md`, `R CMD build` still touches a `.codegraph` socket, the main ML-NMR external validation is assertion-free and contradicts the progress log, and the tests do not prove the strongest public claims.

## CRITICAL

No CRITICAL findings in this lane. The blocking issues below are HIGH severity release and validation risks.

## HIGH

### 1. Package metadata advertises public URLs that currently 404

- Files: `DESCRIPTION:22-23`, `man/cpaic-package.Rd:14-16`
- Evidence: `R CMD check --as-cran --no-manual` on the built tarball reported 404 for `https://choxos.github.io/cpaic/`, `https://github.com/choxos/cpaic`, and `https://github.com/choxos/cpaic/issues`. A direct `curl -L` check also returned 404 for all three.
- Impact: CRAN incoming produces a NOTE, users cannot reach the repository or issue tracker, and the pkgdown workflow is advertising a site that is not live.
- Recommended fix: make the GitHub repo and Pages URL reachable before release, or remove/replace the URLs until they exist.

### 2. A vignette DOI is invalid and triggers CRAN URL checks

- Files: `vignettes/cpaic.bib:12-21`, `documentation/manuscript/cpaic.bib:8-12`
- Evidence: CRAN incoming reported `https://doi.org/10.1002/bimj.202000136` as 404 from `inst/doc/cpaic-methods.html`; direct `curl -L` also returned 404.
- Impact: release checks carry a public reference-quality defect. This is especially bad for a methodology package.
- Recommended fix: correct the DOI if mistyped, or remove the DOI field and cite stable bibliographic metadata.

### 3. The claimed `cmlnmr()` vs `multinma` validation is not an executable PASS

- Files: `documentation/validation/VALIDATION.md:57-73`, `documentation/validation/reproduce_multinma.R:49-62`, `documentation/PROGRESS.md:198-202`
- Evidence: `VALIDATION.md` labels `cmlnmr() vs multinma ML-NMR` as PASS and reports a table. `PROGRESS.md` says the same comparison was partial because the installed `multinma::add_integration()` API differed. The script prints estimates and says they "should be close", but has no `stopifnot()`, `testthat::expect_*()`, tolerance, Rhat check, or exit failure after `multinma` runs.
- Impact: the strongest external Bayesian validation can silently drift and still exit 0. The summary overstates the evidence.
- Recommended fix: either downgrade the validation summary to "informal comparison" or add assertions with explicit tolerances, convergence checks, and a failing exit path.

### 4. The validation and mathematical documentation are local-only despite reproducibility claims

- Files: `README.md:72-76`, `vignettes/cpaic-intro.Rmd:121-126`, `.gitignore:1-2`, `.Rbuildignore:1`, `documentation/validation/VALIDATION.md:3-6`
- Evidence: README and vignette tell users the technical docs and validation study live in `documentation/` or are provided with development sources. The entire `documentation/` tree is ignored by Git and excluded from the package tarball.
- Impact: external users, reviewers, and CI cannot reproduce the validation study from the source distribution or the repository unless the local ignored tree is separately provided.
- Recommended fix: either version the validation scripts and summary, or change public docs to say the validation artifacts are local/private and not part of the development source.

### 5. The generated pkgdown site publishes internal `CLAUDE.md`

- Files: `CLAUDE.md:1-4`, `_pkgdown.yml:40-45`, generated inventory `docs/CLAUDE.md`, `docs/CLAUDE.html`, `docs/search.json`
- Evidence: the generated `docs/` site contains `CLAUDE.md`, `CLAUDE.html`, and search entries for the developer reference. The pkgdown workflow deploys `folder: docs` to `gh-pages` in `.github/workflows/pkgdown.yaml:39-45`.
- Impact: internal agent/developer instructions will be published to the public documentation site if the workflow deploys this tree.
- Recommended fix: prevent pkgdown from copying root `CLAUDE.md`, remove existing generated `docs/CLAUDE*` artifacts, and use `clean: true` or equivalent to remove stale public files.

### 6. CI does not exercise the exported Bayesian flagship

- Files: `tests/testthat/test-cmlnmr.R:1-6`, `tests/testthat/test-cmlnmr.R:46-51`, `.github/workflows/R-CMD-check.yaml:44-52`
- Evidence: both `cmlnmr()` tests are skipped on CRAN and also skipped when CmdStan is absent. The R-CMD-check workflow installs R package dependencies but does not install CmdStan. Local `R CMD check --as-cran` confirmed 2 skipped tests, both from `test-cmlnmr.R`.
- Impact: CI can be green while `cmlnmr()` is broken. The package exports `cmlnmr()` and documents it as the Bayesian flagship, so this is a serious release confidence gap.
- Recommended fix: add a non-CRAN CI lane that installs CmdStan and runs the `cmlnmr()` tests, or mark ML-NMR as experimental and keep its validation outside release claims.

### 7. Tests mirror the implementation instead of validating behavior independently

- Files: `tests/testthat/test-paic.R:15-35`, `tests/testthat/test-paic.R:37-53`, `tests/testthat/test-cmlnmr.R:70`, `tests/testthat/test-cmlnmr.R:86`, `tests/testthat/test-cmlnmr.R:102`
- Evidence: `test-paic.R` manually repeats the same weighted GLM and STC regression formulas used by the implementation, then compares the package to that mirrored implementation. The multi-family ML-NMR test only asserts finite component estimates for gaussian, poisson, and survival.
- Impact: the tests provide false confidence. They do not prove anchoredness, population adjustment, survival interpretation, or agreement with an independent reference.
- Recommended fix: add independent truth-based tests with known data-generating parameters, explicit anchored comparator checks, failure tests for bad inputs, and meaningful posterior/convergence assertions for ML-NMR.

### 8. `R CMD build` still touches `.codegraph/daemon.sock`

- Files: `.Rbuildignore:3`, build log `/var/folders/2h/yqztgsf96gsbkqn60cymgsjr0000gn/T/tmp.1GLggwvbRI/build.log`
- Evidence: `R CMD build` printed `cp: cpaic/.codegraph/daemon.sock is a socket (not copied).` The `.Rbuildignore` pattern is `^\.codegraph$`, which does not match files under `.codegraph/`.
- Impact: hidden local tooling leaks into the build process. The tarball did not include `.codegraph` in this run, but the build warning is a release hygiene defect.
- Recommended fix: use a recursive pattern such as `^\.codegraph(/|$)` and similarly verify `.omo/` and other local tool directories cannot enter build staging.

## MEDIUM

### 9. Default strict check is not reproducible on this machine without missing Suggests

- Files: `DESCRIPTION:42-58`
- Evidence: the built tarball check with default `_R_CHECK_FORCE_SUGGESTS_` failed with `Packages suggested but not available: 'spelling', 'viscomp'`. Both are on CRAN, so this is local library state, not an unavailable dependency.
- Impact: the local release command is not one-command reproducible unless dependencies are installed first.
- Recommended fix: document a release setup command that installs all Suggests, or run checks in a clean dependency manager environment.

### 10. The release check still carries a CRAN incoming NOTE after missing Suggests are relaxed

- Files: `DESCRIPTION:3`, `DESCRIPTION:44`, `DESCRIPTION:59-60`
- Evidence: `_R_CHECK_FORCE_SUGGESTS_=false R CMD check --as-cran --no-manual cpaic_0.0.0.9000.tar.gz` completed with `Status: 1 NOTE`. The NOTE includes new submission, development version `0.0.0.9000`, non-mainstream `cmdstanr` in Suggests with `Additional_repositories`, and invalid URLs.
- Impact: this is not CRAN-ready as-is. Some NOTE items are expected for a development package, but not for a release.
- Recommended fix: use a release version, resolve URLs, and decide whether exporting a function that depends on non-CRAN `cmdstanr` is acceptable for the target release channel.

### 11. Documentation status for `cmlnmr()` is inconsistent

- Files: `README.md:31-32`, `vignettes/cpaic-intro.Rmd:27-31`, `documentation/manual/cpaic-manual.qmd:268-276`, `NAMESPACE:19-31`
- Evidence: README says Phase 2 `cmlnmr()` is in development, while the package exports `cmlnmr()` and the manual documents all four families as supported.
- Impact: users cannot tell whether ML-NMR is production functionality or experimental.
- Recommended fix: choose one status and apply it consistently across README, vignettes, manual, Rd, and validation docs.

### 12. CLAUDE developer docs reference nonexistent public API and stale source layout

- Files: `CLAUDE.md:42`, `CLAUDE.md:106`, `CLAUDE.md:109`, `NAMESPACE:19-31`
- Evidence: CLAUDE mentions `marginal_effects()` and `predict()`, but neither is exported or present in the namespace. It also lists `utils.R, zzz.R`, which are not in the tracked `R/` files.
- Impact: future maintainers and agents will chase APIs/files that do not exist.
- Recommended fix: either implement the documented APIs/files or remove those references.

### 13. Mathematical documentation overstates the integration structure

- Files: `documentation/mathematical-foundations/cpaic-math.qmd:395-405`, `R/cmlnmr.R:5-10`
- Evidence: the math doc says aggregate likelihood integration uses quasi-Monte Carlo with a Gaussian copula. The implementation comment says a Gaussian-copula correlation from IPD is a planned extension; current integration maps independent Sobol margins through normal means/SDs.
- Impact: the math doc claims dependence modeling that the code does not implement.
- Recommended fix: change the math doc to independent normal margins, or implement and test the Gaussian copula.

### 14. Manuscript validation section omits the ML-NMR comparison that `VALIDATION.md` claims

- Files: `documentation/manuscript/cpaic.tex:321-333`, `documentation/validation/VALIDATION.md:57-73`
- Evidence: manuscript says three checks accompany the package and lists bridge, MAIC/Bucher, and simulation. `VALIDATION.md` lists a fourth `cmlnmr()` vs `multinma` check.
- Impact: the manuscript and validation summary disagree on the evidence package.
- Recommended fix: either add the ML-NMR comparison with caveats and assertions, or keep it out of `VALIDATION.md` until it is enforceable.

### 15. Validation scripts depend on undeclared local development tooling

- Files: `documentation/validation/reproduce_netmeta.R:6-8`, `documentation/validation/reproduce_maicplus.R:11-13`, `documentation/validation/reproduce_multinma.R:15-18`, `documentation/validation/simulation.R:15`
- Evidence: scripts use `devtools::load_all()` but `devtools` is not declared in `DESCRIPTION`. The scripts are also ignored by Git, so this is a local-only dependency.
- Impact: a user following `Rscript documentation/validation/<file>.R` may fail before validation starts.
- Recommended fix: either declare a validation environment, use installed-package validation, or keep these scripts explicitly local/private.

### 16. Pkgdown deploy uses `clean: false`

- File: `.github/workflows/pkgdown.yaml:39-45`
- Evidence: the deploy action pushes `folder: docs` with `clean: false`.
- Impact: removed or renamed generated files can remain live on `gh-pages`, including the already-generated `CLAUDE` pages.
- Recommended fix: deploy with cleanup enabled after explicitly excluding private/local files.

## LOW

### 17. `.gitignore` documents `docs/` as ignored, but current state tracks `docs/`

- Files: `.gitignore:14`, tracked generated inventory `docs/`
- Evidence: `git ls-files docs | wc -l` returned 89. The `.gitignore` line contains an inline comment, which Git does not treat as a comment, and in any case tracked files are not ignored.
- Impact: the repository policy for pkgdown output is unclear.
- Recommended fix: decide whether `docs/` is source-controlled or deployment-only, then make `.gitignore`, workflow, and repository state match.

### 18. Generated/local artifacts are large and mostly outside version control

- Files/inventory: `documentation/refs`, `documentation/manual/cpaic-manual.pdf`, `documentation/manuscript/*.aux/*.bbl/*.blg/*.log/*.out/*.pdf`, `.DS_Store`, `data/*.rda`
- Evidence: `documentation/` is 393M; `documentation/refs` is 390M with 273 files; 230 generated/binary-like files were found under `documentation`; `docs/` has 89 files and is 3.2M.
- Impact: acceptable for local research, but not for a clean source release unless explicitly excluded and reproducibility claims are scoped.
- Recommended fix: keep bulk refs and compiled outputs ignored, but version any scripts and small text summaries needed to reproduce claims.

### 19. Generated dates make local documents non-reproducible byte-for-byte

- Files: `documentation/manual/cpaic-manual.qmd:6`, `documentation/mathematical-foundations/cpaic-math.qmd:6`, `documentation/manuscript/cpaic.tex:32`
- Evidence: Quarto docs use `date: today`; manuscript uses `\today`.
- Impact: regenerated PDFs/HTML differ by date even when content is unchanged.
- Recommended fix: use fixed release dates in release artifacts.

## Positive Evidence

- `R CMD build` completed and built `cpaic_0.0.0.9000.tar.gz`.
- With `_R_CHECK_FORCE_SUGGESTS_=false`, `R CMD check --as-cran --no-manual` reached the end: examples OK, `--run-donttest` OK, tests OK, vignettes OK, status 1 NOTE.
- Testthat in that check reported `FAIL 0 | WARN 0 | SKIP 2 | PASS 24`; both skips were `test-cmlnmr.R` on CRAN.
- Validation scripts that do assert conditions passed: `reproduce_netmeta.R`, `reproduce_maicplus.R`, and `simulation.R`.
- `reproduce_multinma.R` ran locally and printed estimates close to the validation summary, but it remains assertion-free.

## Command Evidence

- `rtk zsh -lc 'R CMD build /Users/choxos/Documents/GitHub/cpaic ...; R CMD check --as-cran --no-manual cpaic_0.0.0.9000.tar.gz ...'`
  - Build OK.
  - Default check: `1 ERROR, 1 NOTE`, stopped at missing local Suggests.
- `rtk zsh -lc '_R_CHECK_FORCE_SUGGESTS_=false R CMD check --as-cran --no-manual cpaic_0.0.0.9000.tar.gz ...'`
  - Completed with `Status: 1 NOTE`.
- `rtk Rscript documentation/validation/reproduce_netmeta.R`
  - PASS, max bridge vs `discomb` diff 0.
- `rtk Rscript documentation/validation/reproduce_maicplus.R`
  - PASS, ESS and estimates match manual MAIC/Bucher.
- `rtk Rscript documentation/validation/simulation.R`
  - PASS, output matched `VALIDATION.md`.
- `rtk Rscript documentation/validation/reproduce_multinma.R`
  - Exited 0 and printed estimates, but did not assert any acceptance criteria.
- `rtk Rscript -e 'available.packages(...)'`
  - `spelling`, `viscomp`, `maicplus`, `multinma`, `netmeta`, and `mlumr` are currently on CRAN.
- `rtk zsh -lc 'curl -L ...'`
  - Metadata URLs and DOI above returned 404.

## Coverage Reviewed

Line-audited: `DESCRIPTION`, `.Rbuildignore`, `.github/workflows/*.yaml`, `README.md`, `CLAUDE.md`, `_pkgdown.yml`, `tests/testthat.R`, `tests/testthat/*.R`, `vignettes/*.Rmd`, `vignettes/cpaic.bib`, `documentation/validation/*.R`, `documentation/validation/VALIDATION.md`, `documentation/PROGRESS.md`, `documentation/manual/*.qmd`, `documentation/mathematical-foundations/*.qmd`, `documentation/manuscript/cpaic.tex`, `documentation/manuscript/cpaic.bib`, and `man/*.Rd`.

Inventory-only: `documentation/refs`, generated `docs/`, PDFs, `.rda`, compiled LaTeX outputs, `.DS_Store`, and other bulk/generated/binary artifacts.

## Blockers Before Approval

1. Fix or remove the 404 repository, issue tracker, pkgdown, and DOI links.
2. Stop publishing `CLAUDE.md` in the generated pkgdown site and clean existing generated copies.
3. Make the `cmlnmr()` vs `multinma` validation either executable with assertions or clearly non-PASS/informal.
4. Add a CI path that actually runs the CmdStan-backed `cmlnmr()` tests, or downgrade ML-NMR release claims.
5. Replace implementation-mirroring and finite-only tests with independent behavior/truth tests for core claims.
6. Fix `.Rbuildignore` so `.codegraph/` and other local tool artifacts are never touched by `R CMD build`.
7. Decide whether validation docs are public/reproducible or local/private, then align README, vignettes, `.gitignore`, and `.Rbuildignore`.
