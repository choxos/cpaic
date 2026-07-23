# cpaic 0.0.0.9000

Development software; results may change between versions.

## Data integrity and correctness

* The two-stage methods now fail closed. `cstc()` and `cmaic()` reject a study
  whose regression or weight fit did not converge, is separated, is rank
  deficient, has a degenerate treatment-coefficient covariance, or (for
  `cmaic()`) did not achieve moment balance or whose bootstrap fell below a
  success threshold (`min_boot_success`, default 0.8). An invalid edge is no
  longer passed silently into the additive bridge.
* `cmlnmr()` rejects a study present in both `ipd` and `agd` (which would stack
  two outcome likelihoods for one trial), rejects protected sampler arguments
  such as `data` passed through `...` (which would fit a different dataset than
  the returned object describes), and validates that a supplied covariate
  correlation matrix carrying names is reordered to the effect-modifier order
  before use.
* `.cpaic_replace_contrasts()` enforces a unique {study, treatment-pair} key and
  builds an appended contrast from a typed prototype rather than cloning an
  unrelated aggregate row, so adjusted edges cannot be double-counted or inherit
  stale metadata.
* Front-door validation is completed: missing or non-finite individual outcomes
  and non-finite aggregate covariate means are rejected before reaching Stan;
  `n_int` and `n_boot` must be whole numbers; a supplied `seed` must be a
  non-negative integer within range, and an unset seed is drawn and recorded
  rather than fixed to a constant.
* Construction (`cpaic_network()`) rejects self-comparisons, duplicate
  {study, treatment-pair} rows, and missing treatment labels; `build_C_matrix()`
  rejects empty component tokens and an `inactive` label that matches no
  component. `cnma_bridge()` rejects the empty `common = random = FALSE` model,
  asserts a unique key before fitting, and no longer suppresses substantive
  warnings on non-degenerate networks.
* The Gaussian residual standard deviation has its own prior scale
  (`prior_sigma_sd`), separate from the component-effect prior.
* `cmaic()` restores the caller's random-number state on exit.

## Reporting and diagnostics

* `relative_effects()` for `cmlnmr()` fits reports a `basis` per contrast
  (`"exact"`, `"first-order screen"`, or `"not identified"`), and `cpaic_ranks()`
  excludes elements identified only by aggregate arms unless
  `include_screen_only = TRUE`.
* The automatic MCMC check covers every sampled parameter block and now also
  flags low effective sample size and low E-BFMI, not only divergences,
  tree-depth, and the beta/mu R-hat.
* `prior_sensitivity(prior = "all")` includes the residual and baseline-smoothing
  scales; `component_effects()` accepts a `level`; `loo()` is documented as an
  observation-level diagnostic that does not validate cross-gap prediction.
* Print methods are more precise: a single-covariate contrast is labeled a
  conditional effect at a covariate profile (not a target population), the
  two-stage bridges state that only IPD edges are adjusted, and the bridge output
  notes that the additivity Q statistic cannot test cross-sub-network constancy.

## Documentation

* The README no longer describes an obsolete engine: the survival likelihood is
  analytic with several censoring types and delayed entry, `cmlnmr()` supports
  fixed or random treatment effects, and the simulation figures carry an explicit
  reproducibility caveat.
