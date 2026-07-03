# Set up a (possibly disconnected) component network for cpaic

Builds the network object consumed by
[`cnma_bridge()`](https://choxos.github.io/cpaic/reference/cnma_bridge.md),
[`cmaic()`](https://choxos.github.io/cpaic/reference/cmaic.md),
[`cstc()`](https://choxos.github.io/cpaic/reference/cstc.md) and (Phase
2) [`cmlnmr()`](https://choxos.github.io/cpaic/reference/cmlnmr.md).
Aggregate data are supplied at the **contrast level** (one row per
pairwise comparison: `treat1` vs `treat2`, with a treatment effect `TE`
and standard error `seTE`), matching the input convention of
[`netmeta::discomb()`](https://rdrr.io/pkg/netmeta/man/discomb.html).
Individual patient data (IPD) are optional and used only by the
population-adjusted methods to replace a study's unadjusted contrast
with an adjusted one.

## Usage

``` r
cpaic_network(
  agd,
  ipd = NULL,
  treat1 = "treat1",
  treat2 = "treat2",
  TE = "TE",
  seTE = "seTE",
  studlab = "studlab",
  sm,
  inactive = NULL,
  sep.comps = "+",
  reference = NULL,
  family = NULL,
  ipd_study = ".study",
  ipd_trt = ".trt",
  ipd_outcome = ".y",
  ipd_time = NULL,
  ipd_status = NULL,
  ipd_exposure = NULL,
  ipd_covariates = NULL
)
```

## Arguments

- agd:

  Aggregate (contrast-level) data frame.

- ipd:

  Optional individual patient data frame (one row per patient). Must
  contain a study column, a treatment column, an outcome, and the
  effect-modifier / prognostic covariates.

- treat1, treat2, TE, seTE, studlab:

  Column names in `agd`.

- sm:

  Summary measure (e.g. `"OR"`, `"RR"`, `"MD"`, `"HR"`), passed to
  [`netmeta::discomb()`](https://rdrr.io/pkg/netmeta/man/discomb.html).

- inactive:

  Name of the inactive component / reference (e.g. `"placebo"`).

- sep.comps:

  Component separator in treatment labels. Default `"+"`.

- reference:

  Reference treatment for reported relative effects. Defaults to
  `inactive` when available, otherwise the first treatment.

- family:

  Outcome family for the IPD model, one of `"binomial"`, `"gaussian"`,
  `"poisson"`, `"survival"` (required when `ipd` is given).

- ipd_study, ipd_trt, ipd_outcome:

  Column names in `ipd`. For survival outcomes use
  `ipd_time`/`ipd_status` instead of (or alongside) `ipd_outcome`; for
  Poisson rates an optional `ipd_exposure` offset.

- ipd_time, ipd_status:

  Time and event-indicator column names in `ipd` (survival family).

- ipd_exposure:

  Optional exposure/person-time column in `ipd` (Poisson family); used
  as a log offset.

- ipd_covariates:

  Character vector of covariate column names in `ipd` (the candidate
  effect modifiers / prognostic factors).

## Value

An object of class `cpaic_network`: a list with elements `agd`, `ipd`,
`treatments`, `components`, `C.matrix`, `sm`, `inactive`, `reference`,
`sep.comps`, `family`, and the column-name mappings.

## Details

Treatment labels encode components via `sep.comps` (e.g. `"A + B"`), so
that sub-networks sharing components can be bridged.

## See also

[`cpaic_connectivity()`](https://choxos.github.io/cpaic/reference/cpaic_connectivity.md),
[`cnma_bridge()`](https://choxos.github.io/cpaic/reference/cnma_bridge.md),
[`cmaic()`](https://choxos.github.io/cpaic/reference/cmaic.md),
[`cstc()`](https://choxos.github.io/cpaic/reference/cstc.md)

## Examples

``` r
# Aggregate-only, disconnected network
net <- cpaic_network(cpaic_bin_agd, sm = "OR", inactive = "Placebo")
net
#> cpaic component network
#>   Summary measure:   OR
#>   Treatments:        6
#>   Components:        4 (A, B, C, D)
#>   AgD comparisons:   5
#>   Reference:         Placebo
#>   Inactive:          Placebo
#>   IPD studies:        none
#>   Connected:         FALSE | components bridgeable: TRUE

# With individual patient data for the population-adjusted methods
net_ipd <- cpaic_network(cpaic_bin_agd, ipd = cpaic_bin_ipd, sm = "OR",
                         family = "binomial", ipd_covariates = "x1",
                         inactive = "Placebo")
```
