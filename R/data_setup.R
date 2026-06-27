# Network construction and component coding -----------------------------------

#' Build a component-coded treatment-by-component matrix
#'
#' Splits multi-component treatment labels (e.g. `"A + B"`) on `sep.comps`
#' and returns the binary treatment-by-component design matrix `C`, where
#' `C[t, c] = 1` if treatment `t` contains component `c`. The `inactive`
#' component (e.g. placebo) is represented by an all-zero row, matching the
#' convention of [netmeta::netcomb()].
#'
#' @param treatments Character vector of (unique) treatment labels.
#' @param sep.comps Single character separating components in a treatment
#'   label. Default `"+"`.
#' @param inactive Optional name of the inactive/reference treatment or
#'   component (mapped to a zero row / dropped as a column).
#'
#' @return A binary matrix with one row per treatment and one column per
#'   component (treatments as row names, components as column names).
#' @export
#' @examples
#' build_C_matrix(c("A", "B", "A + B", "placebo"), inactive = "placebo")
build_C_matrix <- function(treatments, sep.comps = "+", inactive = NULL) {
  treatments <- unique(as.character(treatments))
  if (anyNA(treatments) || any(treatments == "")) {
    stop("`treatments` must not contain NA or empty labels.", call. = FALSE)
  }
  split_one <- function(x) trimws(strsplit(x, sep.comps, fixed = TRUE)[[1]])
  comp_list <- lapply(treatments, split_one)

  comps <- unique(unlist(comp_list))
  # Drop the inactive component from the column set; its treatments become
  # zero rows (no active component).
  if (!is.null(inactive)) {
    comps <- setdiff(comps, inactive)
  }
  comps <- comps[order(comps)]

  C <- matrix(0L, nrow = length(treatments), ncol = length(comps),
              dimnames = list(treatments, comps))
  for (i in seq_along(treatments)) {
    active <- intersect(comp_list[[i]], comps)
    if (length(active)) C[i, active] <- 1L
  }
  C
}

#' Set up a (possibly disconnected) component network for cpaic
#'
#' Builds the network object consumed by [cnma_bridge()], [cmaic()],
#' [cstc()] and (Phase 2) `cmlnmr()`. Aggregate data are supplied at the
#' **contrast level** (one row per pairwise comparison: `treat1` vs
#' `treat2`, with a treatment effect `TE` and standard error `seTE`),
#' matching the input convention of [netmeta::discomb()]. Individual
#' patient data (IPD) are optional and used only by the population-adjusted
#' methods to replace a study's unadjusted contrast with an adjusted one.
#'
#' Treatment labels encode components via `sep.comps` (e.g. `"A + B"`), so
#' that sub-networks sharing components can be bridged.
#'
#' @param agd Aggregate (contrast-level) data frame.
#' @param ipd Optional individual patient data frame (one row per patient).
#'   Must contain a study column, a treatment column, an outcome, and the
#'   effect-modifier / prognostic covariates.
#' @param treat1,treat2,TE,seTE,studlab Column names in `agd`.
#' @param sm Summary measure (e.g. `"OR"`, `"RR"`, `"MD"`, `"HR"`), passed
#'   to [netmeta::discomb()].
#' @param inactive Name of the inactive component / reference (e.g.
#'   `"placebo"`).
#' @param sep.comps Component separator in treatment labels. Default `"+"`.
#' @param reference Reference treatment for reported relative effects.
#'   Defaults to `inactive` when available, otherwise the first treatment.
#' @param family Outcome family for the IPD model, one of `"binomial"`,
#'   `"gaussian"`, `"poisson"`, `"survival"` (required when `ipd` is given).
#' @param ipd_study,ipd_trt,ipd_outcome Column names in `ipd`. For survival
#'   outcomes use `ipd_time`/`ipd_status` instead of (or alongside)
#'   `ipd_outcome`; for Poisson rates an optional `ipd_exposure` offset.
#' @param ipd_time,ipd_status Time and event-indicator column names in
#'   `ipd` (survival family).
#' @param ipd_exposure Optional exposure/person-time column in `ipd`
#'   (Poisson family); used as a log offset.
#' @param ipd_covariates Character vector of covariate column names in
#'   `ipd` (the candidate effect modifiers / prognostic factors).
#'
#' @return An object of class `cpaic_network`: a list with elements
#'   `agd`, `ipd`, `treatments`, `components`, `C.matrix`, `sm`,
#'   `inactive`, `reference`, `sep.comps`, `family`, and the column-name
#'   mappings.
#' @seealso [cpaic_connectivity()], [cnma_bridge()], [cmaic()], [cstc()]
#' @examples
#' # Aggregate-only, disconnected network
#' net <- cpaic_network(cpaic_bin_agd, sm = "OR", inactive = "Placebo")
#' net
#'
#' # With individual patient data for the population-adjusted methods
#' net_ipd <- cpaic_network(cpaic_bin_agd, ipd = cpaic_bin_ipd, sm = "OR",
#'                          family = "binomial", ipd_covariates = "x1",
#'                          inactive = "Placebo")
#' @export
cpaic_network <- function(agd, ipd = NULL,
                          treat1 = "treat1", treat2 = "treat2",
                          TE = "TE", seTE = "seTE", studlab = "studlab",
                          sm, inactive = NULL, sep.comps = "+",
                          reference = NULL, family = NULL,
                          ipd_study = ".study", ipd_trt = ".trt",
                          ipd_outcome = ".y", ipd_time = NULL,
                          ipd_status = NULL, ipd_exposure = NULL,
                          ipd_covariates = NULL) {
  if (missing(sm) || is.null(sm)) {
    stop("`sm` (summary measure, e.g. \"OR\", \"MD\", \"HR\") is required.",
         call. = FALSE)
  }
  agd <- as.data.frame(agd)
  req <- c(treat1, treat2, TE, seTE, studlab)
  miss <- setdiff(req, names(agd))
  if (length(miss)) {
    stop("`agd` is missing column(s): ", paste(miss, collapse = ", "),
         call. = FALSE)
  }

  treatments <- sort(unique(c(as.character(agd[[treat1]]),
                              as.character(agd[[treat2]]))))
  C <- build_C_matrix(treatments, sep.comps = sep.comps, inactive = inactive)
  components <- colnames(C)

  if (is.null(reference)) {
    reference <- if (!is.null(inactive) && inactive %in% treatments) {
      inactive
    } else {
      treatments[1]
    }
  } else if (!reference %in% treatments) {
    stop("`reference` (\"", reference, "\") is not a treatment in `agd`.",
         call. = FALSE)
  }

  ipd_info <- NULL
  if (!is.null(ipd)) {
    ipd <- as.data.frame(ipd)
    if (is.null(family)) {
      stop("`family` is required when `ipd` is supplied.", call. = FALSE)
    }
    family <- match.arg(family,
                        c("binomial", "gaussian", "poisson", "survival"))
    req_ipd <- c(ipd_study, ipd_trt)
    miss_ipd <- setdiff(req_ipd, names(ipd))
    if (length(miss_ipd)) {
      stop("`ipd` is missing column(s): ", paste(miss_ipd, collapse = ", "),
           call. = FALSE)
    }
    if (is.null(ipd_covariates)) {
      stop("`ipd_covariates` must name the covariate columns in `ipd`.",
           call. = FALSE)
    }
    miss_cov <- setdiff(ipd_covariates, names(ipd))
    if (length(miss_cov)) {
      stop("`ipd` is missing covariate column(s): ",
           paste(miss_cov, collapse = ", "), call. = FALSE)
    }
    if (family == "survival") {
      if (is.null(ipd_time) || is.null(ipd_status)) {
        stop("survival family requires `ipd_time` and `ipd_status`.",
             call. = FALSE)
      }
      miss_surv <- setdiff(c(ipd_time, ipd_status), names(ipd))
      if (length(miss_surv)) {
        stop("`ipd` is missing survival column(s): ",
             paste(miss_surv, collapse = ", "), call. = FALSE)
      }
    }
    ipd_studies <- unique(as.character(ipd[[ipd_study]]))
    ipd_info <- list(study = ipd_study, trt = ipd_trt, outcome = ipd_outcome,
                     covariates = ipd_covariates, studies = ipd_studies)
  }

  structure(
    list(
      agd = agd,
      ipd = ipd,
      treatments = treatments,
      components = components,
      C.matrix = C,
      sm = sm,
      inactive = inactive,
      reference = reference,
      sep.comps = sep.comps,
      family = family,
      cols = list(treat1 = treat1, treat2 = treat2, TE = TE, seTE = seTE,
                  studlab = studlab, ipd_time = ipd_time,
                  ipd_status = ipd_status, ipd_exposure = ipd_exposure),
      ipd_info = ipd_info
    ),
    class = "cpaic_network"
  )
}

#' @export
print.cpaic_network <- function(x, ...) {
  cat("cpaic component network\n")
  cat("  Summary measure:   ", x$sm, "\n", sep = "")
  cat("  Treatments:        ", length(x$treatments), "\n", sep = "")
  cat("  Components:        ", length(x$components),
      " (", paste(utils::head(x$components, 8), collapse = ", "),
      if (length(x$components) > 8) ", ..." else "", ")\n", sep = "")
  cat("  AgD comparisons:   ", nrow(x$agd), "\n", sep = "")
  cat("  Reference:         ", x$reference, "\n", sep = "")
  if (!is.null(x$inactive)) cat("  Inactive:          ", x$inactive, "\n", sep = "")
  if (!is.null(x$ipd)) {
    cat("  IPD studies:       ", length(x$ipd_info$studies),
        " (", x$family, "; ", nrow(x$ipd), " patients)\n", sep = "")
  } else {
    cat("  IPD studies:        none\n")
  }
  conn <- tryCatch(cpaic_connectivity(x), error = function(e) NULL)
  if (!is.null(conn)) {
    cat("  Connected:         ", conn$connected,
        " | components bridgeable: ", conn$identifiable, "\n", sep = "")
  }
  invisible(x)
}
