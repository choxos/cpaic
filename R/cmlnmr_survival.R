# Exact survival likelihood support for component-additive ML-NMR
#
# The M-spline hazard basis, integrated basis, censoring contributions, delayed
# entry adjustment, and likelihood-level AgD integration are ported from
# multinma (Phillippo et al. 2020). Both packages are licensed under GPL-3.

#' Define a piecewise or continuous M-spline baseline basis
#' @noRd
.cpaic_survival_basis_spec <- function(observed_times,
                                       baseline = c("piecewise", "mspline"),
                                       cut_points = numeric(), n_basis = 6L) {
  baseline <- match.arg(baseline)
  observed_times <- as.numeric(observed_times)
  if (!length(observed_times) || any(!is.finite(observed_times)) ||
      any(observed_times <= 0)) {
    stop("Survival outcome times must be positive and finite.", call. = FALSE)
  }

  if (baseline == "piecewise") {
    cut_points <- as.numeric(cut_points)
    if (length(cut_points) &&
        (is.unsorted(cut_points, strictly = TRUE) ||
         any(!is.finite(cut_points)) || any(cut_points <= 0) ||
         any(cut_points >= max(observed_times)))) {
      stop("Piecewise `cut_points` must be strictly increasing, positive, ",
           "finite, and below the largest observed time.", call. = FALSE)
    }
    return(structure(
      list(baseline = baseline, cut_points = cut_points,
           n_basis = length(cut_points) + 1L),
      class = "cpaic_survival_basis"
    ))
  }

  if (!requireNamespace("splines2", quietly = TRUE)) {
    stop("Package 'splines2' is required for baseline = 'mspline'.",
         call. = FALSE)
  }
  if (length(n_basis) != 1L || !is.finite(n_basis) ||
      n_basis != as.integer(n_basis) || n_basis < 4L) {
    stop("`n_basis` must be an integer of at least 4 for a cubic M-spline.",
         call. = FALSE)
  }
  boundary <- c(0, max(observed_times))
  basis <- tryCatch(
    splines2::mSpline(
      observed_times, df = as.integer(n_basis), degree = 3L,
      intercept = TRUE, Boundary.knots = boundary
    ),
    error = function(e) {
      stop("Could not construct the continuous M-spline basis: ",
           conditionMessage(e), call. = FALSE)
    }
  )
  structure(
    list(baseline = baseline, basis = basis,
         n_basis = ncol(basis), boundary = boundary),
    class = "cpaic_survival_basis"
  )
}

#' Evaluate hazard and cumulative-hazard bases
#' @noRd
.cpaic_survival_basis_eval <- function(spec, time, start_time, entry_time) {
  stopifnot(inherits(spec, "cpaic_survival_basis"))
  time <- as.numeric(time)
  start_time <- rep_len(as.numeric(start_time), length(time))
  entry_time <- rep_len(as.numeric(entry_time), length(time))
  if (any(!is.finite(time)) || any(time <= 0) ||
      any(!is.finite(start_time)) || any(start_time < 0) ||
      any(!is.finite(entry_time)) || any(entry_time < 0)) {
    stop("Survival times must be finite; outcome times must be positive and ",
         "start or entry times must be nonnegative.", call. = FALSE)
  }

  if (spec$baseline == "piecewise") {
    cuts <- spec$cut_points
    K <- spec$n_basis
    hazard_at <- function(x) {
      idx <- findInterval(x, cuts) + 1L
      out <- matrix(0, nrow = length(x), ncol = K)
      out[cbind(seq_along(x), idx)] <- 1
      out
    }
    integrated_at <- function(x) {
      lower <- c(0, cuts)
      upper <- c(cuts, Inf)
      out <- vapply(seq_len(K), function(k) {
        pmax(pmin(x, upper[k]) - lower[k], 0)
      }, numeric(length(x)))
      matrix(out, nrow = length(x), ncol = K)
    }
  } else {
    as_basis_matrix <- function(x, integral) {
      out <- stats::update(spec$basis, x = x, integral = integral)
      matrix(as.numeric(out), nrow = length(x), ncol = spec$n_basis)
    }
    hazard_at <- function(x) as_basis_matrix(x, integral = FALSE)
    integrated_at <- function(x) as_basis_matrix(x, integral = TRUE)
  }

  list(
    time = hazard_at(time),
    itime = integrated_at(time),
    start_itime = integrated_at(start_time),
    entry_itime = integrated_at(entry_time),
    delayed = as.integer(entry_time > 0)
  )
}

#' R reference implementation of the exact survival log likelihood
#' @noRd
.cpaic_survival_loglik <- function(basis, status, eta, coefficients) {
  status <- as.integer(status)
  eta <- rep_len(as.numeric(eta), length(status))
  coefficients <- as.numeric(coefficients)
  if (any(!status %in% 0:3) || any(!is.finite(eta)) ||
      any(!is.finite(coefficients)) || any(coefficients <= 0)) {
    stop("Invalid status, linear predictor, or baseline coefficients.",
         call. = FALSE)
  }
  cumulative <- as.numeric(basis$itime %*% coefficients) * exp(eta)
  start_cumulative <- as.numeric(basis$start_itime %*% coefficients) * exp(eta)
  entry_cumulative <- as.numeric(basis$entry_itime %*% coefficients) * exp(eta)
  hazard <- as.numeric(basis$time %*% coefficients) * exp(eta)

  # Every contribution is conditioned on survival to the entry time `a` (left
  # truncation), so each is expressed through the ENTRY-CONDITIONED cumulative
  # hazard H(t) - H(a). Conditioning must be applied inside the censoring
  # probabilities, not added afterwards: writing the left-censoring term as
  # (1 - S(t)) / S(a) instead of 1 - S(t)/S(a) is not a probability at all and
  # can exceed one (with H(t) = 3 and H(a) = 2 it returns 7.0).
  cum <- cumulative - entry_cumulative           # H(t) - H(a) >= 0
  start_cum <- start_cumulative - entry_cumulative

  observed <- status == 1L
  left <- status == 2L
  interval <- status == 3L

  out <- -cum                                    # right-censored: S(t)/S(a)
  out[observed] <- -cum[observed] + log(hazard[observed])
  out[left] <- log1p(-exp(-cum[left]))           # 1 - S(t)/S(a)
  out[interval] <- log(exp(-start_cum[interval]) - exp(-cum[interval]))
  out
}

#' Parse exact event and censoring information from a data source
#' @noRd
.cpaic_survival_outcomes <- function(data, time_col, status_col, start_col,
                                     entry_col, source) {
  missing <- setdiff(c(time_col, status_col), names(data))
  if (length(missing)) {
    stop("Exact survival likelihood requires `", source,
         "` event or censoring rows. Missing column(s): ",
         paste(missing, collapse = ", "), ". Aggregate event counts and ",
         "person-time cannot identify the exact individual likelihood.",
         call. = FALSE)
  }
  time <- as.numeric(data[[time_col]])
  raw_status <- data[[status_col]]
  status <- as.integer(raw_status)
  start_time <- if (start_col %in% names(data)) {
    as.numeric(data[[start_col]])
  } else {
    rep(0, nrow(data))
  }
  entry_time <- if (entry_col %in% names(data)) {
    as.numeric(data[[entry_col]])
  } else {
    rep(0, nrow(data))
  }

  if (any(!is.finite(time)) || any(time <= 0)) {
    stop("Survival outcome times in `", source,
         "` must be positive and finite.", call. = FALSE)
  }
  if (any(is.na(raw_status)) || any(raw_status != status) ||
      any(!status %in% 0:3)) {
    stop("Survival status in `", source,
         "` must use 0 right-censored, 1 observed, 2 left-censored, or 3 ",
         "interval-censored.", call. = FALSE)
  }
  if (any(!is.finite(start_time)) || any(start_time < 0) ||
      any(status == 3L & start_time >= time)) {
    stop("Interval-censored rows in `", source,
         "` need a finite start time satisfying 0 <= start < time.",
         call. = FALSE)
  }
  if (any(!is.finite(entry_time)) || any(entry_time < 0) ||
      any(entry_time >= time)) {
    stop("Delayed-entry times in `", source,
         "` must satisfy 0 <= entry < time.", call. = FALSE)
  }
  # An interval-censored subject is known to fail inside (start, time]. That
  # knowledge is impossible if the subject was not under observation until after
  # the interval had already begun, so the entry time cannot exceed the interval
  # start. Without this constraint the likelihood forms (S(start) - S(time)) /
  # S(entry) with entry > start, which is not a probability and can exceed one.
  if (any(status == 3L & entry_time > start_time)) {
    stop("Interval-censored rows in `", source,
         "` need entry <= start: a subject cannot be known to fail after a ",
         "time it was not yet under observation for.", call. = FALSE)
  }

  list(time = time, status = status, start_time = start_time,
       entry_time = entry_time)
}

