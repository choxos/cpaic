# Random effects for component-additive ML-NMR
#
# The contrast-based multi-arm correlation follows the ML-NMR construction in
# multinma (Phillippo et al. 2020). Each study uses one baseline arm. The
# remaining arm deviations have unit variance and pairwise correlation 0.5.

#' Baseline-shift correlation of the study-arm deviations
#'
#' The standard network meta-analysis structure: each non-baseline arm deviation
#' has unit variance, and two deviations from the same study correlate at 0.5.
#' This reproduces `multinma::RE_cor(type = "blshift")` (Phillippo et al. 2020,
#' GPL-3) for the rows this package keeps, which is the whole of what cpaic used
#' that package for at run time. `study` must list the studies of the retained
#' (non-baseline) arms in `re_idx` order; `tests/testthat/test-cmlnmr-random.R`
#' checks the two agree exactly.
#' @noRd
.cpaic_re_cor_blshift <- function(study) {
  study <- as.character(study)
  R <- outer(study, study, "==") * 0.5
  diag(R) <- 1
  R
}

#' Build the random-effects correlation and row mappings
#' @noRd
.cpaic_random_effects <- function(ipd, agd, study_col, trt_col) {
  arm_rows <- unique(rbind(
    data.frame(
      study = as.character(ipd[[study_col]]),
      trt = as.character(ipd[[trt_col]]),
      stringsAsFactors = FALSE
    ),
    data.frame(
      study = as.character(agd[[study_col]]),
      trt = as.character(agd[[trt_col]]),
      stringsAsFactors = FALSE
    )
  ))
  studies <- sort(unique(arm_rows$study))
  arm_rows <- arm_rows[order(match(arm_rows$study, studies), arm_rows$trt), ]

  mapping <- arm_rows
  mapping$re_idx <- 0L
  offset <- 0L
  for (ss in studies) {
    rows <- which(mapping$study == ss)
    m <- length(rows) - 1L
    if (m < 1L) next
    ids <- offset + seq_len(m)
    mapping$re_idx[rows[-1L]] <- ids
    offset <- offset + m
  }

  L_delta <- if (offset > 1L) {
    # The retained rows are exactly the non-baseline arms, already ordered by
    # re_idx, so their studies are all the correlation structure needs.
    keep <- order(mapping$re_idx[mapping$re_idx > 0L])
    corr <- .cpaic_re_cor_blshift(mapping$study[mapping$re_idx > 0L][keep])
    t(chol(corr))
  } else if (offset == 1L) {
    matrix(1, nrow = 1L, ncol = 1L)
  } else {
    matrix(numeric(), nrow = 0L, ncol = 0L)
  }
  keys <- paste(mapping$study, mapping$trt, sep = "\r")
  row_index <- function(data) {
    match(paste(as.character(data[[study_col]]),
                as.character(data[[trt_col]]), sep = "\r"), keys)
  }

  list(
    N_delta = offset,
    L_delta = L_delta,
    re_idx_ipd = as.integer(mapping$re_idx[row_index(ipd)]),
    re_idx_agd = as.integer(mapping$re_idx[row_index(agd)]),
    arms = mapping
  )
}
