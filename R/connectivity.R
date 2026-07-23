# Network connectivity and component-bridge identifiability -------------------

#' Edge-incidence (contrast) matrix B for a set of comparisons
#'
#' One row per comparison, with `+1` in the `treat1` column and `-1` in the
#' `treat2` column. Columns follow `treatments`.
#' @noRd
build_B_matrix <- function(treat1, treat2, treatments) {
  treat1 <- as.character(treat1)
  treat2 <- as.character(treat2)
  B <- matrix(0L, nrow = length(treat1), ncol = length(treatments),
              dimnames = list(NULL, treatments))
  idx1 <- match(treat1, treatments)
  idx2 <- match(treat2, treatments)
  B[cbind(seq_along(treat1), idx1)] <- 1L
  B[cbind(seq_along(treat2), idx2)] <- -1L
  B
}

#' Orthonormal basis of the null space of a design matrix
#'
#' Columns span `null(X)`. A contrast vector `l` (in component space) is
#' uniquely estimable if and only if `l` lies in the row space of `X`,
#' equivalently `l` is orthogonal to every column of this basis.
#' @noRd
.cpaic_null_space <- function(X, tol = 1e-8) {
  X <- as.matrix(X)
  s <- svd(X, nu = 0L, nv = ncol(X))
  d <- if (length(s$d)) s$d else 0
  rk <- sum(d > tol * max(1, max(d)))
  if (rk >= ncol(X)) {
    return(matrix(numeric(0), nrow = ncol(X), ncol = 0L))
  }
  s$v[, seq.int(rk + 1L, ncol(X)), drop = FALSE]
}

#' Is each row of `L` in the row space of the design matrix?
#'
#' Implements the estimability criterion of Wigle et al. (2026): the set of
#' uniquely estimable relative effects of a component NMA is exactly the row
#' space of the design matrix `X = B C`. A contrast is estimable if and only
#' if it is orthogonal to the null space of `X`.
#' @noRd
.cpaic_in_rowspace <- function(L, N, tol = 1e-8) {
  L <- as.matrix(L)
  if (!nrow(L)) return(logical(0))
  if (!ncol(N)) return(rep(TRUE, nrow(L)))          # full column rank
  proj <- L %*% N                                    # component along null(X)
  scale <- pmax(1, sqrt(rowSums(L^2)))
  sqrt(rowSums(proj^2)) / scale < tol
}

#' Estimability table of every treatment versus a reference
#' @noRd
.cpaic_estimable_table <- function(C, N, reference, tol = 1e-8) {
  trts <- rownames(C)
  others <- setdiff(trts, reference)
  L <- C[others, , drop = FALSE] -
    matrix(C[reference, ], nrow = length(others), ncol = ncol(C), byrow = TRUE)
  data.frame(treatment = others, comparator = reference,
             estimable = .cpaic_in_rowspace(L, N, tol = tol),
             row.names = NULL, stringsAsFactors = FALSE)
}

#' Which relative effects of a component network are uniquely estimable?
#'
#' The additive component model identifies a relative effect
#' `theta_i - theta_j = (C_i - C_j)' beta` only when the contrast vector
#' `C_i - C_j` lies in the **row space** of the design matrix `X = B C`
#' (Wigle et al. 2026). Full column rank of `X` (rank equal to the number of
#' components) is *sufficient* for every contrast to be estimable, but it is
#' **not necessary**: a disconnected, rank-deficient component network can
#' still identify many cross-sub-network treatment contrasts.
#'
#' Checking this matters because both engines otherwise return a
#' finite-looking answer for a contrast that carries no information: the
#' frequentist weighted least squares through the Moore-Penrose pseudoinverse,
#' and the Bayesian model through the prior.
#'
#' @param object A [cpaic_network()], [cpaic_connectivity()], `cpaic_bridge`
#'   or `cpaic_mlnmr` object.
#' @param reference Reference treatment. Defaults to the network reference.
#' @param ... Unused.
#'
#' @return A data frame with one row per treatment, giving the `treatment`,
#'   the `comparator` (the reference), and `estimable` (logical).
#' @references
#' Wigle A, Beliveau A, Nikolakopoulou A, Lin L (2026). Creating Treatment and
#' Component Hierarchies in Component Network Meta-Analysis.
#' @examples
#' net <- cpaic_network(cpaic_bin_agd, sm = "OR", inactive = "Placebo")
#' estimable_effects(net)
#' @export
estimable_effects <- function(object, reference = NULL, ...) {
  UseMethod("estimable_effects")
}

#' @export
estimable_effects.cpaic_network <- function(object, reference = NULL, ...) {
  estimable_effects(cpaic_connectivity(object), reference = reference, ...)
}

#' @export
estimable_effects.cpaic_bridge <- function(object, reference = NULL, ...) {
  if (is.null(reference)) reference <- object$reference
  estimable_effects(object$connectivity, reference = reference, ...)
}

#' @export
estimable_effects.cpaic_connectivity <- function(object, reference = NULL,
                                                 ...) {
  C <- object$C
  if (is.null(reference)) reference <- object$reference
  if (!reference %in% rownames(C)) {
    stop("`reference` must be one of the network treatments.", call. = FALSE)
  }
  .cpaic_estimable_table(C, object$null_space, reference)
}

#' @export
estimable_effects.cpaic_mlnmr <- function(object, reference = NULL, ...) {
  if (is.null(reference)) reference <- object$reference
  C <- object$C.matrix
  if (!reference %in% rownames(C)) {
    stop("`reference` must be one of the network treatments.", call. = FALSE)
  }
  .cpaic_estimable_table(C, object$null_space, reference)
}

#' Assess connectivity and component-bridge identifiability of a network
#'
#' Reports whether the treatment network is connected and, for a disconnected
#' network, which relative effects the additive component structure makes
#' estimable.
#'
#' Two distinct questions are answered, and they are not the same (Wigle et
#' al. 2026):
#'
#' * **Are all component effects identified?** Yes if and only if the
#'   component design matrix `X = B C` has full column rank (`rank(X)` equal to
#'   `n_components`). Reported as `identifiable`.
#' * **Is a particular relative effect estimable?** Yes if and only if its
#'   contrast vector lies in the row space of `X`. Full column rank is
#'   sufficient but **not necessary**, so a rank-deficient network can still
#'   identify useful cross-sub-network contrasts. Reported per treatment in
#'   `estimable` (see [estimable_effects()]).
#'
#' @param network A [cpaic_network()] object.
#' @param tol Numerical tolerance for the rank and null-space computations.
#'
#' @return An object of class `cpaic_connectivity`: a list with
#'   `connected` (logical), `n_subnetworks`, `subnetworks` (list of
#'   treatment-label vectors), `bridging_components` (components shared
#'   across sub-networks), `rank` and `n_components`, `identifiable`
#'   (logical: `rank == n_components`), `null_space`, `estimable_components`,
#'   `estimable` (a data frame of estimable relative effects versus the
#'   reference), and the `B`/`C`/`X` matrices.
#' @references
#' Wigle A, Beliveau A, Nikolakopoulou A, Lin L (2026). Creating Treatment and
#' Component Hierarchies in Component Network Meta-Analysis.
#' @seealso [estimable_effects()]
#' @examples
#' net <- cpaic_network(cpaic_bin_agd, sm = "OR", inactive = "Placebo")
#' cpaic_connectivity(net)
#' @export
cpaic_connectivity <- function(network, tol = 1e-8) {
  stopifnot(inherits(network, "cpaic_network"))
  trts <- network$treatments
  C <- network$C.matrix
  cols <- network$cols
  t1 <- network$agd[[cols$treat1]]
  t2 <- network$agd[[cols$treat2]]

  # Treatment graph connectivity via igraph.
  edges <- unique(data.frame(from = as.character(t1), to = as.character(t2),
                             stringsAsFactors = FALSE))
  g <- igraph::graph_from_data_frame(edges, directed = FALSE,
                                     vertices = data.frame(name = trts))
  comp <- igraph::components(g)
  connected <- comp$no == 1L
  membership <- comp$membership[trts]
  subnetworks <- split(trts, membership)
  names(subnetworks) <- paste0("subnetwork_", seq_along(subnetworks))

  # Components present in each sub-network; bridging components appear in
  # treatments belonging to more than one sub-network.
  comp_names <- colnames(C)
  comp_subnet <- lapply(comp_names, function(cc) {
    trts_with <- rownames(C)[C[, cc] == 1L]
    sort(unique(membership[trts_with]))
  })
  names(comp_subnet) <- comp_names
  bridging <- comp_names[vapply(comp_subnet, function(s) length(s) > 1L,
                                logical(1))]

  # Component design matrix X = B C, its rank, and its null space. The null
  # space is what decides estimability of an individual contrast.
  B <- build_B_matrix(t1, t2, trts)
  X <- B %*% C
  rk <- as.integer(Matrix::rankMatrix(X, tol = tol))
  n_comp <- ncol(C)
  identifiable <- isTRUE(rk == n_comp)
  N <- .cpaic_null_space(X, tol = tol)

  est_comp <- stats::setNames(
    .cpaic_in_rowspace(diag(n_comp), N, tol = tol), comp_names)
  est_tbl <- .cpaic_estimable_table(C, N, network$reference, tol = tol)

  structure(
    list(
      connected = connected,
      n_subnetworks = comp$no,
      subnetworks = subnetworks,
      bridging_components = bridging,
      rank = rk,
      n_components = n_comp,
      identifiable = identifiable,
      null_space = N,
      estimable_components = est_comp,
      estimable = est_tbl,
      reference = network$reference,
      B = B, C = C, X = X
    ),
    class = "cpaic_connectivity"
  )
}

#' @export
print.cpaic_connectivity <- function(x, ...) {
  cat("cpaic connectivity\n")
  cat("  Connected network: ", x$connected, "\n", sep = "")
  cat("  Sub-networks:      ", x$n_subnetworks, "\n", sep = "")
  if (!x$connected) {
    for (i in seq_along(x$subnetworks)) {
      cat("    [", i, "] ", length(x$subnetworks[[i]]), " treatments\n",
          sep = "")
    }
    cat("  Bridging components: ",
        if (length(x$bridging_components))
          paste(x$bridging_components, collapse = ", ") else "(none)",
        "\n", sep = "")
    cat("    (components that OCCUR in more than one sub-network; occurrence is\n",
        "     not identifiability, homogeneity, or influence for any contrast.)\n",
        sep = "")
  }
  cat("  Component design:  rank(X) = ", x$rank, " / ", x$n_components,
      " components -> ",
      if (x$identifiable) "all component effects identified"
      else "some component effects NOT identified", "\n", sep = "")
  if (!x$identifiable) {
    bad <- names(x$estimable_components)[!x$estimable_components]
    cat("    Not identified: ", paste(bad, collapse = ", "), "\n", sep = "")
  }
  est <- x$estimable
  n_ok <- sum(est$estimable)
  cat("  Estimable effects: ", n_ok, " / ", nrow(est), " vs ", x$reference,
      "\n", sep = "")
  if (n_ok < nrow(est)) {
    cat("    Not estimable: ",
        paste(est$treatment[!est$estimable], collapse = ", "), "\n", sep = "")
    cat("  ! Rank deficiency alone does not rule out every contrast: only the\n",
        "    listed effects lie outside the row space of X. See\n",
        "    estimable_effects().\n", sep = "")
  }
  invisible(x)
}
