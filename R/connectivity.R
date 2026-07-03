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

#' Assess connectivity and component-bridge identifiability of a network
#'
#' Reports whether the treatment network is connected, and -- crucially for
#' a disconnected network -- whether the additive component structure makes
#' all component effects estimable. A disconnected network is *bridgeable*
#' only when sub-networks share enough components that the component design
#' matrix `X = B C` has full column rank (`rank(X) = number of
#' components`). Otherwise component effects cannot be uniquely identified
#' and the bridge would produce arbitrary estimates.
#'
#' @param network A [cpaic_network()] object.
#' @param tol Numerical tolerance for the rank computation.
#'
#' @return An object of class `cpaic_connectivity`: a list with
#'   `connected` (logical), `n_subnetworks`, `subnetworks` (list of
#'   treatment-label vectors), `bridging_components` (components shared
#'   across sub-networks), `rank` and `n_components`, `identifiable`
#'   (logical: `rank == n_components`), and the `B`/`C`/`X` matrices.
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

  # Component design matrix X = B C and its rank.
  B <- build_B_matrix(t1, t2, trts)
  X <- B %*% C
  rk <- as.integer(Matrix::rankMatrix(X, tol = tol))
  n_comp <- ncol(C)
  identifiable <- isTRUE(rk == n_comp)

  structure(
    list(
      connected = connected,
      n_subnetworks = comp$no,
      subnetworks = subnetworks,
      bridging_components = bridging,
      rank = rk,
      n_components = n_comp,
      identifiable = identifiable,
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
  }
  cat("  Component identifiability: rank(X) = ", x$rank, " / ",
      x$n_components, " components -> ",
      if (x$identifiable) "IDENTIFIABLE" else "NOT identifiable",
      "\n", sep = "")
  if (!x$identifiable) {
    cat("  ! Component effects are not uniquely estimable. The network\n",
        "    cannot be bridged as specified; check that sub-networks share\n",
        "    enough components.\n", sep = "")
  }
  invisible(x)
}
