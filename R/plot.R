# Plots: network graph and forest --------------------------------------------

#' Plot the component network
#'
#' Draws the treatment network, coloring nodes by sub-network so that
#' disconnection is visible at a glance. Edges with population-adjusted IPD
#' evidence are not distinguished here (see the vignettes).
#'
#' @param x A [cpaic_network()] object.
#' @param ... Passed to [igraph::plot.igraph()].
#' @return The `igraph` object, invisibly.
#' @export
plot.cpaic_network <- function(x, ...) {
  cols <- x$cols
  edges <- unique(data.frame(from = as.character(x$agd[[cols$treat1]]),
                             to = as.character(x$agd[[cols$treat2]]),
                             stringsAsFactors = FALSE))
  g <- igraph::graph_from_data_frame(edges, directed = FALSE,
                                     vertices = data.frame(name = x$treatments))
  memb <- igraph::components(g)$membership
  palette <- grDevices::hcl.colors(max(memb), "Set 2")
  igraph::V(g)$color <- palette[memb]
  args <- list(...)
  if (is.null(args$vertex.label.cex)) args$vertex.label.cex <- 0.8
  if (is.null(args$vertex.size)) args$vertex.size <- 18
  do.call(graphics::plot, c(list(g), args))
  invisible(g)
}

#' Forest plot of relative effects
#'
#' @param x A `cpaic_effects` data frame from [relative_effects()], or a
#'   fit object (then [relative_effects()] is called on it).
#' @param ... Passed to [relative_effects()] when `x` is a fit.
#' @return `x`, invisibly.
#' @export
forest <- function(x, ...) {
  if (!inherits(x, "cpaic_effects")) {
    x <- relative_effects(x, ...)
  }
  bt <- isTRUE(attr(x, "backtransf"))
  null_line <- if (bt) 1 else 0
  est <- x$estimate; lo <- x$lower; hi <- x$upper
  labs <- paste0(x$treatment, " vs ", x$comparator)
  n <- nrow(x)
  yy <- rev(seq_len(n))
  op <- graphics::par(mar = c(4, 12, 2, 2))
  on.exit(graphics::par(op))
  xlim <- range(c(lo, hi, null_line), na.rm = TRUE)
  graphics::plot(NA, xlim = xlim, ylim = c(0.5, n + 0.5), yaxt = "n",
                 xlab = paste0(attr(x, "sm"),
                               if (bt) "" else " (link scale)"),
                 ylab = "")
  graphics::abline(v = null_line, lty = 2, col = "grey50")
  graphics::segments(lo, yy, hi, yy)
  graphics::points(est, yy, pch = 19)
  graphics::axis(2, at = yy, labels = labs, las = 1, cex.axis = 0.8)
  invisible(x)
}
