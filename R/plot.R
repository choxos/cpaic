# Plots: network, forest, and edge influence -----------------------------------
#
# The plotting surface of cpaic is modeled on that of multinma (Phillippo et
# al. 2020, <doi:10.1111/rssa.12579>): network plots, forest plots of relative
# effects, rankograms, deviance and leverage plots, prior-versus-posterior
# overlays, integration-error plots, and survival curves. The plot LOGIC is
# ported and re-implemented here on top of ggplot2 alone; multinma's public
# exports are used where useful, and no `multinma:::` internal is touched. Both
# packages are licensed under GPL-3.
#
# Plots specific to cpaic, with no counterpart in multinma, live in
# R/plot-mlnmr.R: the population-dependent rank curve, the estimability map, and
# the edge-influence plot.

# ggplot2 is a Suggests dependency; every plot guards on it. Bare column names
# used inside aes() are declared here so R CMD check sees them bound.
utils::globalVariables(c(
  ".angle", ".label", ".nx", ".ny", ".x", ".xend", ".y", ".yend",
  "bridges", "contrast", "deviance_x", "deviance_y", "edge_type", "element",
  "estimate", "identified_by", "influence", "leverage", "lower", "n_int",
  "n_studies", "probability", "rank_position", "ssrd", "subnetwork", "surv",
  "target", "time", "type", "upper", "value", "zero_influence"
))

#' Is ggplot2 available?
#' @noRd
.cpaic_need_ggplot <- function(what = "This plot") {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop(what, " needs the 'ggplot2' package. Install it with ",
         "install.packages(\"ggplot2\").", call. = FALSE)
  }
  invisible(TRUE)
}

#' A light ggplot2 theme shared by the cpaic plots
#'
#' Adapted from `multinma::theme_multinma()` (GPL-3).
#' @noRd
.cpaic_theme <- function() {
  ggplot2::theme_light() +
    ggplot2::theme(
      panel.border = ggplot2::element_rect(colour = "grey70", fill = NA),
      panel.grid.major = ggplot2::element_line(colour = "grey95"),
      panel.grid.minor = ggplot2::element_line(colour = "grey95"),
      strip.background = ggplot2::element_rect(colour = "grey70",
                                               fill = "grey90"),
      strip.text = ggplot2::element_text(colour = "black")
    )
}

#' Axis label for a summary measure on its reporting scale
#' @noRd
.cpaic_sm_label <- function(sm, backtransf) {
  if (is.null(sm) || !nzchar(sm)) sm <- "Effect"
  if (isTRUE(backtransf)) sm else paste0(sm, " (link scale)")
}

# Network plot ----------------------------------------------------------------

#' Layout a (possibly disconnected) treatment network
#'
#' Each sub-network is laid out on its own circle and the circles are placed
#' side by side, so a disconnected network *looks* disconnected. Deterministic,
#' unlike a force-directed layout.
#' @noRd
.cpaic_net_layout <- function(membership, treatments) {
  memb <- membership[treatments]
  ks <- sort(unique(memb))
  radius <- vapply(ks, function(k) max(1, sum(memb == k) / (2 * pi) * 1.6),
                   numeric(1))
  # Centre each circle so that neighboring circles do not overlap.
  gap <- 1.2
  centres <- numeric(length(ks))
  running <- 0
  for (i in seq_along(ks)) {
    centres[i] <- running + radius[i]
    running <- running + 2 * radius[i] + gap
  }
  out <- data.frame(name = treatments, .nx = NA_real_, .ny = NA_real_,
                    subnetwork = NA_character_, stringsAsFactors = FALSE)
  for (i in seq_along(ks)) {
    idx <- which(memb == ks[i])
    n <- length(idx)
    ang <- if (n == 1L) 0 else seq(pi / 2, by = 2 * pi / n, length.out = n)
    out$.nx[idx] <- centres[i] + if (n == 1L) 0 else radius[i] * cos(ang)
    out$.ny[idx] <- if (n == 1L) 0 else radius[i] * sin(ang)
    out$subnetwork[idx] <- paste("Sub-network", ks[i])
  }
  # Label direction: outward from the centre of the node's own circle.
  out$.angle <- atan2(out$.ny, out$.nx - centres[match(memb, ks)])
  out
}

#' Plot the component network
#'
#' Draws the treatment network with the disconnection made explicit: nodes are
#' colored by sub-network, edges are colored by whether the study carries
#' individual patient data (IPD) or aggregate data (AgD) only, and the
#' components that bridge the sub-networks are named in the subtitle. Treatments
#' that contain a bridging component are outlined, because those are the nodes
#' through which the additive component model reconnects the network.
#'
#' Each sub-network is laid out on its own circle, so a disconnected network
#' looks disconnected. Ported in spirit from `multinma::plot.nma_data()`
#' (Phillippo et al. 2020), re-implemented on ggplot2 without a `ggraph`
#' dependency.
#'
#' @param x A [cpaic_network()] object.
#' @param ... Unused.
#' @param weight_edges Scale edge width by the number of studies contributing
#'   to a comparison? Default `TRUE`.
#' @param show_bridges Outline treatments that contain a bridging component,
#'   and name the bridging components in the subtitle? Default `TRUE`.
#' @param nudge Distance by which treatment labels are pushed away from their
#'   node. Default `0.25`.
#'
#' @return A `ggplot` object, so it can be modified with the usual ggplot2
#'   verbs.
#' @seealso [cpaic_connectivity()], [forest()]
#' @examplesIf requireNamespace("ggplot2", quietly = TRUE)
#' net <- cpaic_network(cpaic_bin_agd, ipd = cpaic_bin_ipd, sm = "OR",
#'                      family = "binomial", ipd_covariates = "x1",
#'                      inactive = "Placebo")
#' plot(net)
#' @export
plot.cpaic_network <- function(x, ..., weight_edges = TRUE,
                               show_bridges = TRUE, nudge = 0.25) {
  .cpaic_need_ggplot("plot.cpaic_network()")
  stopifnot(inherits(x, "cpaic_network"))
  if (!is.logical(weight_edges) || length(weight_edges) != 1L ||
      is.na(weight_edges)) {
    stop("`weight_edges` must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.logical(show_bridges) || length(show_bridges) != 1L ||
      is.na(show_bridges)) {
    stop("`show_bridges` must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.numeric(nudge) || length(nudge) != 1L || !is.finite(nudge)) {
    stop("`nudge` must be a single finite number.", call. = FALSE)
  }

  conn <- cpaic_connectivity(x)
  cols <- x$cols
  trts <- x$treatments
  ipd_studies <- if (is.null(x$ipd_info)) character(0) else x$ipd_info$studies

  edges <- data.frame(
    from = as.character(x$agd[[cols$treat1]]),
    to = as.character(x$agd[[cols$treat2]]),
    studlab = as.character(x$agd[[cols$studlab]]),
    stringsAsFactors = FALSE)
  edges$is_ipd <- edges$studlab %in% ipd_studies
  key <- paste(pmin(edges$from, edges$to), pmax(edges$from, edges$to),
               sep = "\r")
  agg <- data.frame(
    key = unique(key),
    stringsAsFactors = FALSE)
  agg$n_studies <- vapply(agg$key,
                          function(k) length(unique(edges$studlab[key == k])),
                          numeric(1))
  agg$is_ipd <- vapply(agg$key, function(k) any(edges$is_ipd[key == k]),
                       logical(1))
  agg$from <- vapply(agg$key, function(k) edges$from[key == k][1], character(1))
  agg$to <- vapply(agg$key, function(k) edges$to[key == k][1], character(1))
  agg$edge_type <- ifelse(agg$is_ipd, "IPD", "AgD")

  # Sub-network membership from the treatment graph (via igraph, an Import).
  g <- igraph::graph_from_data_frame(
    unique(edges[, c("from", "to")]), directed = FALSE,
    vertices = data.frame(name = trts))
  memb <- igraph::components(g)$membership

  nodes <- .cpaic_net_layout(memb, trts)
  bridges <- conn$bridging_components
  C <- x$C.matrix
  nodes$bridges <- if (length(bridges)) {
    rowSums(C[nodes$name, bridges, drop = FALSE]) > 0
  } else {
    rep(FALSE, nrow(nodes))
  }
  pos <- match(agg$from, nodes$name)
  pos2 <- match(agg$to, nodes$name)
  agg$.x <- nodes$.nx[pos]; agg$.y <- nodes$.ny[pos]
  agg$.xend <- nodes$.nx[pos2]; agg$.yend <- nodes$.ny[pos2]

  nodes$.label <- nodes$name
  nodes$.x <- nodes$.nx + nudge * cos(nodes$.angle)
  nodes$.y <- nodes$.ny + nudge * sin(nodes$.angle)

  p <- ggplot2::ggplot()
  if (weight_edges) {
    p <- p + ggplot2::geom_segment(
      data = agg,
      ggplot2::aes(x = .x, y = .y, xend = .xend, yend = .yend,
                   colour = edge_type, linewidth = n_studies),
      lineend = "round") +
      ggplot2::scale_linewidth_continuous(
        "Number of studies", range = c(0.4, 2.2),
        breaks = function(z) unique(round(pretty(z))))
  } else {
    p <- p + ggplot2::geom_segment(
      data = agg,
      ggplot2::aes(x = .x, y = .y, xend = .xend, yend = .yend,
                   colour = edge_type),
      linewidth = 0.8, lineend = "round")
  }
  p <- p + ggplot2::scale_colour_manual(
    "Evidence", values = c(AgD = "#113259", IPD = "#55A480"),
    drop = FALSE)

  if (show_bridges && length(bridges)) {
    p <- p + ggplot2::geom_point(
      data = nodes,
      ggplot2::aes(x = .nx, y = .ny, fill = subnetwork, shape = bridges),
      size = 5, colour = "grey20", stroke = 1) +
      ggplot2::scale_shape_manual(
        "Contains a bridging component",
        values = c(`FALSE` = 21, `TRUE` = 24))
  } else {
    p <- p + ggplot2::geom_point(
      data = nodes, ggplot2::aes(x = .nx, y = .ny, fill = subnetwork),
      size = 5, shape = 21, colour = "grey20", stroke = 0.6)
  }

  subtitle <- if (conn$connected) {
    "Connected network"
  } else if (length(bridges)) {
    paste0(conn$n_subnetworks, " sub-networks; bridging components: ",
           paste(bridges, collapse = ", "))
  } else {
    paste0(conn$n_subnetworks,
           " sub-networks; NO bridging component, the gap cannot be closed")
  }

  subnets <- sort(unique(nodes$subnetwork))
  pal <- stats::setNames(
    grDevices::hcl.colors(max(length(subnets), 2L), "Set 2")[
      seq_along(subnets)], subnets)

  p +
    ggplot2::geom_text(
      data = nodes,
      ggplot2::aes(x = .x, y = .y, label = .label),
      hjust = "outward", vjust = "outward", size = 3.2) +
    ggplot2::scale_fill_manual("Sub-network", values = pal) +
    ggplot2::labs(subtitle = subtitle) +
    ggplot2::coord_equal(clip = "off") +
    .cpaic_theme() +
    ggplot2::theme(
      axis.title = ggplot2::element_blank(),
      axis.text = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank()) +
    ggplot2::guides(
      fill = ggplot2::guide_legend(override.aes = list(shape = 21)))
}

# Forest plots ----------------------------------------------------------------

#' Forest plot of relative or component effects
#'
#' A ggplot2 forest plot of the estimates in a cpaic fit. Ported from
#' `multinma::plot.nma_summary()` (Phillippo et al. 2020) and re-implemented on
#' ggplot2 alone.
#'
#' Contrasts that the component design cannot identify are **shown**, labelled
#' `not estimable`, rather than silently dropped. Dropping them would leave the
#' reader with a plot that looks complete when it is not; see
#' [estimable_effects()] and [estimable_effects_at()].
#'
#' @param x A `cpaic_effects` data frame (from [relative_effects()]), a fitted
#'   cpaic object (`cpaic_bridge`, `cpaic_maic`, `cpaic_stc`, `cpaic_mlnmr`),
#'   or a component-effect data frame from [component_effects()].
#' @param ... Passed to [relative_effects()] / [component_effects()] when `x` is
#'   a fit (for example `newdata` for a [cmlnmr()] fit).
#' @param what `"relative"` (default) for relative effects, or `"component"` for
#'   the incremental effect of each component.
#' @param order Row ordering: `"estimate"` (default, most to least favorable),
#'   `"alphabetical"`, or `"none"` (the order in the input).
#' @param ref_line Position of the vertical reference line. Defaults to the null
#'   value of the summary measure (`1` on a back-transformed ratio scale, `0`
#'   otherwise); `NA` draws none.
#' @param fatten Relative size of the point estimate marker. Default `2`.
#' @param show_na Show non-estimable contrasts as labelled empty rows? Default
#'   `TRUE`. Setting this to `FALSE` hides evidence that the network cannot
#'   answer part of your question, so leave it on unless you have a reason.
#'
#' @return A `ggplot` object.
#' @seealso [relative_effects()], [component_effects()], [league_table()]
#' @examplesIf requireNamespace("ggplot2", quietly = TRUE)
#' net <- cpaic_network(cpaic_bin_agd, sm = "OR", inactive = "Placebo")
#' br <- cnma_bridge(net)
#' forest(br)
#' forest(br, what = "component")
#' @export
forest <- function(x, ..., what = c("relative", "component"),
                   order = c("estimate", "alphabetical", "none"),
                   ref_line = NULL, fatten = 2, show_na = TRUE) {
  .cpaic_need_ggplot("forest()")
  what <- match.arg(what)
  order <- match.arg(order)

  is_effects <- inherits(x, "cpaic_effects")
  is_comp_df <- is.data.frame(x) && !is_effects &&
    all(c("component", "estimate", "lower", "upper") %in% names(x))

  if (is_effects) {
    tab <- x
    what <- "relative"
  } else if (is_comp_df) {
    tab <- x
    what <- "component"
  } else if (what == "component") {
    tab <- component_effects(x, ...)
  } else {
    tab <- relative_effects(x, ...)
  }
  sm <- attr(tab, "sm")
  target <- attr(tab, "target")
  # Component effects are always reported on the link (log) scale, matching
  # print.cpaic_bridge(); only relative_effects() back-transforms.
  backtransf <- what == "relative" && isTRUE(attr(tab, "backtransf"))
  tab <- as.data.frame(tab)

  if (what == "component") {
    tab$.label <- as.character(tab$component)
    if (is.null(sm) && inherits(x, "cpaic_fit")) sm <- x$sm
    ylab <- "Component"
  } else {
    tab$.label <- paste0(tab$treatment, " vs ", tab$comparator)
    ylab <- "Contrast"
  }

  tab$estimable <- is.finite(tab$estimate) & is.finite(tab$lower) &
    is.finite(tab$upper)
  if (!show_na) tab <- tab[tab$estimable, , drop = FALSE]
  if (!nrow(tab)) {
    stop("Nothing to plot: no estimate is available.", call. = FALSE)
  }

  if (is.null(ref_line)) ref_line <- if (backtransf) 1 else 0
  if (length(ref_line) != 1L ||
      !(is.numeric(ref_line) || (is.logical(ref_line) && is.na(ref_line)))) {
    stop("`ref_line` must be a single number, or NA for no reference line.",
         call. = FALSE)
  }
  ref_line <- as.numeric(ref_line)

  ord <- switch(
    order,
    # Non-estimable rows sink to the bottom of the plot.
    estimate = order(!tab$estimable, tab$estimate, na.last = TRUE),
    alphabetical = order(!tab$estimable, tab$.label),
    none = seq_len(nrow(tab))
  )
  tab <- tab[ord, , drop = FALSE]
  # In ggplot2 the first factor level is drawn at the bottom, so reverse.
  tab$.label <- factor(tab$.label, levels = rev(unique(tab$.label)))

  ok <- tab[tab$estimable, , drop = FALSE]
  bad <- tab[!tab$estimable, , drop = FALSE]

  xr <- range(c(ok$lower, ok$upper, ref_line), na.rm = TRUE)
  if (!all(is.finite(xr))) xr <- if (backtransf) c(0.5, 2) else c(-1, 1)
  if (diff(xr) == 0) xr <- xr + c(-0.5, 0.5) * max(abs(xr), 1)
  # Park the "not estimable" labels in the middle of the plotted range.
  if (nrow(bad)) {
    bad$estimate <- if (backtransf) exp(mean(log(xr))) else mean(xr)
  }

  # The full table is the plot's default data, so a caller can add layers with
  # ggplot2::aes() and see the non-estimable rows too.
  p <- ggplot2::ggplot(
    tab, ggplot2::aes(x = estimate, y = .label, xmin = lower, xmax = upper))
  if (is.finite(ref_line)) {
    p <- p + ggplot2::geom_vline(xintercept = ref_line, colour = "grey60",
                                 linetype = 2)
  }
  if (nrow(ok)) {
    p <- p + ggplot2::geom_pointrange(
      data = ok, orientation = "y", fatten = fatten)
  }
  if (nrow(bad)) {
    p <- p + ggplot2::geom_text(
      data = bad, ggplot2::aes(x = estimate, y = .label),
      label = "not estimable", colour = "grey45", fontface = "italic",
      size = 3, inherit.aes = FALSE)
  }

  caption <- if (nrow(bad)) {
    paste0("Not estimable from this component design: ",
           paste(as.character(bad$.label), collapse = "; "),
           ". See estimable_effects().")
  } else {
    NULL
  }
  subtitle <- if (!is.null(target) && length(target)) {
    paste0("Target population: ",
           paste(names(target), signif(target, 3), sep = " = ",
                 collapse = ", "))
  } else {
    NULL
  }

  p <- p +
    ggplot2::labs(x = .cpaic_sm_label(sm, backtransf), y = ylab,
                  subtitle = subtitle, caption = caption) +
    .cpaic_theme()
  # A ratio measure is symmetric on the log scale, so plot it there.
  if (backtransf) p <- p + ggplot2::scale_x_log10()
  p
}

#' @rdname forest
#' @param y Unused, for compatibility with the [plot()] generic.
#' @export
plot.cpaic_effects <- function(x, y, ...) forest(x, ...)

#' @rdname forest
#' @export
plot.cpaic_bridge <- function(x, y, ...) forest(x, ...)

#' @rdname forest
#' @export
plot.cpaic_fit <- function(x, y, ...) forest(x, ...)

# Edge influence ---------------------------------------------------------------

#' Plot how much each edge informs a chosen contrast
#'
#' Visualizes [edge_influence()]: the weight with which each observed edge
#' enters the estimate of one relative effect. This plot exists because the
#' usual population-adjustment diagnostics cannot see the failure it detects.
#' An IPD edge with **zero influence** on your contrast contributes nothing to
#' it, so reweighting that edge cannot move the answer, however healthy its
#' effective sample size looks. Such edges are drawn in red and labelled.
#'
#' There is no counterpart to this plot in multinma; it is specific to the
#' component bridge, where a contrast is a weighted combination of edges chosen
#' by the component design rather than by the network path.
#'
#' @param object A `cpaic_bridge`, `cpaic_maic`, or `cpaic_stc` object.
#' @param treatment,comparator The contrast of interest. `comparator` defaults
#'   to the network reference.
#' @param ... Passed to [edge_influence()] (for example `tol`).
#'
#' @return A `ggplot` object.
#' @seealso [edge_influence()], [effective_sample_size()]
#' @examplesIf requireNamespace("ggplot2", quietly = TRUE)
#' net <- cpaic_network(cpaic_bin_agd, ipd = cpaic_bin_ipd, sm = "OR",
#'                      family = "binomial", ipd_covariates = "x1",
#'                      inactive = "Placebo")
#' br <- cnma_bridge(net)
#' plot_edge_influence(br, treatment = "A+B+C")
#' @export
plot_edge_influence <- function(object, treatment, comparator = NULL, ...) {
  .cpaic_need_ggplot("plot_edge_influence()")
  infl <- suppressWarnings(
    edge_influence(object, treatment = treatment, comparator = comparator, ...))
  if (is.null(comparator)) {
    comparator <- if (inherits(object, "cpaic_fit")) object$bridge$reference
                  else object$reference
  }

  infl$.label <- paste0(infl$studlab, ": ", infl$treat1, " vs ", infl$treat2)
  scale <- max(abs(infl$influence), 1)
  infl$zero_influence <- infl$has_ipd & abs(infl$influence) < 1e-8 * scale
  infl$edge_type <- ifelse(infl$zero_influence, "IPD, no influence",
                           ifelse(infl$has_ipd, "IPD", "AgD"))
  infl$edge_type <- factor(infl$edge_type,
                           levels = c("AgD", "IPD", "IPD, no influence"))
  infl <- infl[order(infl$influence), , drop = FALSE]
  infl$.label <- factor(infl$.label, levels = unique(infl$.label))

  dead <- infl$.label[infl$zero_influence]
  caption <- if (length(dead)) {
    paste0("Zero influence on this contrast, so adjusting it changes nothing: ",
           paste(as.character(dead), collapse = "; "), ".")
  } else {
    NULL
  }

  ggplot2::ggplot(infl,
                  ggplot2::aes(x = influence, y = .label, fill = edge_type)) +
    ggplot2::geom_vline(xintercept = 0, colour = "grey60") +
    ggplot2::geom_col(width = 0.65, colour = "grey30", linewidth = 0.2) +
    ggplot2::scale_fill_manual(
      "Edge",
      values = c(AgD = "#113259", IPD = "#55A480",
                 `IPD, no influence` = "#B2182B"),
      drop = FALSE) +
    ggplot2::labs(
      x = "Influence weight on the contrast",
      y = "Edge",
      subtitle = paste0("Contribution of each edge to ", treatment, " vs ",
                        comparator),
      caption = caption) +
    .cpaic_theme()
}
