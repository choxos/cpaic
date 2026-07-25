#!/usr/bin/env Rscript
# Precompile cpaic's Stan-fitting vignettes (the multinma / mlumr pattern). Run
# locally from the package root:
#
#   Rscript vignettes/precompile.R
#
# For each "<stem>.Rmd.orig" this:
#   1. knits it to a static "<stem>.Rmd", fitting the Stan models ONCE, here;
#   2. renders a self-contained "<stem>.html" (output: rmarkdown::html_vignette);
#   3. writes "<stem>.html.asis" so R CMD build and CRAN register and serve the
#      pre-rendered HTML through the R.rsp::asis engine and never run Stan.
#
# Re-run whenever a "<stem>.Rmd.orig" changes. This script, the *.Rmd.orig
# sources, the knitted *.Rmd intermediates and figure/ are all build-ignored (see
# .Rbuildignore); only *.html and *.html.asis ship.
#
# Requires: rstan (the default backend, whose models are compiled into the
# installed package), ggplot2 (figures), bayesplot, knitr, rmarkdown, R.rsp.
# The vignettes fit through cmlnmr()'s default backend, so they exercise what a
# reader gets from install.packages() rather than a toolchain most readers do
# not have.

for (pkg in c("knitr", "rmarkdown", "rstan", "ggplot2", "bayesplot")) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("precompile.R needs the '", pkg, "' package installed.", call. = FALSE)
  }
}

# rstan takes its core count from `mc.cores`, which is 1 unless it is set, so a
# default run puts every chain of every fit end to end and the whole pass takes
# roughly four times as long as it needs to. Set it here rather than in the
# vignettes: this is a local build script, and a vignette that sets a global
# option would be telling readers to do the same.
if (is.null(getOption("mc.cores"))) {
  options(mc.cores = max(1L, min(4L, parallel::detectCores() - 1L)))
}
message("Fitting with mc.cores = ", getOption("mc.cores"), ".")

# The vignettes call library(cpaic), so they run against the INSTALLED package,
# not this source tree. If the installed build is stale the vignettes are built
# from the wrong code, and the failure is not always loud: a plot method that is
# not yet registered falls through to the base-R generic and produces a WRONG
# FIGURE with no error at all. Refuse to run unless the installed package
# matches the source.
.assert_install_current <- function() {
  if (!requireNamespace("cpaic", quietly = TRUE)) {
    stop("cpaic is not installed. Run R CMD INSTALL . from the package root ",
         "before precompiling.", call. = FALSE)
  }
  root <- if (basename(getwd()) == "vignettes") ".." else "."
  stale <- character()
  for (f in list.files(file.path(root, "inst", "stan"), pattern = "[.]stan$")) {
    src <- readLines(file.path(root, "inst", "stan", f), warn = FALSE)
    dst <- system.file("stan", f, package = "cpaic")
    if (!nzchar(dst) || !identical(src, readLines(dst, warn = FALSE))) {
      stale <- c(stale, f)
    }
  }
  ns <- getNamespaceExports("cpaic")
  missing_fns <- setdiff(c("plot_estimability", "plot_rank_curve", "rank_probs",
                           "plot_edge_influence", "plot_prior_posterior"), ns)
  unreg <- !any(grepl("^plot",
                      as.character(utils::.S3methods(class = "cpaic_ranks"))))
  if (length(stale) || length(missing_fns) || unreg) {
    stop("The installed cpaic is STALE, so the vignettes would be built from ",
         "the wrong code.\n",
         if (length(stale)) paste0("  Stan models differing from source: ",
                                   paste(stale, collapse = ", "), "\n"),
         if (length(missing_fns)) paste0("  Functions not exported: ",
                                         paste(missing_fns, collapse = ", "),
                                         "\n"),
         if (unreg) "  plot.cpaic_ranks is not registered, so plot() would ",
         if (unreg) "silently draw the WRONG figure.\n",
         "Run R CMD INSTALL . from the package root, then precompile again.",
         call. = FALSE)
  }
  message("Installed cpaic matches the source. Proceeding.")
}
.assert_install_current()

all_stems <- c(
  "binary-outcomes",
  "continuous-outcomes",
  "count-outcomes",
  "survival-outcomes",
  "cpaic-disconnected-myeloma"
)

# Naming stems on the command line re-renders only those, which is what you want
# after editing one vignette: a full pass fits every Stan model in the package
# and takes hours.
#
#   Rscript vignettes/precompile.R survival-outcomes
stems <- commandArgs(trailingOnly = TRUE)
if (!length(stems)) {
  stems <- all_stems
} else if (length(unknown <- setdiff(stems, all_stems))) {
  stop("No such vignette stem: ", paste(unknown, collapse = ", "),
       ". Known stems: ", paste(all_stems, collapse = ", "), call. = FALSE)
}

# Operate inside vignettes/ so figure paths and the bibliography path resolve
# exactly as they will for pkgdown.
if (basename(getwd()) != "vignettes") setwd("vignettes")

precompile_one <- function(stem) {
  orig <- paste0(stem, ".Rmd.orig")
  rmd  <- paste0(stem, ".Rmd")
  message("\n=== precompiling ", orig, " ===")
  knitr::knit(orig, output = rmd)            # runs the chunks -> fits Stan here
  rmarkdown::render(rmd, quiet = TRUE)       # output format taken from the YAML
  title <- rmarkdown::yaml_front_matter(orig)$title
  writeLines(
    c(sprintf("%%\\VignetteIndexEntry{%s}", title),
      "%\\VignetteEngine{R.rsp::asis}",
      "%\\VignetteEncoding{UTF-8}"),
    paste0(stem, ".html.asis")
  )
  message("    -> ", rmd, ", ", stem, ".html, ", stem, ".html.asis")
}

for (s in stems) precompile_one(s)
message("\nAll vignettes precompiled.")
