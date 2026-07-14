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
# Requires (Suggests): cmdstanr (fits the models), ggplot2 (figures), knitr,
# rmarkdown, R.rsp.

for (pkg in c("knitr", "rmarkdown", "cmdstanr", "ggplot2", "bayesplot")) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("precompile.R needs the '", pkg, "' package installed.", call. = FALSE)
  }
}
if (inherits(try(cmdstanr::cmdstan_path(), silent = TRUE), "try-error")) {
  stop("precompile.R needs a working CmdStan installation.", call. = FALSE)
}

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

stems <- c(
  "binary-outcomes",
  "continuous-outcomes",
  "count-outcomes",
  "survival-outcomes",
  "cpaic-disconnected-myeloma"
)

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
