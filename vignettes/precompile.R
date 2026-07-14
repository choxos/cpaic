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

for (pkg in c("knitr", "rmarkdown", "cmdstanr", "ggplot2")) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("precompile.R needs the '", pkg, "' package installed.", call. = FALSE)
  }
}
if (inherits(try(cmdstanr::cmdstan_path(), silent = TRUE), "try-error")) {
  stop("precompile.R needs a working CmdStan installation.", call. = FALSE)
}

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
