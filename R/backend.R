# Sampler backends: rstan (default) and cmdstanr -------------------------------
#
# cpaic fits the same five Stan models through either engine. rstan is the
# default because it is on CRAN and its models are compiled when the package is
# installed, so `cmlnmr()` works out of the box and CRAN can actually run the
# examples and tests. cmdstanr tracks Stan releases more closely and is often
# faster, but it is not on CRAN and needs a separate CmdStan toolchain, so it is
# opt-in.
#
# Everything downstream of the fit goes through the accessors below rather than
# touching an engine's API directly. They dispatch on the fit's class, not on a
# recorded flag, so a saved object still resolves correctly.

#' Is a sampler backend usable right now?
#'
#' rstan needs only the package. cmdstanr needs the package AND a working
#' CmdStan installation, which is a separate download.
#' @noRd
.cpaic_backend_ready <- function(backend) {
  if (identical(backend, "rstan")) {
    return(requireNamespace("rstan", quietly = TRUE))
  }
  requireNamespace("cmdstanr", quietly = TRUE) &&
    !is.null(tryCatch(cmdstanr::cmdstan_path(), error = function(e) NULL))
}

#' Explain why a backend cannot be used
#' @noRd
.cpaic_backend_stop <- function(backend) {
  if (identical(backend, "rstan")) {
    stop("backend = \"rstan\" needs the 'rstan' package. Install it with ",
         "install.packages(\"rstan\").", call. = FALSE)
  }
  stop("backend = \"cmdstanr\" needs the 'cmdstanr' package and a working ",
       "CmdStan installation. Install the package from ",
       "https://stan-dev.r-universe.dev , then run ",
       "cmdstanr::install_cmdstan(). Or use backend = \"rstan\", which needs ",
       "no separate toolchain.", call. = FALSE)
}

#' What `rstan::sampling()` will accept from a caller
#'
#' `...` reaches the sampler unchanged, so anything cpaic lets through must be
#' something rstan understands. cmdstanr's `$sample()` takes about twenty
#' arguments rstan has no equivalent for (`step_size`, `metric`, `inv_metric`,
#' `adapt_engaged`, `init_buffer`, `fixed_param`, `save_latent_dynamics`, ...),
#' and handing one to `rstan::sampling()` kills every chain with "no applicable
#' method for `@`", which says nothing about the cause.
#'
#' This is a whitelist rather than a list of known-bad names, because the
#' known-bad list can only ever be a guess at cmdstanr's surface, and it also
#' silently lets a misspelling (`adapt_dleta`) through to be ignored. The names
#' are the formals of the `.local()` inside rstan's `sampling()` method for
#' `stanmodel`, plus the dot-arguments `?sampling` documents, minus the ones
#' cpaic derives from its own arguments (`object`, `data`, `chains`, `iter`,
#' `warmup`, `seed`). `control` and `refresh` are here because cpaic supplies a
#' default for each and then steps aside if the caller sets one.
#' @noRd
.cpaic_rstan_sample_args <- c(
  "pars", "thin", "init", "check_data", "sample_file", "diagnostic_file",
  "verbose", "algorithm", "control", "include", "cores", "open_progress",
  "show_messages", "refresh",
  "chain_id", "init_r", "test_grad", "append_samples")

#' rstan arguments cmlnmr() computes for itself
#'
#' These are real `rstan::sampling()` arguments, so they need a different
#' message from a name rstan has never heard of: the caller is not wrong about
#' rstan, only about where to set it. `iter` is the trap worth naming, since
#' rstan counts warmup inside it and cpaic does not.
#' @noRd
.cpaic_rstan_derived_args <- c("object", "data", "chains", "iter", "warmup",
                               "seed")

#' Draw from the posterior with the requested backend
#'
#' `iter_warmup` and `iter_sampling` are cpaic's (and cmdstanr's) names. rstan
#' counts `iter` as the TOTAL including warmup, so the translation happens here
#' rather than being left to the caller.
#'
#' `adapt_delta` and `max_treedepth` are named arguments because the two engines
#' take them in different places: top level for cmdstanr, inside `control` for
#' rstan. Anything else in `sample_args` is passed through untouched and is
#' therefore backend-specific.
#'
#' Both branches drop their own defaults for anything the caller supplied, so
#' `refresh` (and cmdstanr's `show_messages`) behave the same way on either
#' engine. Building the argument list by concatenation instead would hand
#' `do.call()` two elements of the same name and fail with "formal argument
#' matched by multiple actual arguments".
#'
#' The rstan branch deliberately does not set `cores`. rstan's own default reads
#' `getOption("mc.cores")`, which is the R convention and is what CRAN's
#' two-core limit is enforced through; hard-coding `cores = chains` here would
#' override a user's setting and break that limit during checks.
#' @noRd
.cpaic_sample <- function(backend, family, standata, chains, iter_warmup,
                          iter_sampling, seed, adapt_delta = NULL,
                          max_treedepth = NULL, sample_args = list()) {
  if (identical(backend, "rstan")) {
    # `iter`, `warmup` and friends ARE rstan arguments; cpaic just computes
    # them, so saying they are not would be false and would send the caller
    # looking in the wrong place.
    derived <- intersect(names(sample_args), .cpaic_rstan_derived_args)
    if (length(derived)) {
      stop("`", paste(derived, collapse = "`, `"), "` ",
           if (length(derived) == 1L) "is" else "are",
           " derived by cmlnmr() and cannot be set through `...`. Use ",
           "`chains`, `iter_warmup`, `iter_sampling`, and `seed`, which work ",
           "on both backends; rstan's single `iter` counts warmup, cpaic's ",
           "does not.", call. = FALSE)
    }
    clash <- setdiff(names(sample_args), .cpaic_rstan_sample_args)
    if (length(clash)) {
      stop("`", paste(clash, collapse = "`, `"), "` ",
           if (length(clash) == 1L) "is not an argument of rstan::sampling()"
           else "are not arguments of rstan::sampling()",
           ", so the rstan backend cannot take ",
           if (length(clash) == 1L) "it" else "them",
           ". Several cmdstanr sampler arguments (`step_size`, `metric`, ",
           "`adapt_engaged`, `parallel_chains`, ...) have no rstan equivalent; ",
           "drop ", if (length(clash) == 1L) "it" else "them",
           " or fit with backend = \"cmdstanr\". Note that `adapt_delta` and ",
           "`max_treedepth` are arguments of cmlnmr() itself and work on both.",
           call. = FALSE)
    }
    control <- sample_args$control %||% list()
    if (!is.null(adapt_delta)) control$adapt_delta <- adapt_delta
    if (!is.null(max_treedepth)) control$max_treedepth <- max_treedepth
    sample_args$control <- NULL
    defaults <- list(
      object = .cpaic_rstan_model(family), data = standata,
      chains = chains, warmup = iter_warmup,
      iter = iter_warmup + iter_sampling, seed = seed, refresh = 0)
    if (length(control)) defaults$control <- control
    if (length(sample_args)) defaults[names(sample_args)] <- NULL
    return(do.call(rstan::sampling, c(defaults, sample_args)))
  }

  # .cpaic_stan_model() prepends "cpaic_" itself, so it takes the bare family
  # stem. Handing it the already-prefixed name looks for cpaic_cpaic_*.stan.
  mod <- .cpaic_stan_model(family)
  defaults <- list(
    data = standata, chains = chains, parallel_chains = chains,
    iter_warmup = iter_warmup, iter_sampling = iter_sampling,
    seed = seed, refresh = 0, show_messages = FALSE
  )
  if (!is.null(adapt_delta)) defaults$adapt_delta <- adapt_delta
  if (!is.null(max_treedepth)) defaults$max_treedepth <- max_treedepth
  if (length(sample_args)) defaults[names(sample_args)] <- NULL
  do.call(mod$sample, c(defaults, sample_args))
}

#' The precompiled rstan model object for a family
#' @noRd
.cpaic_rstan_model <- function(family) {
  name <- .cpaic_stan_family_file(family)
  mod <- stanmodels[[name]]
  if (is.null(mod)) {
    stop("No compiled Stan model named '", name, "'. The package was ",
         "installed without its Stan models; reinstall cpaic from source.",
         call. = FALSE)
  }
  mod
}

#' Map a cpaic family (plus survival baseline) to its Stan file stem
#' @noRd
.cpaic_stan_family_file <- function(family) paste0("cpaic_", family)

#' Posterior draws as a plain numeric matrix, draws x parameters
#'
#' The column names are the Stan names (`beta[1]`, `gamma[1,2]`), which both
#' engines agree on, so every caller can index by name.
#' @noRd
.cpaic_draws_matrix <- function(fit, variables) {
  if (inherits(fit, "stanfit")) {
    out <- as.matrix(fit, pars = variables)
  } else {
    out <- as.matrix(fit$draws(variables, format = "draws_matrix"))
  }
  res <- matrix(as.numeric(out), nrow = nrow(out), ncol = ncol(out))
  colnames(res) <- colnames(out)
  res
}

#' Posterior draws as an iterations x chains x parameters array
#' @noRd
.cpaic_draws_array <- function(fit, variables) {
  if (inherits(fit, "stanfit")) {
    return(as.array(fit, pars = variables))
  }
  fit$draws(variables, format = "draws_array")
}

#' Convergence summary with identical definitions across backends
#'
#' Both engines are summarized through `posterior`, so `rhat`, `ess_bulk`, and
#' `ess_tail` mean exactly the same thing whichever backend produced the fit.
#' Taking rstan's own `Rhat`/`n_eff` instead would compare different quantities
#' against the same thresholds.
#' @noRd
.cpaic_fit_summary <- function(fit, variables) {
  draws <- tryCatch(
    posterior::as_draws_array(.cpaic_draws_array(fit, variables)),
    error = function(e) NULL)
  if (is.null(draws)) return(NULL)
  as.data.frame(posterior::summarise_draws(
    draws, "rhat", "ess_bulk", "ess_tail"))
}

#' Sampler diagnostics: divergences, tree-depth saturation, and E-BFMI
#' @noRd
.cpaic_sampler_diagnostics <- function(fit) {
  errors <- character()
  capture <- function(name, expr, missing) {
    out <- tryCatch(
      expr,
      error = function(e) {
        errors <<- c(errors, paste0(name, ": ", conditionMessage(e)))
        missing
      }
    )
    if (is.null(out) || !length(out)) {
      errors <<- c(errors, paste0(name, ": returned no values"))
      missing
    } else {
      out
    }
  }
  finish <- function(num_divergent, num_max_treedepth, ebfmi) {
    list(
      num_divergent = num_divergent,
      num_max_treedepth = num_max_treedepth,
      ebfmi = ebfmi,
      unavailable = length(errors) > 0L,
      error = if (length(errors)) paste(errors, collapse = "; ") else NA_character_
    )
  }
  if (inherits(fit, "stanfit")) {
    return(finish(
      capture("get_num_divergent()", rstan::get_num_divergent(fit), NA_integer_),
      capture("get_num_max_treedepth()",
              rstan::get_num_max_treedepth(fit), NA_integer_),
      capture("get_bfmi()", rstan::get_bfmi(fit), NA_real_)
    ))
  }
  d <- capture("diagnostic_summary()", fit$diagnostic_summary(quiet = TRUE), NULL)
  if (is.null(d)) {
    return(finish(NA_integer_, NA_integer_, NA_real_))
  }
  value <- function(name, missing) {
    out <- capture(paste0("diagnostic_summary()$", name), d[[name]], missing)
    if (is.null(out)) {
      errors <<- c(errors, paste0("diagnostic_summary()$", name,
                                  " is unavailable"))
      return(missing)
    }
    out
  }
  finish(
    value("num_divergent", NA_integer_),
    value("num_max_treedepth", NA_integer_),
    value("ebfmi", NA_real_)
  )
}

#' Names of the sampled and generated variables in a fit
#' @noRd
.cpaic_stan_variables <- function(fit) {
  if (inherits(fit, "stanfit")) {
    return(tryCatch(fit@model_pars, error = function(e) character(0)))
  }
  tryCatch(fit$metadata()$stan_variables, error = function(e) character(0))
}

#' Which backend produced this fit?
#' @noRd
.cpaic_fit_backend <- function(fit) {
  if (inherits(fit, "stanfit")) "rstan" else "cmdstanr"
}

#' Every parameter block a cpaic model can sample
#'
#' One list serves both the convergence check and the default of
#' `posterior_summary()`, so the two cannot drift apart and start reporting
#' different parameters. Which of these a given fit actually has depends on the
#' family, the baseline, and `trt_effects`, so callers intersect it with
#' `.cpaic_stan_variables()`.
#' @noRd
.cpaic_param_blocks <- c("mu", "beta", "gamma", "breg", "tau", "sigma",
                         "bsmooth", "bshape_raw", "delta_aux")

#' Posterior summary of a component ML-NMR fit
#'
#' @description
#' Summarizes the posterior draws of a [cmlnmr()] fit: the component effects
#' `beta`, the component by effect-modifier interactions `gamma`, the study
#' baselines `mu`, the heterogeneity `tau`, and any other sampled block, with
#' the usual convergence quantities alongside.
#'
#' Use this rather than reaching into `fit$fit`. That object is whatever the
#' sampler backend returned: an S4 `stanfit` under `backend = "rstan"` and an R6
#' object under `backend = "cmdstanr"`. The two share no accessors, so
#' `fit$fit$summary(...)` works on one and fails on the other. This function
#' works on both and returns the same columns either way.
#'
#' @param x A `cpaic_mlnmr` fit from [cmlnmr()].
#' @param variables Character vector of Stan variable names to summarize, for
#'   example `"tau"` or `c("beta", "gamma")`. Naming a block returns one row per
#'   element of it. The default summarizes every sampled block the fit has.
#' @param ... Further summary functions passed to
#'   `posterior::summarise_draws()`, for example `"quantile2"` or a function.
#'   With none given, the default set is returned.
#'
#' @return A data frame with one row per scalar parameter. With the default
#'   summaries the columns are `variable`, `mean`, `median`, `sd`, `mad`, `q5`,
#'   `q95`, `rhat`, `ess_bulk`, and `ess_tail`. Because both backends are
#'   summarized through the same code, `rhat`, `ess_bulk`, and `ess_tail` are
#'   the same quantities whichever engine produced the fit.
#'
#' @seealso [cmlnmr()] for the fit, [relative_effects()] and
#'   [component_effects()] for effects on the outcome scale rather than the
#'   parameters themselves, and [redact_fit()], which strips the draws.
#'
#' @examples
#' \dontrun{
#' fit <- cmlnmr(ipd, agd, effect_modifiers = "x1", inactive = "Placebo")
#' posterior_summary(fit, "tau")
#' min(posterior_summary(fit)$ess_bulk)
#' }
#' @export
posterior_summary <- function(x, variables = NULL, ...) {
  if (!inherits(x, "cpaic_mlnmr")) {
    stop("`x` must be a cmlnmr() fit.", call. = FALSE)
  }
  fit <- x$fit
  if (is.null(fit)) {
    stop("This fit carries no posterior draws, so it cannot be summarized. ",
         "redact_fit() removes them; refit to get them back.", call. = FALSE)
  }
  have <- .cpaic_stan_variables(fit)
  if (is.null(variables)) {
    variables <- intersect(.cpaic_param_blocks, have)
    if (!length(variables)) variables <- intersect(c("beta", "mu"), have)
  } else {
    variables <- as.character(variables)
    miss <- setdiff(variables, have)
    if (length(miss)) {
      stop("This fit has no parameter named ",
           paste0("`", miss, "`", collapse = ", "),
           ". It has: ", paste(have, collapse = ", "), ".", call. = FALSE)
    }
  }
  draws <- posterior::as_draws_array(.cpaic_draws_array(fit, variables))
  as.data.frame(posterior::summarise_draws(draws, ...))
}
