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
#' @noRd
.cpaic_sample <- function(backend, family, standata, chains, iter_warmup,
                          iter_sampling, seed, adapt_delta = NULL,
                          max_treedepth = NULL, sample_args = list()) {
  if (identical(backend, "rstan")) {
    # `...` reaches the sampler unchanged, so an argument that only cmdstanr
    # understands would be handed to rstan::sampling() and kill every chain
    # with "no applicable method for `@`", which says nothing about the cause.
    # Name the problem instead.
    cmdstanr_only <- c("show_exceptions", "show_messages", "parallel_chains",
                       "threads_per_chain", "output_dir", "output_basename",
                       "sig_figs", "iter_warmup", "iter_sampling",
                       "save_warmup", "opencl_ids")
    clash <- intersect(names(sample_args), cmdstanr_only)
    if (length(clash)) {
      stop("`", paste(clash, collapse = "`, `"), "` ",
           if (length(clash) == 1L) "is a cmdstanr argument" else
             "are cmdstanr arguments",
           " and cannot be passed to the rstan backend. Drop ",
           if (length(clash) == 1L) "it" else "them",
           ", or fit with backend = \"cmdstanr\". Note that `adapt_delta` and ",
           "`max_treedepth` are arguments of cmlnmr() itself and work on both.",
           call. = FALSE)
    }
    control <- sample_args$control %||% list()
    if (!is.null(adapt_delta)) control$adapt_delta <- adapt_delta
    if (!is.null(max_treedepth)) control$max_treedepth <- max_treedepth
    sample_args$control <- NULL
    args <- c(
      list(object = .cpaic_rstan_model(family), data = standata,
           chains = chains, warmup = iter_warmup,
           iter = iter_warmup + iter_sampling, seed = seed, refresh = 0),
      if (length(control)) list(control = control),
      sample_args
    )
    return(do.call(rstan::sampling, args))
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
  if (inherits(fit, "stanfit")) {
    return(list(
      num_divergent = tryCatch(rstan::get_num_divergent(fit),
                               error = function(e) 0L),
      num_max_treedepth = tryCatch(rstan::get_num_max_treedepth(fit),
                                   error = function(e) 0L),
      ebfmi = tryCatch(rstan::get_bfmi(fit), error = function(e) NA_real_)))
  }
  d <- tryCatch(fit$diagnostic_summary(quiet = TRUE), error = function(e) NULL)
  if (is.null(d)) return(NULL)
  list(num_divergent = d$num_divergent, num_max_treedepth = d$num_max_treedepth,
       ebfmi = d$ebfmi)
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
