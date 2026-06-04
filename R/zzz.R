## Package-level state and startup hooks.
##
## `.pkgenv` is a private environment used to hold mutable package options
## (e.g. the logging threshold). It must be a real environment so that
## `assign()`/`get()` work and the binding survives package load. See
## `set_log_level()`/`get_log_level()` in utils.R.

.pkgenv <- new.env(parent = emptyenv())

.onLoad <- function(libname, pkgname) {
  # Initialise mutable package options on load.
  if (is.null(.pkgenv$log_level)) {
    assign("log_level", "INFO", envir = .pkgenv)
  }
}

.onAttach <- function(libname, pkgname) {
  # Single attach hook: announce the package and report Stan availability.
  msg <- paste0(
    "convergenceDFM ", utils::packageVersion("convergenceDFM"),
    " - Dynamic Factor Models for Economic Convergence\n",
    "Type vignette('convergence-analysis') for an introduction"
  )

  stan_status <- if (requireNamespace("cmdstanr", quietly = TRUE)) {
    has_cmdstan <- isTRUE(tryCatch(
      !is.null(cmdstanr::cmdstan_version(error_on_NA = FALSE)),
      error = function(e) FALSE
    ))
    if (has_cmdstan) {
      "Stan backend: CmdStan available for OU estimation."
    } else {
      "Stan backend: cmdstanr installed but CmdStan not configured; using fallback."
    }
  } else if (requireNamespace("rstan", quietly = TRUE)) {
    "Stan backend: RStan available for OU estimation."
  } else {
    "Stan backend: not available; OU estimation uses the discrete fallback."
  }

  packageStartupMessage(msg, "\n", stan_status)
}
