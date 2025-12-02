.onAttach <- function(libname, pkgname) {
  # Check for optional dependencies
  if (requireNamespace("cmdstanr", quietly = TRUE)) {
    tryCatch({
      cmdstanr::cmdstan_version(error_on_NA = FALSE)
      packageStartupMessage("CmdStan is available for advanced OU estimation")
    }, error = function(e) {
      packageStartupMessage("CmdStan not configured, using fallback methods")
    })
  } else if (requireNamespace("rstan", quietly = TRUE)) {
    packageStartupMessage("RStan is available for OU estimation")
  } else {
    packageStartupMessage("Stan not available, using fallback OU methods")
  }
}

.onAttach <- function(libname, pkgname) {
  packageStartupMessage(
    "convergenceDFM ", utils::packageVersion("convergenceDFM"),
    " - Dynamic Factor Models for Economic Convergence\n",
    "Type vignette('convergence-analysis') for an introduction"
  )
}
