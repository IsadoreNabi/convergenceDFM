#' Half-life from AR(1) coefficient (safe)
#'
#' Computes the half-life implied by an AR(1) parameter using either Yule-Walker
#' or a manual estimator; guards against degenerate inputs. The half-life is the
#' decay of the deviation envelope: \code{log(0.5) / log(|phi|)} for
#' \code{0 < |phi| < 1}, and \code{Inf} when \code{|phi| >= 1} (no mean reversion).
#'
#' @param x Numeric vector (time series).
#' @return Half-life in periods (numeric scalar, \code{Inf}, or \code{NA}).
#' @keywords internal
#' @noRd

safe_ar1_half_life <- function(x) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]

  if (length(x) < 3) return(NA_real_)

  xm <- x - mean(x, na.rm = TRUE)
  phi <- NA_real_

  fit <- try(stats::ar(xm, aic = FALSE, order.max = 1, method = "yw"), silent = TRUE)

  if (!inherits(fit, "try-error") && length(fit$ar)) {
    phi <- as.numeric(fit$ar[1])
  } else {
    y <- xm[-1]
    z <- xm[-length(xm)]
    phi <- sum(y * z) / sum(z * z)
  }

  if (!is.finite(phi)) return(NA_real_)
  a <- abs(phi)
  if (a >= 1) return(Inf)
  if (a < .Machine$double.eps) return(0)
  log(0.5) / log(a)
}

#' Discrete AR(1) fallback OU estimation (no Stan)
#'
#' Provides a lightweight OU proxy when Stan backends are unavailable, including
#' half-lives for X/Y factors and a simple coupling-strength measure.
#'
#' @param factors_data List with \code{scores_X}, \code{scores_Y}.
#' @param verbose Logical; print progress.
#' @return A list with method tag, half-lives, and coupling strength.
#' @keywords internal
#' @noRd

estimate_factor_OU_fallback <- function(factors_data, verbose = TRUE) {
  if (verbose) {
    message("========================================")
    message("OU FACTOR (DISCRETE FALLBACK without Stan)")
    message("========================================\n")
  }

  FX <- tryCatch(as.matrix(factors_data$scores_X), error = function(e) NULL)
  FY <- tryCatch(as.matrix(factors_data$scores_Y), error = function(e) NULL)

  hlX <- if (!is.null(FX)) apply(FX, 2, safe_ar1_half_life) else numeric(0)
  hlY <- if (!is.null(FY)) apply(FY, 2, safe_ar1_half_life) else numeric(0)

  coupling_strength <- NA_real_

  if (!is.null(FX) && !is.null(FY)) {
    Tn <- min(nrow(FX), nrow(FY))
    if (Tn >= 3) {
      kx <- ncol(FX)
      ky <- ncol(FY)

      X0 <- scale(FX[1:(Tn - 1), , drop = FALSE])
      Y0 <- scale(FY[1:(Tn - 1), , drop = FALSE])
      Y1 <- scale(FY[2:Tn, , drop = FALSE])

      D <- cbind(Y0, X0)

      Bt <- tryCatch({
        solve(crossprod(D), crossprod(D, Y1))
      }, error = function(e) {
        lambda <- 1e-6
        solve(crossprod(D) + diag(lambda, ncol(D)), crossprod(D, Y1))
      })

      Bxy <- Bt[(ncol(Y0) + 1):nrow(Bt), , drop = FALSE]
      coupling_strength <- sqrt(sum(Bxy^2)) / sqrt(kx * ky)
    }
  }

  list(
    method = "discrete_AR1_fallback",
    half_lives_X = hlX,
    half_lives_Y = hlY,
    coupling_strength = coupling_strength
  )
}

#' Half-life from a (possibly non-stationary) AR(1) persistence
#'
#' @param phi_vec Numeric vector of AR(1) persistences.
#' @return Numeric vector of half-lives; \code{Inf} for \code{|phi| >= 1}.
#' @keywords internal
#' @noRd

ou_half_life <- function(phi_vec) {
  phi_vec <- as.numeric(phi_vec)
  vapply(phi_vec, function(ph) {
    if (!is.finite(ph)) return(NA_real_)
    a <- abs(ph)
    if (a >= 1) return(Inf)
    if (a < .Machine$double.eps) return(0)
    log(0.5) / log(a)
  }, numeric(1))
}

#' Estimate Factor Ornstein-Uhlenbeck / AR(1) model (Stan if available)
#'
#' Estimates a Bayesian multivariate, discrete-time mean-reverting model
#' (a first-order vector autoregression with cross-equation coupling) in which
#' \eqn{Y_t} depends on its own lag and on lagged \eqn{X_{t-1}} through a
#' \eqn{\beta} matrix. This is the discrete-time analogue of a coupled
#' Ornstein-Uhlenbeck system; see the vignette section "Methodological notes" for
#' the exact mapping \eqn{\phi = e^{-\kappa \Delta t}}. Uses \pkg{cmdstanr} when
#' available, otherwise \pkg{rstan}, with a discrete AR(1) fallback.
#'
#' @section Convergence is testable, not assumed:
#' The persistence \eqn{\phi} is \emph{not} constrained to \eqn{(0,1)}: it is
#' given a generous support and a weakly-informative, convergence-neutral prior
#' (\code{normal(0.5, 0.5)}), so the posterior can place mass on unit-root or
#' mildly explosive dynamics. Mean reversion (convergence) is therefore a
#' conclusion supported by the data, not an artefact of the parameterization.
#'
#' @param factors_data List with \code{scores_X}, \code{scores_Y}.
#' @param data_prep Optional preprocessed data (reserved).
#' @param chains,iter,seed Stan sampling controls. \code{iter} is the total per
#'   chain; warmup is \code{floor(iter/2)} with a floor that guarantees at least
#'   100 post-warmup draws.
#' @param adapt1,mtd1 \code{adapt_delta} and \code{max_treedepth} for the first
#'   run.
#' @param adapt2,mtd2 \code{adapt_delta} and \code{max_treedepth} used for an
#'   automatic re-run if the first run shows divergent transitions or poor mixing.
#' @param verbose Logical; print progress/details.
#' @return A list with posterior medians (\eqn{\phi}, \eqn{\mu}, \eqn{\beta}),
#'   half-lives, coupling strength, pseudo-R2, MCMC convergence diagnostics
#'   (\code{diagnostics}: max R-hat, min ESS, divergences) and the fitted Stan
#'   object.
#' @examples
#' \donttest{
#'   set.seed(123)
#'   n <- 60
#'   X_scores <- matrix(rnorm(n * 2), n, 2)
#'   Y_scores <- matrix(rnorm(n * 2), n, 2)
#'   factors_data <- list(scores_X = X_scores, scores_Y = Y_scores)
#'   ou_result <- estimate_factor_OU(factors_data, chains = 2, iter = 600,
#'                                   verbose = FALSE)
#'   print(ou_result$half_lives_Y)
#'   print(ou_result$diagnostics)
#' }
#' @export

estimate_factor_OU <- function(factors_data, data_prep = NULL,
                               chains = 4, iter = 2000,
                               seed = 1234,
                               adapt1 = 0.95, adapt2 = 0.999,
                               mtd1 = 12, mtd2 = 15,
                               verbose = TRUE) {

  if (verbose) {
    message("========================================")
    message("FACTOR ORNSTEIN-UHLENBECK / AR(1) MODEL")
    message("========================================\n")
  }

  FX <- tryCatch(as.matrix(factors_data$scores_X), error = function(e) NULL)
  FY <- tryCatch(as.matrix(factors_data$scores_Y), error = function(e) NULL)

  if (is.null(FX) || is.null(FY)) {
    if (verbose) message("WARNING: No X or Y; using fallback.")
    return(estimate_factor_OU_fallback(factors_data, verbose = verbose))
  }

  if (nrow(FX) != nrow(FY)) {
    Tn <- min(nrow(FX), nrow(FY))
    FX <- FX[seq_len(Tn), , drop = FALSE]
    FY <- FY[seq_len(Tn), , drop = FALSE]
  }

  if (nrow(FX) < 5L) {
    if (verbose) message("WARNING: Very short T; using fallback.")
    return(estimate_factor_OU_fallback(factors_data, verbose = verbose))
  }

  nzv <- function(M) {
    v <- apply(M, 2, stats::var, na.rm = TRUE)
    keep <- which(is.finite(v) & v > 1e-10)
    if (!length(keep)) return(matrix(numeric(0), nrow(M), 0))
    M[, keep, drop = FALSE]
  }

  FX <- nzv(FX)
  FY <- nzv(FY)

  if (ncol(FX) == 0L || ncol(FY) == 0L) {
    if (verbose) message("WARNING: Variance approx 0; using fallback.")
    return(estimate_factor_OU_fallback(factors_data, verbose = verbose))
  }

  zscore <- function(M) {
    M <- as.matrix(M)
    storage.mode(M) <- "numeric"
    M[!is.finite(M)] <- NA_real_

    mu <- colMeans(M, na.rm = TRUE)
    sdv <- apply(M, 2, sd, na.rm = TRUE)
    sdv[!is.finite(sdv) | sdv < 1e-8] <- 1

    for (j in seq_len(ncol(M))) {
      v <- M[, j]
      v[is.na(v)] <- mu[j]
      M[, j] <- (v - mu[j]) / sdv[j]
    }

    list(M = M, mu = mu, sd = sdv)
  }

  zsX <- zscore(FX)
  zsY <- zscore(FY)
  X <- zsX$M
  Y <- zsY$M

  data_list <- list(T = nrow(X), Kx = ncol(X), Ky = ncol(Y), X = X, Y = Y)

  # Discrete-time coupled mean-reverting model. phi is NOT confined to (0,1):
  # it can express unit-root / explosive dynamics so that convergence is a
  # testable conclusion. The initial state carries an implicit flat prior
  # (conditional likelihood from t = 2), which keeps the model well defined for
  # |phi| >= 1 (the stationary variance 1/sqrt(1 - phi^2) is not used).
  stan_code <- "
  data {
    int<lower=2> T;
    int<lower=0> Kx;
    int<lower=0> Ky;
    matrix[T, Kx] X;
    matrix[T, Ky] Y;
  }
  parameters {
    vector[Kx] mu_x;
    vector<lower=-1.5, upper=1.5>[Kx] phi_x;
    vector<lower=1e-6>[Kx] sigma_x;

    vector[Ky] mu_y;
    vector<lower=-1.5, upper=1.5>[Ky] phi_y;
    vector<lower=1e-6>[Ky] sigma_y;

    matrix[Ky, Kx] beta;
  }
  model {
    // Weakly-informative, convergence-neutral priors.
    mu_x ~ normal(0, 1);
    mu_y ~ normal(0, 1);
    phi_x ~ normal(0.5, 0.5);   // does not exclude phi >= 1
    phi_y ~ normal(0.5, 0.5);
    sigma_x ~ normal(0, 1);
    sigma_y ~ normal(0, 1);
    to_vector(beta) ~ normal(0, 0.5);

    for (t in 2:T) {
      row_vector[Kx] mean_x = to_row_vector(mu_x) +
                              (X[t-1] - to_row_vector(mu_x)) .* to_row_vector(phi_x);
      X[t] ~ normal(mean_x, to_row_vector(sigma_x));

      row_vector[Ky] mean_y = to_row_vector(mu_y) +
                              (Y[t-1] - to_row_vector(mu_y)) .* to_row_vector(phi_y) +
                              (X[t-1] - to_row_vector(mu_x)) * beta';
      Y[t] ~ normal(mean_y, to_row_vector(sigma_y));
    }
  }"

  have_cmdstan <- requireNamespace("cmdstanr", quietly = TRUE) &&
    isTRUE(tryCatch(!is.null(cmdstanr::cmdstan_version(error_on_NA = FALSE)),
                    error = function(e) FALSE))
  have_rstan <- requireNamespace("rstan", quietly = TRUE)

  if (!have_cmdstan && !have_rstan) {
    if (verbose) message("Stan unavailable. Using fallback.")
    return(estimate_factor_OU_fallback(factors_data, verbose = verbose))
  }

  # Guaranteed-positive warmup/sampling split (fixes the iter <= 1000 bug that
  # produced zero post-warmup draws).
  warmup <- max(100L, floor(iter / 2))
  sampling <- max(100L, iter - warmup)
  n_cores <- max(1L, parallel::detectCores() - 1L)

  ## ---- cmdstanr diagnostics ------------------------------------------------
  diag_cmdstan <- function(fit) {
    out <- list(backend = "cmdstanr", rhat_max = NA_real_,
                ess_bulk_min = NA_real_, ess_tail_min = NA_real_,
                num_divergent = NA_real_, converged = NA)
    if (requireNamespace("posterior", quietly = TRUE)) {
      s <- tryCatch(fit$summary(), error = function(e) NULL)
      if (!is.null(s)) {
        if (!is.null(s$rhat))     out$rhat_max     <- suppressWarnings(max(s$rhat, na.rm = TRUE))
        if (!is.null(s$ess_bulk)) out$ess_bulk_min <- suppressWarnings(min(s$ess_bulk, na.rm = TRUE))
        if (!is.null(s$ess_tail)) out$ess_tail_min <- suppressWarnings(min(s$ess_tail, na.rm = TRUE))
      }
    }
    ds <- tryCatch(fit$diagnostic_summary(quiet = TRUE), error = function(e) NULL)
    if (!is.null(ds) && !is.null(ds$num_divergent)) {
      out$num_divergent <- sum(ds$num_divergent)
    }
    out$converged <- isTRUE(is.finite(out$rhat_max) && out$rhat_max < 1.01 &&
                              (!is.finite(out$num_divergent) || out$num_divergent == 0))
    out
  }

  run_cmdstan <- function(adapt, mtd) {
    stan_path <- cmdstanr::write_stan_file(stan_code)
    mod <- cmdstanr::cmdstan_model(stan_file = stan_path, pedantic = FALSE)
    mod$sample(
      data = data_list,
      chains = chains,
      parallel_chains = min(chains, n_cores),
      iter_warmup = warmup,
      iter_sampling = sampling,
      refresh = 0,
      seed = seed,
      adapt_delta = adapt,
      max_treedepth = mtd,
      show_messages = FALSE
    )
  }

  ## ---- rstan diagnostics ---------------------------------------------------
  diag_rstan <- function(fit) {
    out <- list(backend = "rstan", rhat_max = NA_real_,
                ess_bulk_min = NA_real_, ess_tail_min = NA_real_,
                num_divergent = NA_real_, converged = NA)
    s <- tryCatch(rstan::summary(fit)$summary, error = function(e) NULL)
    if (!is.null(s)) {
      if ("Rhat" %in% colnames(s))  out$rhat_max     <- suppressWarnings(max(s[, "Rhat"], na.rm = TRUE))
      if ("n_eff" %in% colnames(s)) out$ess_bulk_min <- suppressWarnings(min(s[, "n_eff"], na.rm = TRUE))
    }
    sp <- tryCatch(rstan::get_sampler_params(fit, inc_warmup = FALSE),
                   error = function(e) NULL)
    if (!is.null(sp)) {
      out$num_divergent <- sum(vapply(sp, function(x) sum(x[, "divergent__"]), numeric(1)))
    }
    out$converged <- isTRUE(is.finite(out$rhat_max) && out$rhat_max < 1.01 &&
                              (!is.finite(out$num_divergent) || out$num_divergent == 0))
    out
  }

  run_rstan <- function(adapt, mtd, sm) {
    rstan::sampling(
      sm,
      data = data_list,
      chains = chains,
      iter = warmup + sampling,
      warmup = warmup,
      seed = seed,
      control = list(adapt_delta = adapt, max_treedepth = mtd),
      refresh = 0
    )
  }

  res <- NULL
  diagnostics <- NULL

  ## ---- cmdstanr path -------------------------------------------------------
  if (have_cmdstan) {
    tryCatch({
      if (verbose) message("Compiling Stan model (cmdstanr)...")
      if (verbose) message(sprintf("Running MCMC (chains=%d, warmup=%d, sampling=%d)...",
                                   chains, warmup, sampling))
      fit <- run_cmdstan(adapt1, mtd1)
      if (is.null(fit) || !inherits(fit, "CmdStanMCMC")) {
        stop("Stan model did not return a valid fit object")
      }
      d <- diag_cmdstan(fit)
      if (!isTRUE(d$converged)) {
        if (verbose) {
          message(sprintf("  Diagnostics weak (max R-hat=%.3f, divergences=%s); re-running with adapt_delta=%.3f, max_treedepth=%d...",
                          d$rhat_max,
                          ifelse(is.finite(d$num_divergent), d$num_divergent, "NA"),
                          adapt2, mtd2))
        }
        fit2 <- tryCatch(run_cmdstan(adapt2, mtd2), error = function(e) NULL)
        if (!is.null(fit2) && inherits(fit2, "CmdStanMCMC")) {
          d2 <- diag_cmdstan(fit2)
          if (isTRUE(d2$converged) || (is.finite(d2$rhat_max) &&
              (!is.finite(d$rhat_max) || d2$rhat_max <= d$rhat_max))) {
            fit <- fit2; d <- d2
          }
        }
      }
      res <- list(fit = fit, method = "cmdstanr")
      diagnostics <- d
    }, error = function(e) {
      if (verbose) message("WARNING: cmdstanr error: ", e$message)
    })
  }

  ## ---- rstan path ----------------------------------------------------------
  if (is.null(res) && have_rstan) {
    tryCatch({
      if (verbose) message("Trying with rstan...")
      old_rstan_opts <- rstan::rstan_options()
      old_mc_cores <- getOption("mc.cores")
      on.exit({
        rstan::rstan_options(auto_write = old_rstan_opts$auto_write)
        options(mc.cores = old_mc_cores)
      }, add = TRUE)
      rstan::rstan_options(auto_write = TRUE)
      options(mc.cores = n_cores)

      sm <- rstan::stan_model(model_code = stan_code)
      fit <- run_rstan(adapt1, mtd1, sm)
      d <- diag_rstan(fit)
      if (!isTRUE(d$converged)) {
        if (verbose) message("  Diagnostics weak; re-running rstan with higher adapt_delta/treedepth...")
        fit2 <- tryCatch(run_rstan(adapt2, mtd2, sm), error = function(e) NULL)
        if (!is.null(fit2)) {
          d2 <- diag_rstan(fit2)
          if (isTRUE(d2$converged) || (is.finite(d2$rhat_max) &&
              (!is.finite(d$rhat_max) || d2$rhat_max <= d$rhat_max))) {
            fit <- fit2; d <- d2
          }
        }
      }
      res <- list(fit = fit, method = "rstan")
      diagnostics <- d
    }, error = function(e) {
      if (verbose) message("WARNING: rstan error: ", e$message)
    })
  }

  if (is.null(res)) {
    if (verbose) message("Stan could not be executed. Using fallback.")
    return(estimate_factor_OU_fallback(factors_data, verbose = verbose))
  }

  if (!is.null(diagnostics)) {
    if (!isTRUE(diagnostics$converged)) {
      # Always warn: questionable MCMC convergence is critical information.
      warning(sprintf(
        paste0("MCMC convergence is questionable (max R-hat = %.3f, min ESS = %.0f, ",
               "divergences = %s). Interpret the OU estimates with caution; ",
               "increase 'iter' or 'adapt2'."),
        diagnostics$rhat_max, diagnostics$ess_bulk_min,
        ifelse(is.finite(diagnostics$num_divergent), diagnostics$num_divergent, "NA")),
        call. = FALSE)
    } else if (verbose) {
      message(sprintf("MCMC OK: max R-hat = %.3f, min ESS = %.0f, divergences = %s",
                      diagnostics$rhat_max, diagnostics$ess_bulk_min,
                      ifelse(is.finite(diagnostics$num_divergent), diagnostics$num_divergent, "NA")))
    }
  }

  Ky <- ncol(Y)
  Kx <- ncol(X)

  if (res$method == "cmdstanr") {
    if (!requireNamespace("posterior", quietly = TRUE)) {
      stop("Package 'posterior' is required for cmdstanr output. Install it with install.packages('posterior').")
    }

    draws_obj <- tryCatch(res$fit$draws(format = "df"), error = function(e) {
      if (verbose) message("Error extracting draws: ", e$message)
      NULL
    })

    if (is.null(draws_obj)) {
      if (verbose) message("Could not extract draws, using fallback...")
      return(estimate_factor_OU_fallback(factors_data, verbose = verbose))
    }

    df <- as.data.frame(draws_obj)

    extract_vec <- function(prefix, K) {
      cols <- grep(paste0("^", prefix, "\\[\\d+\\]$"), names(df), value = TRUE)
      out <- rep(NA_real_, K)
      if (length(cols)) {
        idx <- as.integer(sub(paste0("^", prefix, "\\[(\\d+)\\]$"), "\\1", cols))
        meds <- vapply(cols, function(nm) stats::median(df[[nm]]), 0.0)
        valid <- which(is.finite(idx) & idx >= 1 & idx <= K)
        if (length(valid)) out[idx[valid]] <- meds[valid]
      }
      out
    }

    mu_x_med <- extract_vec("mu_x", Kx)
    mu_y_med <- extract_vec("mu_y", Ky)
    phi_x_med <- extract_vec("phi_x", Kx)
    phi_y_med <- extract_vec("phi_y", Ky)

    beta_cols <- grep("^beta\\[\\d+,\\d+\\]$", names(df), value = TRUE)
    beta_med <- matrix(0, Ky, Kx)
    if (length(beta_cols)) {
      for (col_name in beta_cols) {
        indices <- gsub("beta\\[|\\]", "", col_name)
        ij <- as.integer(strsplit(indices, ",")[[1]])
        if (length(ij) == 2 && all(is.finite(ij)) &&
            ij[1] >= 1 && ij[1] <= Ky && ij[2] >= 1 && ij[2] <= Kx) {
          beta_med[ij[1], ij[2]] <- median(df[[col_name]])
        }
      }
    }

  } else {
    mu_x_med <- rstan::summary(res$fit, pars = "mu_x")$summary[, "50%"]
    mu_y_med <- rstan::summary(res$fit, pars = "mu_y")$summary[, "50%"]
    phi_x_med <- rstan::summary(res$fit, pars = "phi_x")$summary[, "50%"]
    phi_y_med <- rstan::summary(res$fit, pars = "phi_y")$summary[, "50%"]

    beta_med <- matrix(NA_real_, Ky, Kx)
    for (l in 1:Ky) {
      for (k in 1:Kx) {
        nm <- paste0("beta[", l, ",", k, "]")
        beta_med[l, k] <- rstan::summary(res$fit, pars = nm)$summary[, "50%"]
      }
    }
  }

  hlX <- ou_half_life(phi_x_med)
  hlY <- ou_half_life(phi_y_med)
  coupling_strength <- sqrt(sum(beta_med^2)) / sqrt(length(beta_med))

  Tn <- nrow(X)
  xhat <- yhat <- matrix(NA_real_, Tn, max(Kx, Ky))

  for (t in 2:Tn) {
    xhat[t, 1:Kx] <- mu_x_med + (X[t - 1, ] - mu_x_med) * phi_x_med
    yhat[t, 1:Ky] <- mu_y_med + (Y[t - 1, ] - mu_y_med) * phi_y_med +
      (X[t - 1, ] - mu_x_med) %*% t(beta_med)
  }

  ex <- X[2:Tn, ] - xhat[2:Tn, 1:Kx]
  ey <- Y[2:Tn, ] - yhat[2:Tn, 1:Ky]

  R2_global_X <- 1 - sum(ex^2, na.rm = TRUE) /
    sum((X[2:Tn, ] - matrix(colMeans(X[2:Tn, ]), Tn - 1, Kx, byrow = TRUE))^2,
        na.rm = TRUE)

  R2_global_Y <- 1 - sum(ey^2, na.rm = TRUE) /
    sum((Y[2:Tn, ] - matrix(colMeans(Y[2:Tn, ]), Tn - 1, Ky, byrow = TRUE))^2,
        na.rm = TRUE)

  if (verbose) {
    message("\nPseudo-R2 OU:")
    message("  Global X: ", round(R2_global_X, 3), " | Global Y: ", round(R2_global_Y, 3))
  }

  list(
    method = res$method,
    half_lives_X = hlX,
    half_lives_Y = hlY,
    coupling_strength = coupling_strength,
    phi_x = as.numeric(phi_x_med),
    phi_y = as.numeric(phi_y_med),
    mu_x = as.numeric(mu_x_med),
    mu_y = as.numeric(mu_y_med),
    beta = beta_med,
    stan_fit = res$fit,
    diagnostics = diagnostics,
    r2_ou_global_X = R2_global_X,
    r2_ou_global_Y = R2_global_Y
  )
}
