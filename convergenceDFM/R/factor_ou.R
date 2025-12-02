#' Half-life from AR(1) coefficient (safe)
#'
#' Computes the half-life implied by an AR(1) parameter using either Yule-Walker
#' or a manual estimator; guards against degenerate inputs.
#'
#' @param x Numeric vector (time series).
#' @return Half-life in periods (numeric scalar or \code{Inf}).
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
  phi <- max(min(phi, 0.9999), 1e-6)
  
  if (phi <= 0) return(Inf)
  log(0.5) / log(phi)
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
    message("OU FACTOR (DISCREET FALLBACK without Stan)")
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
      
      X0 <- scale(FX[1:(Tn-1), , drop = FALSE])
      Y0 <- scale(FY[1:(Tn-1), , drop = FALSE])
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

#' Estimate Factor Ornstein-Uhlenbeck model (Stan if available)
#'
#' Estimates a multivariate OU with cross-equation coupling \eqn{Y_{t}} depending
#' on lagged \eqn{X_{t-1}} via a \eqn{\beta} matrix. Uses \pkg{cmdstanr} when
#' available, otherwise \pkg{rstan}, with a discrete AR(1) fallback.
#'
#' @param factors_data List with \code{scores_X}, \code{scores_Y}.
#' @param data_prep Optional preprocessed data (reserved).
#' @param chains,iter,seed Stan sampling controls.
#' @param adapt1,adapt2,mtd1,mtd2 Advanced Stan controls.
#' @param verbose Logical; print progress/details.
#' @return A list with posterior medians (\eqn{\phi}, \eqn{\mu}, \eqn{\beta}),
#'   half-lives, coupling strength, pseudo-R2, and the fitted Stan object.
#' @examples
#' \donttest{
#'   # Create toy factor data
#'   set.seed(123)
#'   n <- 50
#'   X_scores <- matrix(rnorm(n * 2), n, 2)
#'   Y_scores <- matrix(rnorm(n * 2), n, 2)
#'   factors_data <- list(scores_X = X_scores, scores_Y = Y_scores)
#'   
#'   # Estimate OU model (reduce iterations for speed)
#'   ou_result <- estimate_factor_OU(factors_data, chains = 2, iter = 500, 
#'                                    verbose = FALSE)
#'   
#'   # Check half-lives
#'   print(ou_result$half_lives_Y)
#' }
#' @export

estimate_factor_OU <- function(factors_data, data_prep = NULL,
                               chains = 4, iter = 2000,
                               seed = 1234,
                               adapt1 = 0.98, adapt2 = 0.999,
                               mtd1 = 12, mtd2 = 15,
                               verbose = TRUE) {
  
  if (verbose) {
    message("========================================")
    message("FACTOR ORNSTEIN-UHLENBECK MODEL")
    message("========================================\n")
  }
  
  FX <- tryCatch(as.matrix(factors_data$scores_X), error = function(e) NULL)
  FY <- tryCatch(as.matrix(factors_data$scores_Y), error = function(e) NULL)
  
  if (is.null(FX) || is.null(FY)) {
    if (verbose) message("WARNING: No X or Y; use fallback.")
    return(estimate_factor_OU_fallback(factors_data, verbose = verbose))
  }
  
  if (nrow(FX) != nrow(FY)) {
    Tn <- min(nrow(FX), nrow(FY))
    FX <- FX[seq_len(Tn), , drop = FALSE]
    FY <- FY[seq_len(Tn), , drop = FALSE]
  }
  
  if (nrow(FX) < 5L) {
    if (verbose) message("WARNING: Very short T; use fallback.")
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
    if (verbose) message("WARNING: Var approx 0; use fallback.")
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
    vector[Kx] phi_x_raw;
    vector<lower=1e-6>[Kx] sigma_x;
    
    vector[Ky] mu_y;
    vector[Ky] phi_y_raw;
    vector<lower=1e-6>[Ky] sigma_y;
    
    matrix[Ky, Kx] beta;
  }
  transformed parameters {
    real eps = 1e-6;
    vector[Kx] phi_x = inv_logit(phi_x_raw) * (1 - 2*eps) + eps;
    vector[Ky] phi_y = inv_logit(phi_y_raw) * (1 - 2*eps) + eps;
  }
  model {
    mu_x ~ normal(0, 1);
    mu_y ~ normal(0, 1);
    phi_x_raw ~ normal(2, 1);
    phi_y_raw ~ normal(2, 1);
    sigma_x ~ normal(0, 1);
    sigma_y ~ normal(0, 1);
    to_vector(beta) ~ normal(0, 0.5);
    
    for (k in 1:Kx)
      X[1, k] ~ normal(mu_x[k], sigma_x[k] / sqrt(1 - square(phi_x[k])));
    for (l in 1:Ky)
      Y[1, l] ~ normal(mu_y[l], sigma_y[l] / sqrt(1 - square(phi_y[l])));
    
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
    if (verbose) message("ERROR: Stan unavailable. Using fallback.")
    return(estimate_factor_OU_fallback(factors_data, verbose = verbose))
  }
  
  compute_hl <- function(phi_vec) {
    phi_vec <- as.numeric(phi_vec)
    phi_vec[!is.finite(phi_vec)] <- NA_real_
    ifelse(!is.finite(phi_vec), NA_real_,
           ifelse(phi_vec <= 0, Inf, log(0.5) / log(phi_vec)))
  }
  
  res <- NULL
  
  if (have_cmdstan) {
    tryCatch({
      if (verbose) message("Compiling Stan model (cmdstanr)...")
      stan_path <- cmdstanr::write_stan_file(stan_code)
      mod <- cmdstanr::cmdstan_model(stan_file = stan_path, pedantic = FALSE)
      
      if (verbose) {
        message(sprintf("Running MCMC (chains=%d, iter=%d)...", chains, iter))
      }
      
      warmup <- max(500L, floor(iter / 2))
      sampling <- iter - warmup
      
      fit <- mod$sample(
        data = data_list,
        chains = chains,
        parallel_chains = min(chains, max(1L, parallel::detectCores() - 1L)),
        iter_warmup = warmup,
        iter_sampling = sampling,
        refresh = 0,
        seed = seed,
        adapt_delta = adapt1,
        max_treedepth = mtd1,
        show_messages = FALSE
      )
      
      # CRITICAL FIX: Verify fit is valid before continuing
      if (is.null(fit) || !inherits(fit, "CmdStanMCMC")) {
        stop("Stan model did not return a valid fit object")
      }
      
      # Verify draws() method exists
      if (!is.function(fit$draws)) {
        stop("Fit object does not have draws() method")
      }
      
      res <- list(fit = fit, method = "cmdstanr")
      
    }, error = function(e) {
      if (verbose) message("WARNING: cmdstanr error: ", e$message)
      res <- NULL
    })
  }
  
  if (is.null(res) && have_rstan) {
    tryCatch({
      if (verbose) message("Trying with rstan...")
      
      # Save original options
      old_rstan_opts <- rstan::rstan_options()
      old_mc_cores <- getOption("mc.cores")
      on.exit({
        rstan::rstan_options(auto_write = old_rstan_opts$auto_write)
        options(mc.cores = old_mc_cores)
      }, add = TRUE)
      
      # Set options
      rstan::rstan_options(auto_write = TRUE)
      options(mc.cores = max(1L, parallel::detectCores() - 1L))
      
      sm <- rstan::stan_model(model_code = stan_code)
      
      fit <- rstan::sampling(
        sm,
        data = data_list,
        chains = chains,
        iter = iter,
        warmup = max(500L, floor(iter / 2)),
        seed = seed,
        control = list(adapt_delta = adapt1, max_treedepth = mtd1),
        refresh = 0
      )
      
      res <- list(fit = fit, method = "rstan")
      
    }, error = function(e) {
      if (verbose) message("WARNING: rstan error: ", e$message)
    })
  }
  
  if (is.null(res)) {
    if (verbose) message("ERROR: Stan could not be executed. Using fallback.")
    return(estimate_factor_OU_fallback(factors_data, verbose = verbose))
  }
  
  Ky <- ncol(Y)
  Kx <- ncol(X)
  
  if (res$method == "cmdstanr") {
    # CRITICAL FIX: Verify posterior package is available
    if (!requireNamespace("posterior", quietly = TRUE)) {
      stop("Package 'posterior' is required. Install it with: install.packages('posterior')")
    }
    
    # Extract draws safely
    draws_obj <- tryCatch({
      res$fit$draws(format = "df")
    }, error = function(e) {
      if (verbose) message("Error extracting draws: ", e$message)
      NULL
    })
    
    if (is.null(draws_obj)) {
      if (verbose) message("Could not extract draws, trying fallback...")
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
        if (length(ij) == 2 && all(is.finite(ij))) {
          if (ij[1] >= 1 && ij[1] <= Ky && ij[2] >= 1 && ij[2] <= Kx) {
            beta_med[ij[1], ij[2]] <- median(df[[col_name]])
          }
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
  
  hlX <- compute_hl(phi_x_med)
  hlY <- compute_hl(phi_y_med)
  coupling_strength <- sqrt(sum(beta_med^2)) / sqrt(length(beta_med))
  
  Tn <- nrow(X)
  xhat <- yhat <- matrix(NA_real_, Tn, max(Kx, Ky))
  
  for (t in 2:Tn) {
    xhat[t, 1:Kx] <- mu_x_med + (X[t-1, ] - mu_x_med) * phi_x_med
    yhat[t, 1:Ky] <- mu_y_med + (Y[t-1, ] - mu_y_med) * phi_y_med + 
      (X[t-1, ] - mu_x_med) %*% t(beta_med)
  }
  
  ex <- X[2:Tn, ] - xhat[2:Tn, 1:Kx]
  ey <- Y[2:Tn, ] - yhat[2:Tn, 1:Ky]
  
  R2_global_X <- 1 - sum(ex^2, na.rm = TRUE) / 
    sum((X[2:Tn, ] - matrix(colMeans(X[2:Tn, ]), 
                            Tn - 1, Kx, byrow = TRUE))^2, 
        na.rm = TRUE)
  
  R2_global_Y <- 1 - sum(ey^2, na.rm = TRUE) / 
    sum((Y[2:Tn, ] - matrix(colMeans(Y[2:Tn, ]), 
                            Tn - 1, Ky, byrow = TRUE))^2, 
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
    r2_ou_global_X = R2_global_X,
    r2_ou_global_Y = R2_global_Y
  )
}