#' Project vector to probability simplex
#'
#' Projects a numeric vector onto the probability simplex (non-negative values
#' summing to 1) using Euclidean projection.
#'
#' @param v Numeric vector to be projected.
#'
#' @return Numeric vector on the simplex.
#'
#' @keywords internal
#' @noRd

project_to_simplex <- function(v) {
  n <- length(v)
  mu <- sort(v, decreasing = TRUE)
  cs <- cumsum(mu)
  rho <- max(which((mu - (cs - 1) / (1:n)) > 0))
  theta <- (cs[rho] - 1) / rho
  w <- pmax(v - theta, 0)
  w / sum(w)
}

#' Generate alternative prior schemes
#'
#' Creates multiple alternative prior probability schemes based on different
#' transformations of the base weight matrix for sensitivity analysis.
#'
#' @param mip_weights Numeric matrix of baseline sectoral weights.
#' @param n_sectors Integer. Number of sectors.
#' @param seed Optional integer; seeds the random prior for reproducibility. The
#'   RNG state is restored on exit. Default \code{NULL}.
#' @param verbose Logical; print progress and diagnostic information. Default \code{TRUE}.
#'
#' @return List of alternative prior matrices.
#'
#' @keywords internal
#' @noRd

generate_alternative_priors <- function(mip_weights, n_sectors, seed = NULL,
                                        verbose = TRUE) {
  if (!is.null(seed)) {
    if (exists(".Random.seed", envir = .GlobalEnv)) {
      old_seed <- get(".Random.seed", envir = .GlobalEnv)
      on.exit(assign(".Random.seed", old_seed, envir = .GlobalEnv), add = TRUE)
    }
    set.seed(seed)
  }

  if (verbose) {
    message("\n========================================")
    message("GENERATION OF ALTERNATIVE PRIORS")
    message("========================================\n")
  }

  priors <- list()

  priors$equal <- rep(1 / n_sectors, n_sectors)
  if (verbose) {
    message("OK: Generated equiprobable prior")
  }

  random_vec <- runif(n_sectors)
  random_vec <- random_vec / sum(random_vec)

  # Correct Euclidean projection onto mip_weights: (a.b / ||b||^2) b.
  denom <- sum(mip_weights^2)
  proj_mip <- if (denom > 0) (sum(random_vec * mip_weights) / denom) * mip_weights
              else rep(0, length(mip_weights))
  orthogonal <- random_vec - proj_mip

  priors$orthogonal <- project_to_simplex(orthogonal + 0.1)
  
  cor_check <- cor(priors$orthogonal, mip_weights)
  if (verbose) {
    message(sprintf("OK: Misaligned Prior (cor = %.3f)", cor_check))
  }
  
  for (lambda in c(0.25, 0.5, 0.75)) {
    priors[[paste0("shrink_", lambda)]] <- 
      lambda * mip_weights + (1 - lambda) * priors$equal
    priors[[paste0("shrink_", lambda)]] <- 
      priors[[paste0("shrink_", lambda)]] / sum(priors[[paste0("shrink_", lambda)]])
  }
  if (verbose) {
    message("OK: Shrinkage priors generated")
  }
  
  priors
}

#' Permutation-based robustness test for X->Y coupling
#'
#' Tests whether the estimated X->Y coupling strength is stronger than expected
#' when the temporal alignment between X and Y is broken. The null is generated
#' by circularly time-shifting the Y factors, which preserves Y's own
#' autocorrelation while destroying the cross-series relationship.
#'
#' @section Why not a column permutation:
#' The previous version permuted the \emph{columns} of the Y factors. The coupling
#' strength is a Frobenius norm of the coupling matrix and is invariant to such a
#' relabelling, so the old null was degenerate (every replicate reproduced the
#' baseline). Shuffling time, not factor labels, is what tests coupling.
#'
#' @return List with components:
#'   \describe{
#'     \item{\code{baseline}}{Observed coupling strength.}
#'     \item{\code{permuted_mean}, \code{permuted_median}, \code{permuted_sd}}{
#'       Summary of the null coupling distribution.}
#'     \item{\code{p_value}}{Monte Carlo one-sided p-value
#'       \code{(1 + #{null >= baseline}) / (n + 1)}.}
#'     \item{\code{z_score}}{Standardized effect size.}
#'     \item{\code{robust}}{\code{TRUE} if \code{p_value < 0.05}.}
#'     \item{\code{n_successful}}{Number of usable replicates.}
#'   }
#'
#' @param factors_data List with \code{scores_X}, \code{scores_Y}.
#' @param data_prep Prepared data object (passed to the OU estimator).
#' @param n_perms Number of permutations (default: 100).
#' @param seed Random seed for reproducibility (default: 123). Now honoured.
#' @param use_stan Logical, use Stan for estimation (default: TRUE).
#' @param chains Number of MCMC chains (default: 4).
#' @param iter Number of MCMC iterations (default: 2000).
#' @param verbose Logical; print progress and diagnostic information. Default \code{TRUE}.
#' @export

test_permutation_robustness <- function(factors_data, data_prep,
                                        n_perms = 100, seed = 123,
                                        use_stan = TRUE,
                                        chains = 4, iter = 2000,
                                        verbose = TRUE) {
  if (!is.null(seed)) {
    if (exists(".Random.seed", envir = .GlobalEnv)) {
      old_seed <- get(".Random.seed", envir = .GlobalEnv)
      on.exit(assign(".Random.seed", old_seed, envir = .GlobalEnv), add = TRUE)
    }
    set.seed(seed)
  }

  if (verbose) {
    message("\n========================================")
    message("TIME-SHIFT COUPLING PERMUTATION TEST")
    message("========================================\n")
    message("Computing baseline on the original data...")
  }
  
  if (use_stan) {
    ou_baseline <- tryCatch(
      estimate_factor_OU(factors_data, data_prep,
                         chains = chains, iter = iter,
                         verbose = FALSE),
      error = function(e) {
        if (verbose) message("WARNING: Error in Stan, using fallback")
        estimate_factor_OU_fallback(factors_data, verbose = verbose)
      }
    )
  } else {
    ou_baseline <- estimate_factor_OU_fallback(factors_data, verbose = verbose)
  }
  
  if (is.null(ou_baseline)) {
    if (verbose) message("FAIL: Baseline could not be calculated")
    return(NULL)
  }
  
  baseline_coupling <- ou_baseline$coupling_strength
  if (verbose) {
    message(sprintf("Baseline coupling: %.4f", baseline_coupling))
  }
  
  if (verbose) {
    message(sprintf("\nRunning %d permutations", n_perms))
    if (use_stan) {
      message(" (It may take several minutes)...")
    } else {
      message("...")
    }
  }
  
  perm_couplings <- numeric(n_perms)
  successful_perms <- 0
  
  if (verbose) {
    pb <- txtProgressBar(min = 0, max = n_perms, style = 3)
  }
  
  Ty <- nrow(as.matrix(factors_data$scores_Y))
  for (i in 1:n_perms) {
    # Break the X->Y temporal alignment by circularly shifting Y in time.
    tau <- sample.int(max(2L, Ty) - 1L, 1L)
    factors_perm <- factors_data
    factors_perm$scores_Y <- circ_shift_rows(as.matrix(factors_data$scores_Y), tau)

    if (use_stan) {
      ou_perm <- tryCatch(
        suppressMessages(
          estimate_factor_OU(factors_perm, data_prep,
                             chains = chains, iter = iter,
                             verbose = FALSE)
        ),
        error = function(e) {
          if (verbose) message("WARNING: Stan failed in permutation", i, ", using fallback")
          estimate_factor_OU_fallback(factors_perm, verbose = FALSE)
        }
      )
    } else {
      ou_perm <- tryCatch(
        suppressMessages(estimate_factor_OU_fallback(factors_perm, verbose = FALSE)),
        error = function(e) NULL
      )
    }
    
    if (!is.null(ou_perm) && !is.null(ou_perm$coupling_strength)) {
      perm_couplings[i] <- ou_perm$coupling_strength
      successful_perms <- successful_perms + 1
    } else {
      perm_couplings[i] <- NA
    }
    
    if (verbose) {
      setTxtProgressBar(pb, i)
    }
  }
  
  if (verbose) {
    close(pb)
  }
  
  valid_perms <- perm_couplings[!is.na(perm_couplings)]
  if (verbose) {
    message(sprintf("\n\nSuccessful permutations: %d/%d", successful_perms, n_perms))
  }

  if (length(valid_perms) < 10) {
    if (verbose) message("FAIL: Insufficient successful permutations")
    return(NULL)
  }

  mean_perm <- mean(valid_perms)
  sd_perm <- sd(valid_perms)
  median_perm <- median(valid_perms)

  # Monte Carlo one-sided p-value with the (1 + .)/(n + 1) correction.
  p_value <- (1 + sum(valid_perms >= baseline_coupling)) / (length(valid_perms) + 1)
  z_score <- (baseline_coupling - mean_perm) / (sd_perm + 1e-10)
  robust <- isTRUE(p_value < 0.05)

  if (verbose) {
    message("\n=====================================")
    message("RESULTS:")
    message(sprintf("  Original coupling: %.4f", baseline_coupling))
    message(sprintf("  Null coupling: mean %.4f (SD %.4f), median %.4f",
                    mean_perm, sd_perm, median_perm))
    message(sprintf("  p-value: %.4f | z-score: %.2f", p_value, z_score))
    message(if (robust) "  -> Coupling exceeds the time-shift null (p < 0.05)"
            else        "  -> Coupling NOT distinguishable from the null (p >= 0.05)")
  }
  
  list(
    baseline = baseline_coupling,
    permuted_mean = mean_perm,
    permuted_median = median_perm,
    permuted_sd = sd_perm,
    p_value = p_value,
    z_score = z_score,
    robust = robust,
    n_successful = successful_perms
  )
}

#' Reweighting-based robustness test
#'
#' Tests sensitivity of the coupling result to alternative sectoral weighting
#' schemes by re-running the weight-blend disaggregation and the (fallback)
#' factor-OU pipeline with perturbed weights. The coefficient of variation of the
#' coupling across schemes summarizes the sensitivity.
#'
#' @return List with components:
#'   \describe{
#'     \item{\code{results}}{Per-scheme list with \code{coupling} and
#'       \code{half_life_mean}.}
#'     \item{\code{cv_coupling}}{Coefficient of variation of coupling across
#'       schemes.}
#'     \item{\code{robust}}{\code{TRUE} if \code{cv_coupling < 0.3}.}
#'   }
#'
#' @param path_cpi Path to CPI data file
#' @param path_weights Path to weights data file
#' @param X_matrix Matrix of variables
#' @param max_comp Maximum number of components (default: 3)
#' @param seed Optional integer for reproducible priors (default: 123).
#' @param verbose Logical; print progress and diagnostic information. Default \code{TRUE}.
#' @export

test_reweighting_robustness <- function(path_cpi, path_weights, X_matrix,
                                        max_comp = 3, seed = 123, verbose = TRUE) {
  if (verbose) {
    message("\n========================================")
    message("ROBUSTNESS TEST BY REWEIGHTING")
    message("========================================\n")
  }

  if (!file.exists(path_cpi) || !file.exists(path_weights)) {
    if (verbose) message("FAIL: Files not found")
    return(NULL)
  }

  weights_data <- tryCatch(
    read_weights_matrix(path_weights),
    error = function(e) {
      if (verbose) message("FAIL: Error reading weights: ", e$message)
      NULL
    }
  )

  if (is.null(weights_data)) {
    if (verbose) message("FAIL: The weights could not be read")
    return(NULL)
  }
  
  mip_weights <- colMeans(weights_data$P)
  n_sectors <- length(mip_weights)
  
  if (verbose) {
    message("Sectors:", n_sectors)
    message("Years:", length(weights_data$years), "\n")
  }
  
  alt_priors <- generate_alternative_priors(mip_weights, n_sectors,
                                            seed = seed, verbose = verbose)

  results_reweight <- list()
  
  for (prior_name in names(alt_priors)) {
    if (verbose) {
      message("\n----------------------------------------")
      message("Processing prior:", prior_name)
    }
    
    prior_weights <- alt_priors[[prior_name]]
    
    tryCatch({
      disagg_result <- run_disaggregation_custom_prior(
        path_cpi, path_weights,
        custom_prior = prior_weights,
        method = "weighted",
        lambda = 0.7,
        verbose = verbose
      )
      
      if (is.null(disagg_result$Y_matrix)) {
        if (verbose) message("  FAIL: Disaggregation failed")
        results_reweight[[prior_name]] <- list(coupling = NA)
        next
      }
      
      Y_reweighted <- disagg_result$Y_matrix
      
      min_rows <- min(nrow(Y_reweighted), nrow(X_matrix))
      min_cols <- min(ncol(Y_reweighted), ncol(X_matrix))
      
      Y_reweighted <- Y_reweighted[1:min_rows, 1:min_cols]
      X_matrix_adj <- X_matrix[1:min_rows, 1:min_cols]
      
      data_clean <- diagnose_data(X_matrix_adj, Y_reweighted, verbose = verbose)
      
      data_prep <- list(
        X_scaled = scale(data_clean$X_matrix),
        Y_scaled = scale(data_clean$Y_matrix)
      )
      
      selection <- suppressWarnings(
        select_optimal_components_safe(
          data_prep$X_scaled,
          data_prep$Y_scaled,
          max_comp = max_comp,
          verbose = FALSE
        )
      )
      
      factors_data <- list(
        scores_X = selection$pls_model$scores[, 1:selection$optimal_ncomp, drop = FALSE],
        scores_Y = selection$pls_model$Yscores[, 1:selection$optimal_ncomp, drop = FALSE]
      )
      
      ou_result <- suppressMessages(
        estimate_factor_OU_fallback(factors_data, verbose = FALSE)
      )
      
      results_reweight[[prior_name]] <- list(
        coupling = ou_result$coupling_strength,
        half_life_mean = mean(c(ou_result$half_lives_X, 
                                ou_result$half_lives_Y), na.rm = TRUE)
      )
      
      if (verbose) {
        message(sprintf("  OK: Coupling: %.4f", ou_result$coupling_strength))
      }
      
    }, error = function(e) {
      if (verbose) message("  FAIL: Error:", e$message)
      results_reweight[[prior_name]] <- list(coupling = NA)
    })
  }
  
  if (verbose) {
    message("\n========================================")
    message("STABILITY SUMMARY")
    message("========================================\n")
  }
  
  couplings <- sapply(results_reweight, function(x) x$coupling)
  valid_couplings <- couplings[!is.na(couplings)]
  
  if (length(valid_couplings) < 2) {
    if (verbose) message("FAIL: Insufficient valid results")
    return(list(results = results_reweight, robust = FALSE))
  }
  
  if (verbose) {
    message("X->Y Couplings:")
    message(paste(names(couplings), round(couplings, 4), sep = ": ", collapse = "\n"))
  }
  
  cv_coupling <- sd(valid_couplings) / mean(valid_couplings)
  if (verbose) {
    message(sprintf("\nCV (acoplamiento): %.2f%%", cv_coupling * 100))
  }
  
  if (cv_coupling < 0.3) {
    if (verbose) message("OK: ROBUST to reweighting (CV < 30%)")
  } else {
    if (verbose) message("WARNING: Sensitive to reweighting (CV > 30%)")
  }
  
  list(
    results = results_reweight,
    cv_coupling = cv_coupling,
    robust = cv_coupling < 0.3
  )
}

#' Delete-one-sector jackknife of the X->Y coupling
#'
#' Genuine leave-one-sector-out jackknife: each sector (column) is dropped in
#' turn from both X and Y, the coupling is re-estimated, and the delete-one
#' replicates are used to compute the jackknife bias and standard error and to
#' rank sectors by influence. (The previous version dropped the top-\code{k}
#' highest-variance sectors a single time and was not a jackknife.)
#'
#' @return List with components:
#'   \describe{
#'     \item{\code{baseline}}{Full-sample coupling.}
#'     \item{\code{jackknife_estimates}}{Named vector of delete-one coupling
#'       estimates.}
#'     \item{\code{bias}}{Jackknife bias estimate \code{(n-1)(mean(jk) - full)}.}
#'     \item{\code{se}}{Jackknife standard error.}
#'     \item{\code{influence}}{Per-sector change \code{full - jk_i} (signed).}
#'     \item{\code{retention}}{Per-sector \code{jk_i / full}.}
#'     \item{\code{influential_sectors}}{The \code{k_exclude} most influential
#'       sectors.}
#'     \item{\code{robust}}{\code{TRUE} if no single sector changes the coupling
#'       by more than 50 percent.}
#'   }
#'
#' @param X_matrix Matrix of first set of variables (sectors in columns)
#' @param Y_matrix Matrix of second set of variables (sectors in columns)
#' @param sector_names Vector of sector names (default: column names of Y)
#' @param k_exclude Number of most-influential sectors to report (default: 3)
#' @param max_comp Maximum number of components (default: 3)
#' @param seed Optional integer for reproducible component selection/diagnostics.
#' @param verbose Logical; print progress and diagnostic information. Default \code{TRUE}.
#' @export

test_jackknife_sectors <- function(X_matrix, Y_matrix, sector_names = NULL,
                                   k_exclude = 3, max_comp = 3, seed = 123,
                                   verbose = TRUE) {
  if (verbose) {
    message("\n========================================")
    message("DELETE-ONE-SECTOR JACKKNIFE")
    message("========================================\n")
  }

  X_matrix <- as.matrix(X_matrix)
  Y_matrix <- as.matrix(Y_matrix)
  n_sec <- ncol(Y_matrix)
  if (is.null(sector_names)) {
    sector_names <- colnames(Y_matrix)
    if (is.null(sector_names)) sector_names <- paste0("S", seq_len(n_sec))
  }

  coupling_of <- function(Xm, Ym) {
    if (ncol(Xm) < 1 || ncol(Ym) < 1 || nrow(Xm) < 5) return(NA_real_)
    dc <- tryCatch(diagnose_data(Xm, Ym, verbose = FALSE, seed = seed),
                   error = function(e) NULL)
    if (is.null(dc)) return(NA_real_)
    sel <- tryCatch(suppressMessages(
      select_optimal_components_safe(scale(dc$X_matrix), scale(dc$Y_matrix),
                                     max_comp = max_comp, seed = seed, verbose = FALSE)),
      error = function(e) NULL)
    if (is.null(sel)) return(NA_real_)
    k <- sel$optimal_ncomp
    fac <- list(
      scores_X = sel$pls_model$scores[, 1:k, drop = FALSE],
      scores_Y = sel$pls_model$Yscores[, 1:k, drop = FALSE]
    )
    ou <- tryCatch(suppressMessages(estimate_factor_OU_fallback(fac, verbose = FALSE)),
                   error = function(e) NULL)
    if (is.null(ou)) NA_real_ else ou$coupling_strength
  }

  baseline <- coupling_of(X_matrix, Y_matrix)
  if (!is.finite(baseline)) {
    if (verbose) message("FAIL: baseline coupling could not be computed")
    return(NULL)
  }

  jk <- rep(NA_real_, n_sec); names(jk) <- sector_names
  if (verbose) pb <- txtProgressBar(min = 0, max = n_sec, style = 3)
  for (j in seq_len(n_sec)) {
    Xj <- if (j <= ncol(X_matrix)) X_matrix[, -j, drop = FALSE] else X_matrix
    Yj <- Y_matrix[, -j, drop = FALSE]
    jk[j] <- coupling_of(Xj, Yj)
    if (verbose) setTxtProgressBar(pb, j)
  }
  if (verbose) close(pb)

  valid <- is.finite(jk)
  jkv <- jk[valid]
  n <- length(jkv)
  if (n < 2) {
    if (verbose) message("FAIL: too few valid jackknife replicates")
    return(list(baseline = baseline, jackknife_estimates = jk, robust = NA))
  }

  jack_mean <- mean(jkv)
  bias <- (n - 1) * (jack_mean - baseline)
  se   <- sqrt((n - 1) / n * sum((jkv - jack_mean)^2))
  influence <- baseline - jk                      # signed
  retention <- jk / baseline
  rel_change <- abs(influence) / (abs(baseline) + 1e-12)

  ord <- order(abs(influence), decreasing = TRUE, na.last = NA)
  influential_sectors <- sector_names[ord][seq_len(min(k_exclude, length(ord)))]
  robust <- isTRUE(max(rel_change, na.rm = TRUE) < 0.5)

  if (verbose) {
    message(sprintf("\nFull-sample coupling: %.4f", baseline))
    message(sprintf("Jackknife bias: %.4f | SE: %.4f", bias, se))
    message("Most influential sectors: ", paste(influential_sectors, collapse = ", "))
    message(if (robust) "OK: no single sector dominates (<50% change)"
            else        "WARNING: at least one sector changes coupling by >50%")
  }

  list(
    baseline = baseline,
    jackknife_estimates = jk,
    bias = bias,
    se = se,
    influence = influence,
    retention = retention,
    influential_sectors = influential_sectors,
    robust = robust
  )
}