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
#' @param verbose Logical; print progress and diagnostic information. Default \code{TRUE}.
#'
#' @return List of alternative prior matrices.
#'
#' @keywords internal
#' @noRd

generate_alternative_priors <- function(mip_weights, n_sectors, verbose = TRUE) {
  if (verbose) {
    message("\n========================================")
    message("GENERATION OF ALTERNATIVE PRIORS")
    message("========================================\n")
  }
  
  priors <- list()
  
  priors$equal <- rep(1/n_sectors, n_sectors)
  if (verbose) {
    message("OK: Generated equiprobable prior")
  }
  
  random_vec <- runif(n_sectors)
  random_vec <- random_vec / sum(random_vec)
  
  proj_mip <- sum(random_vec * mip_weights) * mip_weights
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

#' Permutation-based robustness test
#'
#' Tests the robustness of factor-OU convergence findings by randomly permuting
#' the Y factor space and re-estimating the model. Generates empirical null
#' distribution for convergence statistics.
#'
#'   (too slow for many iterations).
#'
#' @return List with components:
#'   \describe{
#'     \item{\code{observed_lambda}}{Original mean-reversion speeds.}
#'     \item{\code{null_distribution}}{Matrix of permutation-based lambda values.}
#'     \item{\code{p_values}}{One-sided p-values for each factor.}
#'     \item{\code{significant}}{Logical vector indicating significance at alpha = 0.05.}
#'     \item{\code{effect_size}}{Standardized effect sizes (z-scores).}
#'   }
#'
#'
#' @param factors_data Data frame with factor information
#' @param data_prep Prepared data object
#' @param n_perms Number of permutations (default: 100)
#' @param seed Random seed for reproducibility (default: 123)
#' @param use_stan Logical, use Stan for estimation (default: TRUE)
#' @param chains Number of MCMC chains (default: 4)
#' @param iter Number of MCMC iterations (default: 2000)
#' @param verbose Logical; print progress and diagnostic information. Default \code{TRUE}.
#' @export

test_permutation_robustness <- function(factors_data, data_prep, 
                                        n_perms = 100, seed = 123,
                                        use_stan = TRUE, 
                                        chains = 4, iter = 2000,
                                        verbose = TRUE) {
  if (verbose) {
    message("\n========================================")
    message("CORRESPONDENCE PERMUTATION TEST")
    message("========================================\n")
  }
  
  if (verbose) {
    message("Calculando baseline con datos originales...")
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
    message(sprintf("Acoplamiento baseline: %.4f", baseline_coupling))
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
  
  for (i in 1:n_perms) {
    perm_idx <- sample(ncol(factors_data$scores_Y))
    factors_perm <- factors_data
    factors_perm$scores_Y <- factors_data$scores_Y[, perm_idx]
    
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
    message(sprintf("\n\nPermutaciones exitosas: %d/%d", successful_perms, n_perms))
  }
  
  if (length(valid_perms) < 10) {
    if (verbose) message("FAIL: Insufficient successful permutations")
    return(NULL)
  }
  
  mean_perm <- mean(valid_perms)
  sd_perm <- sd(valid_perms)
  median_perm <- median(valid_perms)
  
  if (verbose) {
    message("\n=====================================")
    message("RESULTS:")
    message(sprintf("  Original coupling: %.4f", baseline_coupling))
    message(sprintf("  Permuted coupling:"))
    message(sprintf("    - Mean: %.4f (SD: %.4f)", mean_perm, sd_perm))
    message(sprintf("    - Median: %.4f", median_perm))
  }
  
  p_value <- mean(valid_perms >= baseline_coupling)
  z_score <- (baseline_coupling - mean_perm) / (sd_perm + 1e-10)
  
  if (verbose) {
    message(sprintf("  p-valor: %.4f", p_value))
    message(sprintf("  z-score: %.2f", z_score))
  }
  
  if (z_score > 2) {
    if (verbose) message("\nOK: ROBUST signal (z > 2)")
    robust <- TRUE
  } else if (z_score > 1) {
    if (verbose) message("\nWARNING: Moderate signal (1 < z < 2)")
    robust <- FALSE
  } else {
    if (verbose) message("\nFAIL: Spurious signal (z < 1)")
    robust <- FALSE
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
#' Tests sensitivity of convergence results to alternative sectoral weighting
#' schemes by re-running Bayesian disaggregation and full pipeline with
#' perturbed weights.
#'
#'
#' @return List with components:
#'   \describe{
#'     \item{\code{lambda_distribution}}{Matrix of lambda estimates across schemes.}
#'     \item{\code{original_lambda}}{Baseline lambda from original weights.}
#'     \item{\code{robust_factors}}{Logical vector indicating robust convergence.}
#'     \item{\code{sensitivity_metrics}}{Summary statistics of sensitivity.}
#'   }
#'
#'
#' @param path_cpi Path to CPI data file
#' @param path_weights Path to weights data file
#' @param X_matrix Matrix of variables
#' @param max_comp Maximum number of components (default: 3)
#' @param verbose Logical; print progress and diagnostic information. Default \code{TRUE}.
#' @export

test_reweighting_robustness <- function(path_cpi, path_weights, X_matrix,
                                        max_comp = 3, verbose = TRUE) {
  if (verbose) {
    message("\n========================================")
    message("ROBUSTNESS TEST BY REWEIGHTING")
    message("========================================\n")
  }
  
  if (!file.exists(path_cpi) || !file.exists(path_weights)) {
    if (verbose) message("FAIL: Files not found")
    return(NULL)
  }
  
  tryCatch({
    weights_data <- read_weights_matrix(path_weights)
  }, error = function(e) {
    if (verbose) message("FAIL: Error reading weights:", e$message)
    return(NULL)
  })
  
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
  
  alt_priors <- generate_alternative_priors(mip_weights, n_sectors, verbose = verbose)
  
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

#' Jackknife robustness test by sector
#'
#' Assesses the influence of individual sectors by systematically dropping each
#' sector and re-estimating convergence parameters. Identifies influential
#' sectors and checks stability.
#'
#'   column names.
#'
#' @return List with components:
#'   \describe{
#'     \item{\code{jackknife_estimates}}{Matrix of lambda estimates (iterations x factors).}
#'     \item{\code{original_estimate}}{Original lambda from full sample.}
#'     \item{\code{bias}}{Estimated jackknife bias.}
#'     \item{\code{se}}{Jackknife standard errors.}
#'     \item{\code{influential_sectors}}{Character vector of highly influential sectors.}
#'     \item{\code{dfbetas}}{Matrix of influence measures (DFBETAS).}
#'   }
#'
#'
#' @param X_matrix Matrix of first set of variables
#' @param Y_matrix Matrix of second set of variables
#' @param sector_names Vector of sector names (default: NULL)
#' @param k_exclude Number of sectors to exclude in each iteration (default: 3)
#' @param max_comp Maximum number of components (default: 3)
#' @param verbose Logical; print progress and diagnostic information. Default \code{TRUE}.
#' @export

test_jackknife_sectors <- function(X_matrix, Y_matrix, sector_names = NULL,
                                   k_exclude = 3, max_comp = 3, verbose = TRUE) {
  if (verbose) {
    message("\n========================================")
    message("JACKKNIFE FOR HEAVY SECTORS")
    message("========================================\n")
  }
  
  var_sectors <- apply(Y_matrix, 2, var)
  top_k <- order(var_sectors, decreasing = TRUE)[1:min(k_exclude, ncol(Y_matrix))]
  
  if (verbose) {
    message("Excluyendo sectores:", paste(top_k, collapse = ", "), "\n")
  }
  
  data_full <- diagnose_data(X_matrix, Y_matrix, verbose = verbose)
  data_prep_full <- list(
    X_scaled = scale(data_full$X_matrix),
    Y_scaled = scale(data_full$Y_matrix)
  )
  
  selection_full <- select_optimal_components_safe(
    data_prep_full$X_scaled,
    data_prep_full$Y_scaled,
    max_comp,
    verbose = verbose
  )
  
  factors_full <- list(
    scores_X = selection_full$pls_model$scores[, 1:selection_full$optimal_ncomp, drop = FALSE],
    scores_Y = selection_full$pls_model$Yscores[, 1:selection_full$optimal_ncomp, drop = FALSE]
  )
  
  ou_full <- estimate_factor_OU_fallback(factors_full, verbose = verbose)
  baseline_coupling <- ou_full$coupling_strength
  
  Y_reduced <- Y_matrix[, -top_k, drop = FALSE]
  X_reduced <- X_matrix[, -top_k, drop = FALSE]
  
  data_red <- diagnose_data(X_reduced, Y_reduced, verbose = verbose)
  data_prep_red <- list(
    X_scaled = scale(data_red$X_matrix),
    Y_scaled = scale(data_red$Y_matrix)
  )
  
  selection_red <- select_optimal_components_safe(
    data_prep_red$X_scaled,
    data_prep_red$Y_scaled,
    max_comp,
    verbose = verbose
  )
  
  factors_red <- list(
    scores_X = selection_red$pls_model$scores[, 1:selection_red$optimal_ncomp, drop = FALSE],
    scores_Y = selection_red$pls_model$Yscores[, 1:selection_red$optimal_ncomp, drop = FALSE]
  )
  
  ou_red <- estimate_factor_OU_fallback(factors_red, verbose = verbose)
  jackknife_coupling <- ou_red$coupling_strength
  
  retention_rate <- jackknife_coupling / baseline_coupling
  
  if (verbose) {
    message(sprintf("Complete coupling: %.4f", baseline_coupling))
    message(sprintf("Reduced coupling: %.4f", jackknife_coupling))
    message(sprintf("Retention: %.0f%%", retention_rate * 100))
  }
  
  list(
    baseline = baseline_coupling,
    jackknife = jackknife_coupling,
    retention_rate = retention_rate,
    robust = retention_rate > 0.7
  )
}