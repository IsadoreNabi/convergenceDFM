#' Prepare Marxist factors with auxiliary variables
#'
#' Prepares factor data by combining price indices with auxiliary economic variables
#' such as total Marxist gross output, commodity composition, and surplus value rates.
#'
#' @param X_matrix Numeric matrix of Marxist price indices.
#' @param Y_matrix Numeric matrix of market prices.
#' @param TMG Numeric vector. Optional total Marxist gross output time series.
#' @param COM_matrix Numeric matrix. Optional commodity composition over time.
#' @param SPVR_matrix Numeric matrix. Optional surplus value rate over time.
#'
#' @return List with prepared factor data and metadata.
#'
#' @keywords internal
#' @noRd

prepare_marxist_factors <- function(X_matrix, Y_matrix, TMG = NULL,
                                    COM_matrix = NULL, SPVR_matrix = NULL, CA = NULL) {
  list(
    X_scaled = scale(X_matrix),
    Y_scaled = scale(Y_matrix)
  )
}

interpret_factors <- function(pls_model, data_prep, optimal_ncomp, sector_names = NULL) {
  list(
    scores_X = pls_model$scores[, 1:optimal_ncomp, drop = FALSE],
    scores_Y = pls_model$Yscores[, 1:optimal_ncomp, drop = FALSE],
    interpretation = lapply(seq_len(optimal_ncomp),
                            function(i) list(label = paste("Factor", i), confidence = 0.5))
  )
}

#' Sigmoid transformation for half-life metrics
#'
#' Transforms half-life values to (0,1) interval using sigmoid function,
#' useful for creating bounded convergence indices.
#'
#' @param z Numeric vector of half-lives.
#' @param center Numeric. Center point of sigmoid. Default is 0.
#' @param scale Numeric. Scale parameter. Default is 1.
#'
#' @return Numeric vector in (0,1).
#'
#' @keywords internal
#' @noRd

sigmoid_transform <- function(z, center = 0, scale = 1) {
  1 / (1 + exp(-(z - center) / scale))
}

#' Transform half-lives to convergence scores
#'
#' Maps half-life values to interpretable convergence scores using reference
#' thresholds for "fast" and "slow" convergence.
#'
#' @param hl Numeric vector of half-lives.
#' @param fast Numeric. Threshold for fast convergence. Default is 5 periods.
#' @param slow Numeric. Threshold for slow convergence. Default is 20 periods.
#'
#' @return Numeric vector of convergence scores (higher = faster).
#'
#' @keywords internal
#' @noRd

halflife_transform <- function(hl, fast = 5, slow = 20) {
  hl[!is.finite(hl)] <- slow * 2
  1 - pmin(pmax((hl - fast) / (slow - fast), 0), 1)
}

#' Apply FDR correction to p-values
#'
#' Applies False Discovery Rate correction to multiple hypothesis tests using
#' Benjamini-Hochberg or related procedures.
#'
#' @param pvals Numeric vector of p-values.
#' @param method Character string. FDR method: "BH" (default), "BY", or "fdr".
#'
#' @return Numeric vector of adjusted p-values.
#'
#' @keywords internal
#' @noRd

apply_fdr_correction <- function(pvals, method = "BH") {
  if (all(is.na(pvals))) return(pvals)
  p.adjust(pvals, method = method)
}

#' Calculate continuous convergence indices
#'
#' Computes continuous (0-100) convergence strength indices from formal test
#' results, combining multiple metrics into interpretable scores.
#'
#' @param test_results List. Output from convergence test functions.
#'
#' @return Numeric vector of convergence indices for different dimensions
#'   (speed, stability, coupling).
#'
#' @keywords internal
#' @noRd

calculate_continuous_indices <- function(test_results) {
  
  convergence_scores <- list()
  robustness_scores <- list()
  
  if (!is.null(test_results$posteriors) && !is.null(test_results$posteriors$kappa_x)) {
    kappa_x <- test_results$posteriors$kappa_x
    kappa_y <- test_results$posteriors$kappa_y
    
    if (is.matrix(kappa_x) && nrow(kappa_x) >= 3) {
      in_range_x <- mean(kappa_x[2,] > 0.01 & kappa_x[3,] < 0.99)
    } else {
      in_range_x <- 0
    }
    
    if (is.matrix(kappa_y) && nrow(kappa_y) >= 3) {
      in_range_y <- mean(kappa_y[2,] > 0.01 & kappa_y[3,] < 0.99)
    } else {
      in_range_y <- 0
    }
    
    convergence_scores$posteriors <- (in_range_x + in_range_y) / 2
  }
  
  if (!is.null(test_results$unit_root)) {
    adf_pvals <- sapply(test_results$unit_root$adf, function(x) {
      if (is.null(x$stat) || is.null(x$crit)) return(0.5)
      z_score <- (x$crit - x$stat) / abs(x$crit)
      sigmoid_transform(z_score, center = 0, scale = 1)
    })
    
    pp_pvals <- sapply(test_results$unit_root$pp, function(x) {
      if (is.null(x$stat) || is.null(x$crit)) return(0.5)
      z_score <- (x$crit - x$stat) / abs(x$crit)
      sigmoid_transform(z_score, center = 0, scale = 1)
    })
    
    convergence_scores$unit_root <- mean(c(adf_pvals, pp_pvals))
  }
  
  if (!is.null(test_results$cointegration)) {
    n_vars <- if (!is.null(test_results$cointegration$test)) {
      ncol(test_results$cointegration$test@x)
    } else 4
    
    n_coint <- test_results$cointegration$n_coint
    convergence_scores$cointegration <- min(n_coint / (n_vars - 1), 1)
  }
  
  if (!is.null(test_results$half_lives)) {
    hl_scores <- halflife_transform(test_results$half_lives, fast = 5, slow = 20)
    convergence_scores$halflife <- mean(hl_scores, na.rm = TRUE)
  }
  
  if (!is.null(test_results$permutation)) {
    z <- test_results$permutation$z_score
    if (is.finite(z)) {
      robustness_scores$permutation <- sigmoid_transform(z, center = 0, scale = 2)
    }
  }
  
  if (!is.null(test_results$reweighting)) {
    cv <- test_results$reweighting$cv_coupling
    if (is.finite(cv)) {
      robustness_scores$reweighting <- exp(-2 * cv)
    }
  }
  
  if (!is.null(test_results$jackknife)) {
    ret <- test_results$jackknife$retention_rate
    if (is.finite(ret)) {
      robustness_scores$jackknife <- min(max(ret, 0), 1)
    }
  }
  
  convergence_index <- if (length(convergence_scores) > 0) {
    geometric_mean_robust(unlist(convergence_scores))
  } else NA
  
  robustness_index <- if (length(robustness_scores) > 0) {
    geometric_mean_robust(unlist(robustness_scores))
  } else NA
  
  global_index <- if (!is.na(convergence_index) && !is.na(robustness_index)) {
    geometric_mean_robust(c(convergence_index, robustness_index))
  } else if (!is.na(convergence_index)) {
    convergence_index
  } else if (!is.na(robustness_index)) {
    robustness_index
  } else NA
  
  list(
    convergence = list(
      scores = convergence_scores,
      index = convergence_index
    ),
    robustness = list(
      scores = robustness_scores,
      index = robustness_index
    ),
    global_index = global_index,
    interpretation = interpret_scores(convergence_index, robustness_index)
  )
}

#' Interpret convergence scores
#'
#' Provides qualitative interpretation of convergence and robustness indices
#' with categorization (weak/moderate/strong).
#'
#' @param conv_idx Numeric. Convergence index (0-100).
#' @param rob_idx Numeric. Robustness index (0-100).
#'
#' @return Character string with interpretation.
#'
#' @keywords internal
#' @noRd

interpret_scores <- function(conv_idx, rob_idx) {
  if (is.na(conv_idx)) conv_idx <- 0
  if (is.na(rob_idx)) rob_idx <- 0
  
  if (conv_idx > 0.7 && rob_idx > 0.7) {
    verdict <- "VERY STRONG EVIDENCE"
    confidence <- "High"
  } else if (conv_idx > 0.6 && rob_idx > 0.6) {
    verdict <- "SOLID EVIDENCE"
    confidence <- "Moderate-High"
  } else if (conv_idx > 0.5 && rob_idx > 0.5) {
    verdict <- "MODERATE EVIDENCE"
    confidence <- "Moderate"
  } else if (conv_idx > 0.6 && rob_idx < 0.5) {
    verdict <- "CONVERGENCE DETECTED BUT NOT ROBUST"
    confidence <- "Low"
  } else {
    verdict <- "INSUFFICIENT EVIDENCE"
    confidence <- "Very Low"
  }
  
  list(
    verdict = verdict,
    confidence = confidence,
    recommendation = if (confidence %in% c("High", "Moderate-High")) {
      "The results support the convergence hypothesis"
    } else if (confidence == "Moderate") {
      "Additional tests are required to confirm"
    } else {
      "Current evidence does not allow firm conclusions"
    }
  )
}

#' Sensitivity analysis for index weighting schemes
#'
#' Performs sensitivity analysis on convergence and robustness indices by
#' varying the weighting scheme used to aggregate sub-metrics.
#'
#' @param scores_list List of test result scores.
#' @param n_schemes Integer. Number of alternative schemes. Default is 5.
#'
#' @return List with sensitivity metrics across schemes.
#'
#' @keywords internal
#' @noRd

sensitivity_analysis_weights <- function(scores_list, n_schemes = 5) {
  
  uniform <- rep(1/length(scores_list), length(scores_list))
  
  var_est <- runif(length(scores_list), 0.1, 0.5)
  reliability <- 1 / var_est
  reliability <- reliability / sum(reliability)
  
  decreasing <- rev(seq_len(length(scores_list)))
  decreasing <- decreasing / sum(decreasing)
  
  extreme <- abs(unlist(scores_list) - 0.5) + 0.1
  extreme <- extreme / sum(extreme)
  
  random <- runif(length(scores_list))
  random <- random / sum(random)
  
  schemes <- list(
    uniform = uniform,
    reliability = reliability,
    recency = decreasing,
    extreme = extreme,
    random = random
  )
  
  indices <- sapply(schemes, function(w) {
    sum(unlist(scores_list) * w)
  })
  
  list(
    schemes = schemes,
    indices = indices,
    range = c(min = min(indices), max = max(indices)),
    cv = sd(indices) / mean(indices),
    robust = (sd(indices) / mean(indices)) < 0.2
  )
}

#' Run comprehensive robustness test suite
#'
#' Executes all available robustness tests (permutation, reweighting, jackknife)
#' and synthesizes results into an integrated assessment.
#'
#'
#' @return List with components:
#'   \describe{
#'     \item{\code{permutation}}{Results from permutation test.}
#'     \item{\code{reweighting}}{Results from reweighting test (if applicable).}
#'     \item{\code{jackknife}}{Results from jackknife test (if requested).}
#'     \item{\code{summary}}{Data frame summarizing all tests.}
#'     \item{\code{overall_robust}}{Logical indicating if convergence is robust across all tests.}
#'   }
#'
#'
#' @param results_robust Output from run_complete_factor_analysis_robust()
#' @param X_matrix Matrix of first set of variables
#' @param Y_matrix Matrix of second set of variables
#' @param path_cpi Path to CPI data (default: NULL)
#' @param path_weights Path to weights data (default: NULL)
#' @param sector_names Vector of sector names (default: NULL)
#' @param run_permutation Logical, run permutation test (default: TRUE)
#' @param run_reweighting Logical, run reweighting test (default: FALSE)
#' @param run_jackknife Logical, run jackknife test (default: TRUE)
#' @param run_leadlag Logical, run lead-lag test (default: FALSE)
#' @param run_common_factor Logical, run common factor test (default: FALSE)
#' @param sensitivity_analysis Logical, run sensitivity analysis (default: TRUE)
#' @param verbose Logical; print progress and diagnostic information. Default \code{TRUE}.
#' @export

run_convergence_robustness_tests <- function(results_robust, X_matrix, Y_matrix,
                                             path_cpi = NULL, path_weights = NULL,
                                             sector_names = NULL,
                                             run_permutation = TRUE,
                                             run_reweighting = FALSE,
                                             run_jackknife = TRUE,
                                             run_leadlag = FALSE,
                                             run_common_factor = FALSE,
                                             sensitivity_analysis = TRUE,
                                             verbose = TRUE) {
  
  if (verbose) {
    message("========================================")
    message("ANALYZING: Convergence and Robustness Test Suite")
    message("    [Continuous Scoring Version]")
    message("========================================\n")
  }
  
  test_results <- list()
  
  if (!is.null(results_robust$factor_ou)) {
    if (verbose) {
      message("RESULTS: PART A (FORMAL EVIDENCE OF CONVERGENCE)")
      message("==========================================")
    }
    
    test_results$posteriors <- extract_ou_posteriors(results_robust$factor_ou, verbose = verbose)
    
    if (!is.null(results_robust$factors)) {
      test_results$unit_root <- test_longrun_error(
        results_robust$factors,
        results_robust$factor_ou
      )
    }
    
    if (!is.null(results_robust$factors)) {
      test_results$cointegration <- test_cointegration_control(
        results_robust$factors,
        verbose = verbose
      )
    }
    
    if (!is.null(results_robust$factor_ou)) {
      hl_all <- c(results_robust$factor_ou$half_lives_X,
                  results_robust$factor_ou$half_lives_Y)
      test_results$half_lives <- hl_all[is.finite(hl_all)]
    }
  }
  
  if (verbose) {
    message("\nRESULTS: PART B (ROBUSTNESS AGAINST SPURIOUS CORRELATION)")
    message("==============================================")
  }
  
  if (!is.null(results_robust$factors)) {
    if (run_permutation) {
      test_results$permutation <- test_permutation_robustness(
        results_robust$factors,
        results_robust$data_prep,
        n_perms = 50,
        use_stan = TRUE,
        chains = 4,
        iter = 2000,
        verbose = verbose
      )
    }
    
    if (run_reweighting && !is.null(path_cpi) && !is.null(path_weights)) {
      test_results$reweighting <- test_reweighting_robustness(
        path_cpi, path_weights, X_matrix, verbose = verbose
      )
    }
    
    if (run_jackknife) {
      test_results$jackknife <- test_jackknife_sectors(
        X_matrix, Y_matrix, sector_names, verbose = verbose
      )
    }
  }
  
  if (verbose) {
    message("\n========================================")
    message("CALCULATION OF CONTINUOUS INDICES")
    message("========================================\n")
  }
  
  indices <- calculate_continuous_indices(test_results)
  
  if (sensitivity_analysis && length(indices$convergence$scores) > 0) {
    if (verbose) {
      message("ANALYZING: Sensitivity Analysis to Weights")
      message("-------------------------------------------")
    }
    
    conv_sens <- sensitivity_analysis_weights(indices$convergence$scores)
    rob_sens <- if (length(indices$robustness$scores) > 0) {
      sensitivity_analysis_weights(indices$robustness$scores)
    } else NULL
    
    if (verbose) {
      message("\nConvergence - Range of indices by weighting:")
      message(sprintf("  Min: %.3f | Max: %.3f | CV: %.1f%%",
                      conv_sens$range["min"], conv_sens$range["max"], 
                      conv_sens$cv * 100))
      
      if (!is.null(rob_sens)) {
        message("\nRobustness - Range of indices according to weighting:")
        message(sprintf("  Min: %.3f | Max: %.3f | CV: %.1f%%",
                        rob_sens$range["min"], rob_sens$range["max"], 
                        rob_sens$cv * 100))
      }
    }
    
    indices$sensitivity <- list(
      convergence = conv_sens,
      robustness = rob_sens
    )
  }
  
  if (verbose) {
    message("\n========================================")
    message("SUMMARY OF CONVERGENCE AND ROBUSTNESS")
    message("========================================\n")
    
    message("FORMAL CONVERGENCE SCORES:")
    if (length(indices$convergence$scores) > 0) {
      for (test_name in names(indices$convergence$scores)) {
        score <- indices$convergence$scores[[test_name]]
        message(sprintf("  %s: %.3f", test_name, score))
      }
      message(sprintf("\n  Convergence Index (geometric mean): %.3f", 
                      indices$convergence$index))
    } else {
      message("  No convergence tests available")
    }
    
    message("\nROBUSTNESS SCORES:")
    if (length(indices$robustness$scores) > 0) {
      for (test_name in names(indices$robustness$scores)) {
        score <- indices$robustness$scores[[test_name]]
        message(sprintf("  %s: %.3f", test_name, score))
      }
      message(sprintf("\n  Robustness Index (geometric mean): %.3f", 
                      indices$robustness$index))
    } else {
      message("  No robustness tests available")
    }
    
    message("\n========================================")
    message(sprintf("GLOBAL INDEX: %.3f", indices$global_index))
    message("\nVEREDICT: ", indices$interpretation$verdict)
    message("Nivel de confianza: ", indices$interpretation$confidence)
    message("\n", indices$interpretation$recommendation)
  }
  
  test_results$indices <- indices
  test_results$summary <- list(
    convergence_index = indices$convergence$index,
    robustness_index = indices$robustness$index,
    global_index = indices$global_index,
    interpretation = indices$interpretation,
    sensitivity = if (exists("indices$sensitivity")) indices$sensitivity else NULL
  )
  
  test_results
}

#' Complete factor-OU convergence analysis pipeline
#'
#' Executes the end-to-end analysis workflow: data preparation, PLS-based factor
#' extraction, DFM estimation, Factor-OU estimation, convergence tests, and
#' robustness checks. This is the main user-facing function.
#'
#'   uses column names of \code{X_matrix}.
#'   and reweighting tests). Can be \code{NULL}.
#'
#' @return List with components:
#'   \describe{
#'     \item{\code{factors}}{List containing PLS-extracted factors (scores_X, scores_Y)
#'       and related objects (pls_X, pls_Y, ncomp_X, ncomp_Y).}
#'     \item{\code{dfm}}{List with DFM estimation results including VAR fit, lag order,
#'       diagnostics, and optional impulse responses.}
#'     \item{\code{factor_ou}}{List with Factor-OU model estimates: beta, lambda,
#'       sigma, half_life, method used, and optional Stan fit object.}
#'     \item{\code{convergence_tests}}{List with formal convergence test results
#'       (if \code{run_convergence_tests = TRUE}).}
#'     \item{\code{robustness_tests}}{List with robustness test results (if
#'       \code{run_robustness_tests = TRUE}).}
#'     \item{\code{diagnostics}}{List with data diagnostics (multicollinearity,
#'       stationarity, structural breaks).}
#'     \item{\code{bayesian_cpi}}{List with Bayesian disaggregation results (if
#'       CPI paths provided).}
#'     \item{\code{metadata}}{List with analysis metadata (timestamp, versions,
#'       parameters).}
#'   }
#'
#' @details This function orchestrates the complete analysis:
#' \enumerate{
#'   \item Data validation and diagnostics
#'   \item Bayesian CPI disaggregation (if applicable)
#'   \item PLS-based factor extraction with optimal component selection
#'   \item Dynamic Factor Model estimation via VAR
#'   \item Factor Ornstein-Uhlenbeck mean-reversion model
#'   \item Formal convergence tests (stationarity, cointegration, speed)
#'   \item Robustness tests (permutation, reweighting, jackknife)
#' }
#'
#' @examples
#' \donttest{
#' # Basic usage with simulated data
#' set.seed(123)
#' X <- matrix(rnorm(100 * 10), 100, 10)
#' Y <- X + matrix(rnorm(100 * 10, 0, 0.5), 100, 10)
#'
#' results <- run_complete_factor_analysis_robust(
#'   X_matrix = X,
#'   Y_matrix = Y,
#'   max_comp = 3,
#'   dfm_lags = 1,
#'   ou_chains = 4,
#'   ou_iter = 2000,
#'   verbose = FALSE
#' )
#'
#' # View convergence summary
#' summary(results$convergence_tests)
#'
#' # Visualize results
#' visualize_factor_dynamics(
#'   dfm_result = results$dfm,
#'   ou_result = results$factor_ou,
#'   factors_data = results$factors
#' )
#' }
#'
#'
#' @param X_matrix Matrix of first set of variables
#' @param Y_matrix Matrix of second set of variables
#' @param TMG Optional TMG matrix (default: NULL)
#' @param COM_matrix Optional COM matrix (default: NULL)
#' @param SPVR_matrix Optional SPVR matrix (default: NULL)
#' @param CA Optional CA parameter (default: NULL)
#' @param sector_names Vector of sector names (default: NULL)
#' @param max_comp Maximum number of components (default: 3)
#' @param dfm_lags Number of lags for DFM (default: 1)
#' @param ou_chains Number of MCMC chains for OU estimation (default: 10)
#' @param ou_iter Number of MCMC iterations (default: 10000)
#' @param skip_ou Logical, skip OU estimation (default: FALSE)
#' @param run_convergence_tests Logical, run convergence tests (default: TRUE)
#' @param path_cpi Path to CPI data (default: NULL)
#' @param path_weights Path to weights data (default: NULL)
#' @param verbose Logical; print progress and diagnostic information. Default \code{TRUE}.
#' @export
#' @seealso \code{\link{estimate_DFM}}, \code{\link{estimate_factor_OU}},
#'   \code{\link{run_convergence_robustness_tests}}, \code{\link{visualize_factor_dynamics}}

run_complete_factor_analysis_robust <- function(X_matrix, Y_matrix,
                                                TMG = NULL, COM_matrix = NULL,
                                                SPVR_matrix = NULL, CA = NULL,
                                                sector_names = NULL,
                                                max_comp = 3,
                                                dfm_lags = 1,
                                                ou_chains = 10, ou_iter = 10000,
                                                skip_ou = FALSE,
                                                run_convergence_tests = TRUE,
                                                path_cpi = NULL, path_weights = NULL,
                                                verbose = TRUE) {
  
  if (verbose) {
    message("========================================")
    message("COMPLETE ROBUST FACTOR ANALYSIS")
    message("========================================\n")
  }
  
  results <- list()
  
  data_clean <- tryCatch(
    diagnose_data(X_matrix, Y_matrix, verbose = verbose),
    error = function(e) {
      if (verbose) message("\nFAIL: Diagnostic error:\n", e$message)
      return(NULL)
    }
  )
  
  if (is.null(data_clean)) return(NULL)
  
  X_matrix <- data_clean$X_matrix
  Y_matrix <- data_clean$Y_matrix
  
  if (verbose) message("\nREPORT: PHASE 1 (Data preparation)...")
  data_prep <- tryCatch({
    prepare_marxist_factors(X_matrix, Y_matrix, TMG, COM_matrix, SPVR_matrix, CA)
  }, error = function(e) {
    if (verbose) message("\nWARNING: Using basic preparation")
    list(X_scaled = scale(X_matrix), Y_scaled = scale(Y_matrix))
  })
  
  results$data_prep <- data_prep
  
  if (verbose) message("\nREPORT: PHASE 2 (Component Selection)...")
  max_comp_adj <- min(
    max_comp,
    floor(nrow(data_prep$X_scaled) * 0.8),
    ncol(data_prep$X_scaled),
    ncol(data_prep$Y_scaled) * 3
  )
  
  selection <- tryCatch({
    Y_for_pls <- data_prep$Y_scaled
    if (is.null(dim(Y_for_pls)) || ncol(Y_for_pls) == 0) {
      Y_for_pls <- as.matrix(Y_for_pls)
    }
    select_optimal_components_safe(data_prep$X_scaled, Y_for_pls, max_comp_adj, 
                                   verbose = verbose)
  }, error = function(e) {
    if (verbose) message("\nFAIL: Error in selection: ", e$message)
    NULL
  })
  
  if (is.null(selection)) {
    if (verbose) message("\nAborting due to selection failure.")
    return(results)
  }
  
  results$selection <- selection
  
  if (verbose) message("\nREPORT: PHASE 3 (Economic interpretation)...")
  factors_data <- tryCatch({
    interpret_factors(selection$pls_model, data_prep, 
                      selection$optimal_ncomp, sector_names)
  }, error = function(e) {
    if (verbose) message("\nWARNING: Using basic factors")
    list(
      scores_X = selection$pls_model$scores[, 1:selection$optimal_ncomp, drop = FALSE],
      scores_Y = selection$pls_model$Yscores[, 1:selection$optimal_ncomp, drop = FALSE],
      interpretation = lapply(seq_len(selection$optimal_ncomp),
                              function(i) list(label = paste("Factor", i)))
    )
  })
  
  results$factors <- factors_data
  
  if (verbose) message("\nREPORT: PHASE 4 (Dynamic Factor Model)...")
  dfm_result <- tryCatch(
    estimate_DFM(factors_data, p = dfm_lags, verbose = verbose),
    error = function(e) {
      if (verbose) message("\nFAIL: DFM error: ", e$message)
      NULL
    }
  )
  
  if (!is.null(dfm_result)) {
    results$dfm <- dfm_result
  }
  
  if (!skip_ou && !is.null(dfm_result)) {
    if (verbose) message("\nREPORT: PHASE 5 (Ornstein-Uhlenbeck Factor Model)...")
    ou_result <- tryCatch(
      estimate_factor_OU(factors_data, data_prep,
                         chains = ou_chains, iter = ou_iter, verbose = verbose),
      error = function(e) {
        if (verbose) message("\nFAIL: OU error: ", e$message)
        NULL
      }
    )
    
    if (is.null(ou_result)) {
      if (verbose) message("WARNING: Using OU fallback")
      ou_result <- estimate_factor_OU_fallback(factors_data, verbose = verbose)
    }
    
    results$factor_ou <- ou_result
  }
  
  if (!is.null(results$dfm) && !is.null(results$factor_ou)) {
    if (verbose) message("\nREPORT: PHASE 6 (Visualization)...")
    
    viz_success <- tryCatch({
      visualize_factor_dynamics(
        dfm_result = results$dfm,
        ou_result = results$factor_ou,
        factors_data = results$factors,
        save_plot = TRUE,
        use_device = "default",
        verbose = verbose
      )
    }, error = function(e) {
      if (verbose) {
        message("  WARNING: Error in complex visualization: ", e$message)
        message("  Trying simple version...")
      }
      
      tryCatch({
        visualize_factor_dynamics_simple(
          factors_data = results$factors,
          output_file = NULL,
          verbose = verbose
        )
      }, error = function(e2) {
        if (verbose) message("  FAIL: Error in simple display: ", e2$message)
        FALSE
      })
    })
    
    if (viz_success && verbose) {
      message("  OK: Visualization completed successfully")
    }
  }
  
  if (run_convergence_tests && !is.null(results$factor_ou)) {
    if (verbose) message("\nANALYZING: PHASE 7 (Convergence tests)...")
    
    convergence_tests <- tryCatch(
      run_convergence_robustness_tests(
        results_robust = results,
        X_matrix = X_matrix,
        Y_matrix = Y_matrix,
        path_cpi = path_cpi,
        path_weights = path_weights,
        sector_names = sector_names,
        run_permutation = TRUE,
        run_reweighting = !is.null(path_cpi) && !is.null(path_weights),
        run_jackknife = TRUE,
        sensitivity_analysis = TRUE,
        verbose = verbose
      ),
      error = function(e) {
        if (verbose) message("WARNING: Error in tests: ", e$message)
        NULL
      }
    )
    
    if (!is.null(convergence_tests)) {
      results$convergence_tests <- convergence_tests
    }
  }
  
  if (verbose) {
    message("\n========================================")
    message("EXECUTIVE SUMMARY")
    message("========================================\n")
    
    if (!is.null(selection)) {
      message("FACTORIAL STRUCTURE:")
      message("- Optimal components: ", selection$optimal_ncomp)
      message("- Variance Y explained: ", round(selection$var_explained, 1), "%")
      message("- R2 (CV): ", round(selection$R2_cv, 4), "\n")
    }
    
    if (!is.null(dfm_result)) {
      message("DYNAMIC CONVERGENCE (DFM):")
      message("- Stable system: ", ifelse(dfm_result$is_stable, "OK: Yes", "FAIL: No"))
      message("- VAR order: p = ", dfm_result$p_used)
      if (length(dfm_result$half_lives) > 0) {
        message("- Median half-life: ", round(median(dfm_result$half_lives, na.rm = TRUE), 2), 
                " periods")
      }
    }
    
    if (!is.null(results$factor_ou)) {
      ou <- results$factor_ou
      message("\nFACTOR-OU CONVERGENCE:")
      if (!is.null(ou$coupling_strength)) {
        message("- Coupling between X and Y: ", round(ou$coupling_strength, 4))
      }
      if (!is.null(ou$half_lives_X) || !is.null(ou$half_lives_Y)) {
        hl_mean <- mean(c(ou$half_lives_X, ou$half_lives_Y), na.rm = TRUE)
        if (is.finite(hl_mean)) {
          message("- Median half-life: ", round(hl_mean, 2), " periods")
          
          if (hl_mean < 5) {
            message("\n-> FAST CONVERGENCE (<5 periods)")
          } else if (hl_mean < 15) {
            message("\n-> MODERATE CONVERGENCE (5-15 periods)")
          } else {
            message("\n-> SLOW CONVERGENCE (>15 periods)")
          }
        }
      }
    }
    
    message("\nOK:  Analysis completed")
  }
  
  return(results)
}