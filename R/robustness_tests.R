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
#' schemes. For each perturbed weight vector the sectoral price levels are
#' obtained from the canonical closed-form conjugate disaggregation
#' ([BayesianDisaggregation::disaggregate_conjugate()], a Kalman/RTS smoother that
#' conditions genuinely on the aggregate CPI), and the (fallback) factor-OU
#' coupling is re-estimated. The coefficient of variation of the coupling across
#' schemes summarizes the sensitivity.
#'
#' @section Canonical disaggregation:
#' Earlier versions used a deterministic weight blend
#' (`run_disaggregation_custom_prior`) that did not condition on the CPI. That
#' method is removed; the single source of truth is now the conjugate engine in
#' \pkg{BayesianDisaggregation}. A constant-in-time custom prior is replicated
#' across periods to form the value-added (VAB) weight matrix `W`, and the point
#' estimate used is the posterior median (`$phi_summary$median`), i.e. the
#' RTS-smoothed sectoral levels, now genuinely conditioned on the CPI.
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
  
  cpi_df <- tryCatch(
    read_cpi(path_cpi),
    error = function(e) {
      if (verbose) message("FAIL: Error reading CPI: ", e$message)
      NULL
    }
  )
  if (is.null(cpi_df)) {
    if (verbose) message("FAIL: The CPI could not be read")
    return(NULL)
  }

  # Align the CPI to the years covered by the weight matrix. The conjugate
  # disaggregation engine consumes a numeric CPI vector and a [T x K] weight
  # matrix (not file paths), so the alignment is computed once here and reused
  # across the perturbed-prior loop.
  common_years <- sort(intersect(weights_data$years, cpi_df$Year))
  if (length(common_years) < 2) {
    if (verbose) message("FAIL: Fewer than 2 years in common between CPI and weights")
    return(NULL)
  }
  cpi_vec   <- cpi_df$CPI[match(common_years, cpi_df$Year)]
  Tn_disagg <- length(common_years)

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
      # A constant-in-time custom prior is replicated across periods to form the
      # VAB weight matrix; the canonical conjugate engine returns the smoothed
      # sectoral levels (posterior median) genuinely conditioned on the CPI.
      W_custom <- matrix(rep(prior_weights, each = Tn_disagg), nrow = Tn_disagg)
      Y_reweighted <- BayesianDisaggregation::disaggregate_conjugate(
        cpi = cpi_vec, W = W_custom,
        years = common_years, industries = weights_data$industries
      )$phi_summary$median

      min_rows <- min(nrow(Y_reweighted), nrow(X_matrix))
      min_cols <- min(ncol(Y_reweighted), ncol(X_matrix))
      
      Y_reweighted <- Y_reweighted[1:min_rows, 1:min_cols]
      X_matrix_adj <- X_matrix[1:min_rows, 1:min_cols]
      
      # Propagate `seed` to the diagnosis and PLS component selection so the
      # whole per-scheme coupling is reproducible. (Previously only the prior
      # generation was seeded; the cross-validated PLS folds drew on the ambient
      # RNG, making the couplings depend on call order -- a reproducibility gap.)
      data_clean <- diagnose_data(X_matrix_adj, Y_reweighted, verbose = verbose,
                                  seed = seed)

      data_prep <- list(
        X_scaled = scale(data_clean$X_matrix),
        Y_scaled = scale(data_clean$Y_matrix)
      )

      selection <- suppressWarnings(
        select_optimal_components_safe(
          data_prep$X_scaled,
          data_prep$Y_scaled,
          max_comp = max_comp,
          seed = seed,
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
      # Use <<- so the failure marker reaches the enclosing results list (the
      # previous `<-` wrote to the handler frame only and was lost).
      results_reweight[[prior_name]] <<- list(coupling = NA)
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

#' Re-estimate the X->Y coupling on a (possibly reduced) sector set
#'
#' Shared coupling pipeline used by the delete-one-sector jackknife and the
#' leave-cluster-out test: diagnose/clean, select PLS components, and read the
#' fallback factor-OU coupling strength. Returns \code{NA_real_} on any
#' degenerate or failing step (too few columns/rows, selection or OU failure) so
#' callers can drop the replicate cleanly. Extracted from a former in-function
#' closure so the two leave-out tests share one implementation.
#'
#' @param Xm,Ym Numeric matrices (sectors in columns).
#' @param max_comp Maximum number of PLS components.
#' @param seed Optional integer for reproducible component selection/diagnostics.
#' @return Numeric scalar coupling strength, or \code{NA_real_}.
#' @keywords internal
#' @noRd
.coupling_estimate <- function(Xm, Ym, max_comp = 3, seed = 123) {
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

  baseline <- .coupling_estimate(X_matrix, Y_matrix, max_comp = max_comp, seed = seed)
  if (!is.finite(baseline)) {
    if (verbose) message("FAIL: baseline coupling could not be computed")
    return(NULL)
  }

  jk <- rep(NA_real_, n_sec); names(jk) <- sector_names
  if (verbose) pb <- txtProgressBar(min = 0, max = n_sec, style = 3)
  for (j in seq_len(n_sec)) {
    Xj <- if (j <= ncol(X_matrix)) X_matrix[, -j, drop = FALSE] else X_matrix
    Yj <- Y_matrix[, -j, drop = FALSE]
    jk[j] <- .coupling_estimate(Xj, Yj, max_comp = max_comp, seed = seed)
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

#' Build a fallback sector-to-cluster map
#'
#' Constructs a sector-to-cluster partition for [test_leave_cluster_out()] when a
#' genuine input-output (Leontief, "MIP") partition is not yet available. This is
#' an explicitly documented **stopgap proxy**: the canonical clusters are value
#' chains defined by inter-industry demand linkages, which neither correlation
#' nor organic composition reproduces. Supply the real partition through the
#' `cluster_map` argument of [test_leave_cluster_out()] once the MIP is at hand.
#'
#' @section Fallback methods:
#' \describe{
#'   \item{`"correlation"`}{Average-linkage hierarchical clustering on the
#'     correlation distance \eqn{1 - \rho_{ij}} between the sectoral series
#'     (columns of `Y_matrix`), cut into `n_clusters` groups. Sectors that move
#'     together over time are grouped; this is a purely statistical proxy for
#'     co-movement, not a demand-linkage cluster.}
#'   \item{`"com"`}{Quantile binning of a per-sector organic-composition vector
#'     `com` into `n_clusters` groups. Motivated by the gravitation reading
#'     (sectors of similar organic composition share a profit-rate neighbourhood),
#'     but again a one-dimensional proxy, not an input-output partition.}
#' }
#'
#' @param Y_matrix Numeric matrix (periods in rows, sectors in columns).
#' @param n_clusters Integer number of clusters (clamped to `1..ncol(Y_matrix)`).
#'   Default `3`.
#' @param method Fallback rule, `"correlation"` (default) or `"com"`.
#' @param com Numeric vector of organic-composition values, one per sector
#'   (length `ncol(Y_matrix)`); required only for `method = "com"`.
#'
#' @return A named integer vector of length `ncol(Y_matrix)`: the cluster label of
#'   each sector, named by sector.
#'
#' @seealso [test_leave_cluster_out()]
#' @export
build_cluster_map <- function(Y_matrix, n_clusters = 3,
                              method = c("correlation", "com"), com = NULL) {
  method <- match.arg(method)
  Y <- as.matrix(Y_matrix)
  K <- ncol(Y)
  labels <- colnames(Y)
  if (is.null(labels)) labels <- paste0("S", seq_len(K))
  n_clusters <- max(1L, min(as.integer(n_clusters), K))

  if (method == "correlation") {
    if (K < 2L) {
      cl <- rep(1L, K)
    } else {
      R <- suppressWarnings(stats::cor(Y, use = "pairwise.complete.obs"))
      R[!is.finite(R)] <- 0
      d <- stats::as.dist(1 - R)
      hc <- stats::hclust(d, method = "average")
      cl <- as.integer(stats::cutree(hc, k = n_clusters))
    }
  } else {
    if (is.null(com)) {
      stop("method = 'com' requires `com` (one organic-composition value per sector).",
           call. = FALSE)
    }
    com <- as.numeric(com)
    if (length(com) != K) {
      stop("`com` must have one value per sector (length ncol(Y_matrix) = ", K, ").",
           call. = FALSE)
    }
    br <- stats::quantile(com, probs = seq(0, 1, length.out = n_clusters + 1),
                          na.rm = TRUE)
    br[1] <- -Inf; br[length(br)] <- Inf
    cl <- as.integer(cut(com, breaks = unique(br), include.lowest = TRUE))
  }
  names(cl) <- labels
  cl
}

#' Normalize a user-supplied cluster map to a per-sector label vector
#'
#' Accepts either a per-sector vector (one cluster label per sector, optionally
#' named) or a named list (cluster label -> sector names or column indices) and
#' returns a character vector of length `K` named by sector. Errors if the map
#' does not cover every sector exactly once.
#' @keywords internal
#' @noRd
.normalize_cluster_map <- function(cluster_map, sector_names) {
  K <- length(sector_names)
  if (is.list(cluster_map) && !is.data.frame(cluster_map)) {
    lab <- rep(NA_character_, K)
    names(lab) <- sector_names
    cl_names <- names(cluster_map)
    if (is.null(cl_names)) cl_names <- paste0("cluster_", seq_along(cluster_map))
    for (i in seq_along(cluster_map)) {
      members <- cluster_map[[i]]
      if (is.numeric(members)) members <- sector_names[members]
      members <- intersect(members, sector_names)
      lab[members] <- cl_names[i]
    }
    if (anyNA(lab)) {
      stop("cluster_map (list) does not cover every sector; unassigned: ",
           paste(sector_names[is.na(lab)], collapse = ", "), call. = FALSE)
    }
    return(lab)
  }
  v <- cluster_map
  if (length(v) != K) {
    stop("cluster_map must have one entry per sector (length ncol(Y_matrix) = ",
         K, ").", call. = FALSE)
  }
  if (anyNA(v)) stop("cluster_map has missing entries.", call. = FALSE)
  v <- as.character(v)
  names(v) <- sector_names
  v
}

#' Leave-cluster-out robustness of the X->Y coupling
#'
#' Generalizes the delete-one-sector jackknife ([test_jackknife_sectors()]) to
#' delete-one-**cluster**: each cluster of sectors is dropped in turn from both
#' `X` and `Y`, and the coupling is re-estimated on the remaining sectors. Under
#' cross-sectional dependence of the input-output ("MIP") kind, a naive
#' leave-one-out is optimistic because the excluded sector's information leaks in
#' through intermediate demand; dropping an entire value chain forces the
#' prediction to rely on the general gravitation rather than on near-collinear
#' neighbours. This is the cross-sectional complement of the temporal
#' time-shift / block-bootstrap nulls already in the package
#' ([rotation_null_test()], [test_permutation_robustness()]); it reuses the same
#' coupling pipeline (`.coupling_estimate`), it does not reimplement it.
#'
#' @section Cluster map (pluggable):
#' Supply `cluster_map` as the genuine input-output partition (a per-sector label
#' vector, or a named list mapping each cluster to its sector names/indices). When
#' it is `NULL`, a documented **fallback** partition is built with
#' [build_cluster_map()] (`"correlation"` or `"com"`); the fallback is a stopgap
#' proxy, not a demand-linkage partition, and a message flags its use.
#'
#' @section Statistical layer (do not over-read):
#' `bias` and `se` are the delete-a-group (block) jackknife estimates with the
#' \eqn{(g-1)/g} factor over the `g` cluster-deletion replicates. They are
#' well-calibrated for roughly balanced clusters; with strongly unequal MIP
#' clusters they are an approximate, conservative summary, and the primary
#' outputs are the per-cluster `influence`/`retention` and the `robust` verdict.
#' The verdict is purely a robustness diagnostic, not a coupling point estimate.
#'
#' @param X_matrix,Y_matrix Numeric matrices (sectors in columns; `Y` carries the
#'   sector identities used by the cluster map).
#' @param cluster_map Optional sector-to-cluster map (see "Cluster map"). `NULL`
#'   triggers the fallback.
#' @param n_clusters Integer; number of clusters for the fallback (ignored when
#'   `cluster_map` is supplied). Default `3`.
#' @param fallback Fallback rule when `cluster_map` is `NULL`: `"correlation"`
#'   (default) or `"com"`.
#' @param com Optional per-sector organic-composition vector for `fallback = "com"`.
#' @param max_comp Maximum number of PLS components (default `3`).
#' @param seed Optional integer for reproducible component selection/diagnostics.
#' @param verbose Logical; print progress and diagnostic information. Default `TRUE`.
#'
#' @return List with components:
#'   \describe{
#'     \item{`baseline`}{Full-sample coupling.}
#'     \item{`cluster_estimates`}{Named vector of leave-one-cluster coupling
#'       estimates.}
#'     \item{`cluster_members`}{Named list of the sectors in each cluster.}
#'     \item{`cluster_sizes`}{Named integer vector of cluster sizes.}
#'     \item{`bias`,`se`}{Delete-a-group jackknife bias and standard error.}
#'     \item{`influence`}{Per-cluster signed change `baseline - lco_c`.}
#'     \item{`retention`}{Per-cluster `lco_c / baseline`.}
#'     \item{`influential_clusters`}{Clusters ordered by absolute influence.}
#'     \item{`n_clusters`}{Number of clusters.}
#'     \item{`cluster_map`}{The resolved per-sector cluster labels.}
#'     \item{`map_source`}{`"user"` or `"fallback:<method>"`.}
#'     \item{`robust`}{`TRUE` if no single cluster changes the coupling by more
#'       than 50 percent.}
#'   }
#'
#' @seealso [test_jackknife_sectors()], [build_cluster_map()]
#' @export
test_leave_cluster_out <- function(X_matrix, Y_matrix, cluster_map = NULL,
                                   n_clusters = 3,
                                   fallback = c("correlation", "com"),
                                   com = NULL, max_comp = 3, seed = 123,
                                   verbose = TRUE) {
  fallback <- match.arg(fallback)
  if (verbose) {
    message("\n========================================")
    message("LEAVE-CLUSTER-OUT")
    message("========================================\n")
  }

  X_matrix <- as.matrix(X_matrix)
  Y_matrix <- as.matrix(Y_matrix)
  n_sec <- ncol(Y_matrix)
  sector_names <- colnames(Y_matrix)
  if (is.null(sector_names)) {
    sector_names <- paste0("S", seq_len(n_sec))
    colnames(Y_matrix) <- sector_names
  }

  map_source <- "user"
  if (is.null(cluster_map)) {
    cluster_labels <- build_cluster_map(Y_matrix, n_clusters = n_clusters,
                                        method = fallback, com = com)
    cluster_labels <- as.character(cluster_labels)
    names(cluster_labels) <- sector_names
    map_source <- paste0("fallback:", fallback)
    if (verbose) {
      message("No cluster_map supplied: using a documented FALLBACK partition (",
              fallback, ").")
      message("This is a stopgap proxy for the real input-output (MIP) clusters; ",
              "supply `cluster_map` once the MIP is available.\n")
    }
  } else {
    cluster_labels <- .normalize_cluster_map(cluster_map, sector_names)
  }

  clusters <- unique(cluster_labels)
  if (verbose) {
    message("Sectors: ", n_sec, " | clusters: ", length(clusters))
  }

  baseline <- .coupling_estimate(X_matrix, Y_matrix, max_comp = max_comp, seed = seed)
  if (!is.finite(baseline)) {
    if (verbose) message("FAIL: baseline coupling could not be computed")
    return(NULL)
  }

  lco <- rep(NA_real_, length(clusters)); names(lco) <- clusters
  members <- vector("list", length(clusters)); names(members) <- clusters
  sizes <- integer(length(clusters)); names(sizes) <- clusters

  if (verbose) pb <- txtProgressBar(min = 0, max = length(clusters), style = 3)
  for (ci in seq_along(clusters)) {
    c_lab <- clusters[ci]
    drop_idx <- which(cluster_labels == c_lab)
    members[[ci]] <- sector_names[drop_idx]
    sizes[ci] <- length(drop_idx)

    Yd <- Y_matrix[, -drop_idx, drop = FALSE]
    drop_idx_X <- drop_idx[drop_idx <= ncol(X_matrix)]
    Xd <- if (length(drop_idx_X) == 0L) X_matrix
          else X_matrix[, -drop_idx_X, drop = FALSE]

    lco[ci] <- .coupling_estimate(Xd, Yd, max_comp = max_comp, seed = seed)
    if (verbose) setTxtProgressBar(pb, ci)
  }
  if (verbose) close(pb)

  valid <- is.finite(lco)
  lcov <- lco[valid]
  g <- length(lcov)
  if (g < 2) {
    if (verbose) message("\nFAIL: too few valid leave-cluster-out replicates")
    return(list(baseline = baseline, cluster_estimates = lco,
                cluster_members = members, cluster_sizes = sizes,
                n_clusters = length(clusters), cluster_map = cluster_labels,
                map_source = map_source, robust = NA))
  }

  jack_mean <- mean(lcov)
  bias <- (g - 1) * (jack_mean - baseline)
  se   <- sqrt((g - 1) / g * sum((lcov - jack_mean)^2))
  influence <- baseline - lco
  retention <- lco / baseline
  rel_change <- abs(influence) / (abs(baseline) + 1e-12)

  ord <- order(abs(influence), decreasing = TRUE, na.last = NA)
  influential_clusters <- clusters[ord]
  robust <- isTRUE(max(rel_change, na.rm = TRUE) < 0.5)

  if (verbose) {
    message(sprintf("\nFull-sample coupling: %.4f", baseline))
    message(sprintf("Delete-a-group jackknife bias: %.4f | SE: %.4f", bias, se))
    message("Most influential cluster: ", influential_clusters[1])
    message(if (robust) "OK: no single value chain dominates (<50% change)"
            else        "WARNING: at least one value chain changes coupling by >50%")
  }

  list(
    baseline = baseline,
    cluster_estimates = lco,
    cluster_members = members,
    cluster_sizes = sizes,
    bias = bias,
    se = se,
    influence = influence,
    retention = retention,
    influential_clusters = influential_clusters,
    n_clusters = length(clusters),
    cluster_map = cluster_labels,
    map_source = map_source,
    robust = robust
  )
}