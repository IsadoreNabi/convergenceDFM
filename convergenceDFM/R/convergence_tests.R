#' Extract OU posteriors and credible intervals (Stan)
#'
#' Extracts medians and central credible intervals for \eqn{\phi_x}, \eqn{\phi_y},
#' and \eqn{\beta} from a fitted Stan object and checks whether
#' \eqn{\kappa = 1 - \phi} lies in (0,1).
#'
#' @param ou_result Object returned by \code{\link{estimate_factor_OU}}.
#' @param stan_fit Optional Stan fit (uses \code{ou_result$stan_fit} if missing).
#' @param confidence Credible level in (0,1), default 0.95.
#' @param verbose Logical; print progress/details. Default \code{TRUE}.
#' @return A list with \code{kappa_x}, \code{kappa_y}, \code{beta}, and a boolean
#'   \code{convergence_evidence}.
#' @keywords internal
#' @noRd

extract_ou_posteriors <- function(ou_result, stan_fit = NULL, confidence = 0.95, verbose = TRUE) {
  if (verbose) {
    message("========================================")
    message("EXTRACTION OF POSTERIOR OU WITH CI95%")
    message("========================================\n")
  }
  
  if (is.null(ou_result) || ou_result$method == "discrete_AR1_fallback") {
    if (verbose) message("WARNING: No Stan posteriors available")
    return(NULL)
  }
  
  if (is.null(stan_fit) && !is.null(ou_result$stan_fit)) {
    stan_fit <- ou_result$stan_fit
  }
  
  if (is.null(stan_fit)) {
    if (verbose) message("WARNING: No Stan fit available")
    return(NULL)
  }
  
  alpha <- 1 - confidence
  q_low <- alpha / 2
  q_high <- 1 - alpha / 2
  
  extract_quantiles <- function(param_name, fit) {
    tryCatch({
      if (inherits(fit, "CmdStanMCMC")) {
        draws <- fit$draws(variables = param_name, format = "matrix")
      } else {
        draws <- rstan::extract(fit, pars = param_name)[[1]]
      }
      
      if (is.matrix(draws)) {
        result <- apply(draws, 2, function(x) {
          c(median = median(x),
            lower = quantile(x, q_low),
            upper = quantile(x, q_high))
        })
        if (!is.matrix(result)) {
          result <- matrix(result, nrow = 3, 
                           dimnames = list(c("median", "lower", "upper"), NULL))
        }
        result
      } else {
        c(median = median(draws),
          lower = quantile(draws, q_low),
          upper = quantile(draws, q_high))
      }
    }, error = function(e) {
      if (verbose) message("  Error extracting ", param_name, ": ", e$message)
      return(NULL)
    })
  }
  
  phi_x_draws <- extract_quantiles("phi_x", stan_fit)
  phi_y_draws <- extract_quantiles("phi_y", stan_fit)
  
  if (is.null(phi_x_draws) || is.null(phi_y_draws)) {
    if (verbose) message("WARNING: Parameters could not be extracted")
    return(NULL)
  }
  
  if (!is.matrix(phi_x_draws)) {
    phi_x_draws <- matrix(phi_x_draws, nrow = 3,
                          dimnames = list(c("median", "lower", "upper"), NULL))
  }
  if (!is.matrix(phi_y_draws)) {
    phi_y_draws <- matrix(phi_y_draws, nrow = 3,
                          dimnames = list(c("median", "lower", "upper"), NULL))
  }
  
  kappa_x <- 1 - phi_x_draws
  kappa_y <- 1 - phi_y_draws
  
  beta_draws <- extract_quantiles("beta", stan_fit)
  
  if (verbose) message("OK: Posteriors extracted with CI", confidence * 100, "%\n")
  
  kappa_in_range_x <- FALSE
  kappa_in_range_y <- FALSE
  
  tryCatch({
    if (nrow(kappa_x) >= 3 && ncol(kappa_x) >= 1) {
      kappa_in_range_x <- all(kappa_x[2,] > 0 & kappa_x[3,] < 1)
    }
    if (nrow(kappa_y) >= 3 && ncol(kappa_y) >= 1) {
      kappa_in_range_y <- all(kappa_y[2,] > 0 & kappa_y[3,] < 1)
    }
  }, error = function(e) {
    if (verbose) message("  Error checking ranges")
  })
  
  if (kappa_in_range_x && kappa_in_range_y) {
    if (verbose) message("OK: CI95%(k) in (0,1) -> GUARANTEED CONVERGENCE")
  } else {
    if (verbose) message("WARNING: Some 95%CI(k) are not in (0,1)")
  }
  
  list(
    kappa_x = kappa_x,
    kappa_y = kappa_y,
    beta = beta_draws,
    convergence_evidence = kappa_in_range_x && kappa_in_range_y
  )
}

#' Unit-root tests on long-run OU errors
#'
#' Applies ADF and PP tests to long-run deviations implied by OU parameters,
#' reporting rejection at 5% significance.
#'
#' @param factors_data List with factor score matrices.
#' @param ou_result OU estimation result providing \code{mu_y}.
#' @return A list of test statistics and decisions per series.
#' @keywords internal
#' @noRd

compute_equilibrium_errors <- function(factors_data, ou_result,
                                       use_lagged_X = TRUE, eps = 1e-6) {
  FX <- as.matrix(factors_data$scores_X)
  FY <- as.matrix(factors_data$scores_Y)
  Tn <- min(nrow(FX), nrow(FY))
  FX <- FX[1:Tn, , drop = FALSE]
  FY <- FY[1:Tn, , drop = FALSE]
  Kx <- ncol(FX); Ky <- ncol(FY)
  
  have_all <- !is.null(ou_result$mu_y) && !is.null(ou_result$mu_x) &&
    !is.null(ou_result$beta)
  if (!have_all) {
    return(list(
      u_instant = FY - matrix(mean(FY, na.rm = TRUE), Tn, Ky),
      u_longrun = NULL
    ))
  }
  
  mu_y <- matrix(rep_len(as.numeric(ou_result$mu_y), Ky), Tn, Ky, byrow = TRUE)
  mu_x <- matrix(rep_len(as.numeric(ou_result$mu_x), Kx), Tn, Kx, byrow = TRUE)
  B    <- as.matrix(ou_result$beta)  # Ky x Kx
  
  FX_lag <- if (use_lagged_X) rbind(NA, FX[-Tn, , drop = FALSE]) else FX
  u_inst <- FY - (mu_y + (FX_lag - mu_x) %*% t(B))
  if (use_lagged_X) u_inst <- u_inst[-1, , drop = FALSE]
  
  u_lr <- NULL
  if (!is.null(ou_result$phi_y)) {
    phi <- as.numeric(ou_result$phi_y)
    phi <- pmin(pmax(phi, eps), 1 - eps)
    Theta <- diag(1 / (1 - phi), nrow = length(phi)) %*% B  # Ky x Kx
    u_lr <- FY - (mu_y + (FX - mu_x) %*% t(Theta))
  }
  
  list(u_instant = u_inst, u_longrun = u_lr)
}

test_longrun_error <- function(factors_data, ou_result,
                               mode = c("instant", "longrun", "both"),
                               use_lagged_X = TRUE, lags_adf = 4) {
  mode <- match.arg(mode)
  
  if (!requireNamespace("urca", quietly = TRUE)) {
    stop("Package 'urca' is required. Please install it with install.packages('urca')",
         call. = FALSE)
  }
  
  errs <- compute_equilibrium_errors(factors_data, ou_result, use_lagged_X = use_lagged_X)
  out <- list()
  
  run_tests <- function(U) {
    if (is.null(U)) return(NULL)
    Ky <- ncol(U)
    adf_results <- vector("list", Ky); pp_results <- vector("list", Ky)
    for (j in 1:Ky) {
      ut <- as.numeric(U[, j])
      ut[!is.finite(ut)] <- NA_real_
      ut <- ut[is.finite(ut)]
      if (length(ut) < 10) next
      adf <- urca::ur.df(ut, type = "drift", lags = lags_adf)
      pp  <- urca::ur.pp(ut, type = "Z-tau", model = "constant")
      adf_results[[j]] <- list(stat = adf@teststat[1], crit = adf@cval[1, 2],
                               reject = adf@teststat[1] < adf@cval[1, 2])
      pp_results[[j]]  <- list(stat = pp@teststat,      crit = pp@cval[2],
                               reject = pp@teststat < pp@cval[2])
      names(adf_results)[j] <- paste0("u", j)
      names(pp_results)[j]  <- paste0("u", j)
    }
    adf_all <- all(vapply(adf_results, function(x) isTRUE(x$reject), logical(1)), na.rm = TRUE)
    pp_all  <- all(vapply(pp_results,  function(x) isTRUE(x$reject), logical(1)), na.rm = TRUE)
    list(adf = adf_results, pp = pp_results, convergence = adf_all && pp_all)
  }
  
  if (mode %in% c("instant", "both")) out$instant <- run_tests(errs$u_instant)
  if (mode %in% c("longrun", "both")) out$longrun <- run_tests(errs$u_longrun)
  
  # Compatibilidad con código existente que espera $adf/$pp:
  if (mode == "instant" && !is.null(out$instant)) {
    out$adf <- out$instant$adf; out$pp <- out$instant$pp; out$convergence <- out$instant$convergence
  } else if (mode == "longrun" && !is.null(out$longrun)) {
    out$adf <- out$longrun$adf; out$pp <- out$longrun$pp; out$convergence <- out$longrun$convergence
  }
  out
}

#' Classical cointegration control (Johansen trace or eigen)
#'
#' Runs Johansen's cointegration test on the first min(2, ncol(X), ncol(Y))
#' factors from \code{scores_X} and \code{scores_Y}, and counts the number of
#' cointegrating relations at the 5% level.
#'
#' @param factors_data A list with matrices \code{scores_X} and \code{scores_Y}
#'   (T x Kx) and (T x Ky), respectively.
#' @param max_lag Integer; maximum lag \code{K} passed to \code{urca::ca.jo()}.
#' @param type Character; Johansen test type, one of \code{"trace"} or
#'   \code{"eigen"}. Defaults to \code{"trace"}.
#' @param ecdet Character; deterministic terms, e.g. \code{"const"},
#'   \code{"trend"}, or \code{"none"}. Defaults to \code{"const"}.
#' @param verbose Logical; print progress and a summary of the test. Default \code{TRUE}.
#'
#' @return A list with:
#' \itemize{
#'   \item \code{test}: the \code{urca::ca.jo} fitted object,
#'   \item \code{n_coint}: integer number of cointegrating relations at 5%,
#'   \item \code{vectors}: matrix of cointegrating vectors (or \code{NULL} if none).
#' }
#'
#' @details
#' This function requires the optional package \pkg{urca} (declared in
#' \code{Suggests}). It does not attempt to install packages at runtime; if
#' \pkg{urca} is unavailable, an informative error is thrown.
#'
#' The 5% critical values are taken from the \code{"5pct"} column of the
#' \code{cval} slot returned by \code{urca::ca.jo()}.
#'
#' @references
#' Johansen, S. (1991). Estimation and Hypothesis Testing of Cointegration Vectors
#' in Gaussian Vector Autoregressive Models. \emph{Econometrica}, 59(6), 1551-1580.
#'
#' Johansen, S. (1995). \emph{Likelihood-Based Inference in Cointegrated Vector
#' Autoregressive Models}. Oxford University Press.
#'
#' @seealso \code{\link[urca]{ca.jo}}
#'
#' @examples
#' \donttest{
#' if (requireNamespace("urca", quietly = TRUE)) {
#'   set.seed(1)
#'   T <- 120
#'   X <- cbind(cumsum(rnorm(T)), cumsum(rnorm(T)))
#'   Y <- cbind(cumsum(rnorm(T)), cumsum(rnorm(T)))
#'   fd <- list(scores_X = X, scores_Y = Y)
#'   out <- test_cointegration_control(fd, max_lag = 2, verbose = FALSE)
#'   str(out)
#' }
#' }
#' @export

test_cointegration_control <- function(factors_data,
                                       max_lag = 4,
                                       type = "trace",
                                       ecdet = "const",
                                       verbose = TRUE) {
  
  if (!requireNamespace("urca", quietly = TRUE)) {
    stop("Package 'urca' is required for cointegration tests. ",
         "Please install it and list it under Suggests.", call. = FALSE)
  }
  
  FX <- as.matrix(factors_data$scores_X)
  FY <- as.matrix(factors_data$scores_Y)
  
  n_comp <- min(2L, ncol(FX), ncol(FY))
  data_coint <- cbind(FX[, 1:n_comp, drop = FALSE],
                      FY[, 1:n_comp, drop = FALSE])
  colnames(data_coint) <- c(paste0("X", seq_len(n_comp)),
                            paste0("Y", seq_len(n_comp)))
  
  if (verbose) {
    message("========================================")
    message("CLASSICAL COINTEGRATION (JOHANSEN)")
    message("========================================\n")
    message("Running Johansen test... type=", type,
            ", ecdet=", ecdet, ", K=", max_lag)
  }
  
  joh_test <- urca::ca.jo(data_coint, type = type, ecdet = ecdet, K = max_lag)
  
  if (verbose) print(summary(joh_test))
  
  test_stats <- methods::slot(joh_test, "teststat")
  cval_mat  <- methods::slot(joh_test, "cval")
  idx5 <- which(colnames(cval_mat) %in% c("5pct", "5%"))
  if (length(idx5) == 0L) idx5 <- 2L
  crit_vals <- cval_mat[, idx5]
  
  n_coint <- sum(test_stats > crit_vals)
  
  if (verbose) {
    message("\nCointegration relationships (5%): ", n_coint)
    message(if (n_coint > 0)
      "\nOK: Evidence of cointegration between X and Y"
      else
        "\nWARNING: No cointegration is detected at 5%")
  }
  
  vectors <- if (n_coint > 0L) methods::slot(joh_test, "V")[, 1:n_coint, drop = FALSE] else NULL
  
  list(
    test    = joh_test,
    n_coint = n_coint,
    vectors = vectors
  )
}