#' Generate a Haar-uniform random orthogonal matrix
#'
#' Generates a random orthogonal matrix of size k x k from the QR decomposition
#' of a Gaussian matrix, with the sign correction of Mezzadri (2007) so that the
#' result is uniform on the orthogonal group (the raw \code{qr.Q} output is not,
#' because LAPACK fixes the signs of the R diagonal).
#'
#' @param k Integer. Matrix dimension.
#'
#' @return k x k orthogonal matrix.
#'
#' @keywords internal
#' @noRd

random_orthogonal <- function(k) {
  Z <- matrix(rnorm(k * k), k, k)
  qrz <- qr(Z)
  Q <- qr.Q(qrz)
  R <- qr.R(qrz)
  d <- diag(R)
  ph <- ifelse(d == 0, 1, sign(d))
  Q <- Q * rep(ph, each = k)        # Haar correction
  if (det(Q) < 0) Q[, 1] <- -Q[, 1] # restrict to SO(k)
  Q
}

#' Circularly shift the rows of a matrix
#'
#' Wraps the rows of \code{M} by \code{tau} positions. Used to build a null that
#' destroys cross-series temporal alignment while preserving each series' own
#' (circular) autocorrelation structure -- the appropriate null for testing
#' coupling with rotation-invariant statistics.
#'
#' @param M Numeric matrix (T x k).
#' @param tau Integer shift in \code{1..T-1}.
#' @return Row-shifted matrix.
#' @keywords internal
#' @noRd

circ_shift_rows <- function(M, tau) {
  Tn <- nrow(M)
  if (Tn < 2L) return(M)
  tau <- ((tau %% Tn) + Tn) %% Tn
  if (tau == 0L) tau <- 1L
  idx <- ((seq_len(Tn) - 1L + tau) %% Tn) + 1L
  M[idx, , drop = FALSE]
}

#' One-sided HAC test that the mean of a loss-differential series is positive
#'
#' @param d Numeric vector of per-period loss differentials.
#' @return List with \code{stat} (asymptotically N(0,1)), one-sided \code{p}, and
#'   the effective length \code{n}.
#' @keywords internal
#' @noRd
hac_mean_test <- function(d) {
  d <- d[is.finite(d)]
  n <- length(d)
  if (n < 8) return(list(stat = NA_real_, p = NA_real_, n = n))
  dbar <- mean(d)
  dc <- d - dbar
  L <- max(1L, floor(4 * (n / 100)^(2 / 9)))   # Bartlett bandwidth
  lrv <- sum(dc^2) / n
  for (l in seq_len(L)) {
    w <- 1 - l / (L + 1)
    gl <- sum(dc[(l + 1):n] * dc[1:(n - l)]) / n
    lrv <- lrv + 2 * w * gl
  }
  if (!is.finite(lrv) || lrv <= 0) return(list(stat = NA_real_, p = NA_real_, n = n))
  stat <- dbar / sqrt(lrv / n)
  list(stat = stat, p = stats::pnorm(stat, lower.tail = FALSE), n = n)
}

#' Diebold-Mariano / Clark-West test for nested out-of-sample forecasts
#'
#' One-sided test that the larger model (with X) improves squared-error loss over
#' the nested baseline. For nested models the Clark-West adjustment is used;
#' otherwise the standard Diebold-Mariano statistic, both via a HAC long-run
#' variance.
#'
#' @param e_base,e_full Numeric vectors of one-step forecast errors (baseline and
#'   full model), same length.
#' @param yhat_base,yhat_full Optional forecast values, required for the
#'   Clark-West adjustment when \code{nested = TRUE}.
#' @param nested Logical; use the Clark-West adjustment.
#' @return List with \code{stat} (asymptotically N(0,1)) and one-sided \code{p}.
#' @keywords internal
#' @noRd

dm_cw_test <- function(e_base, e_full, yhat_base = NULL, yhat_full = NULL,
                       nested = TRUE) {
  ok <- is.finite(e_base) & is.finite(e_full)
  e_base <- e_base[ok]; e_full <- e_full[ok]
  if (nested && !is.null(yhat_base) && !is.null(yhat_full)) {
    yhat_base <- yhat_base[ok]; yhat_full <- yhat_full[ok]
    d <- e_base^2 - e_full^2 + (yhat_base - yhat_full)^2   # Clark-West
  } else {
    d <- e_base^2 - e_full^2                               # Diebold-Mariano
  }
  hac_mean_test(d)
}

#' Orthonormalize matrix via QR
#'
#' Orthonormalizes the columns of a matrix using QR decomposition.
#'
#' @param M Numeric matrix to orthonormalize.
#'
#' @return Matrix with orthonormal columns.
#'
#' @keywords internal
#' @noRd

orthonormalize <- function(M) {
  Q <- tryCatch(qr.Q(qr(as.matrix(M))), error=function(e) NULL)
  if (is.null(Q)) Q <- scale(as.matrix(M), center=TRUE, scale=TRUE)
  Q
}

#' Safe SVD with regularization
#'
#' Computes singular value decomposition with automatic regularization if
#' standard SVD fails.
#'
#' @param A Numeric matrix.
#'
#' @return SVD object.
#'
#' @keywords internal
#' @noRd

safe_svd <- function(A) tryCatch(svd(A), error=function(e) svd(A + 1e-10*diag(nrow(A))))
ridge_solve <- function(X, Y, lambda=1e-8) {
  XtX <- crossprod(X); p <- ncol(X)
  tryCatch(solve(XtX, crossprod(X, Y)),
           error=function(e) solve(XtX + diag(lambda, p), crossprod(X, Y)))
}

#' Procrustes distance metric
#'
#' Computes Procrustes distance between two factor spaces after optimal
#' orthogonal alignment.
#'
#' @param X Numeric matrix (T x K1).
#' @param Y Numeric matrix (T x K2).
#'
#' @return Numeric scalar. Procrustes distance.
#'
#' @keywords internal
#' @noRd

procrustes_metric <- function(X, Y) {
  Xc <- scale(as.matrix(X), center=TRUE, scale=FALSE)
  Yc <- scale(as.matrix(Y), center=TRUE, scale=FALSE)
  S  <- crossprod(Xc, Yc)
  sv <- safe_svd(S)$d
  R  <- sum(sv) / (sqrt(sum(Xc^2)) * sqrt(sum(Yc^2)) + .Machine$double.eps)
  list(R = as.numeric(R), trace = sum(sv))
}

#' Canonical correlation metric
#'
#' Computes maximum canonical correlation between two factor spaces as a
#' measure of coupling strength.
#'
#' @param X Numeric matrix (T x K1).
#' @param Y Numeric matrix (T x K2).
#'
#' @return Numeric scalar between 0 and 1. Maximum canonical correlation.
#'
#' @keywords internal
#' @noRd

cca_metric <- function(X, Y) {
  Xc <- scale(as.matrix(X), center=TRUE, scale=TRUE)
  Yc <- scale(as.matrix(Y), center=TRUE, scale=TRUE)
  kmax <- min(nrow(Xc)-1, ncol(Xc), ncol(Yc))
  if (kmax < 1) return(list(cor=numeric(0), sum_r2=0, max_r=0))
  pca_reduce <- function(M, k) {
    sv <- safe_svd(scale(M, center=TRUE, scale=FALSE))
    U <- sv$u[, seq_len(k), drop=FALSE]; D <- diag(sv$d[seq_len(k)], nrow=k)
    U %*% D
  }
  if (nrow(Xc) <= (ncol(Xc)+ncol(Yc)+1)) {
    Xc <- pca_reduce(Xc, min(kmax, ncol(Xc)))
    Yc <- pca_reduce(Yc, min(kmax, ncol(Yc)))
  }
  cc <- cancor(Xc, Yc)
  cors <- cc$cor
  list(cor=cors, sum_r2=sum(cors^2), max_r=ifelse(length(cors), max(cors), 0))
}

#' Principal angles metric
#'
#' Computes principal (canonical) angles between two subspaces, measuring
#' their geometric alignment.
#'
#' @param X Numeric matrix defining first subspace.
#' @param Y Numeric matrix defining second subspace.
#'
#' @return Numeric vector of principal angles (in radians).
#'
#' @keywords internal
#' @noRd

principal_angles_metric <- function(X, Y) {
  Qx <- orthonormalize(X); Qy <- orthonormalize(Y)
  k <- min(ncol(Qx), ncol(Qy)); if (k==0) return(list(mean_cos=0, sum_cos2=0, max_cos=0))
  sv <- safe_svd(crossprod(Qx, Qy))$d; sv <- pmin(pmax(sv,-1),1)
  list(mean_cos=mean(sv), sum_cos2=sum(sv^2), max_cos=max(sv), cosines=sv)
}

#' Dynamic beta norm for short-run coupling
#'
#' Computes the Frobenius norm of lagged regression coefficients from Y on X
#' as a measure of short-run dynamic coupling.
#'
#' @param Y_scores Numeric matrix of Y factor scores.
#' @param X_scores Numeric matrix of X factor scores.
#' @param lag Integer. Lag order. Default is 1.
#'
#' @return Numeric scalar. Norm of dynamic coefficients.
#'
#' @keywords internal
#' @noRd

dyn_beta_norm <- function(Y_scores, X_scores, lag=1) {
  stopifnot(nrow(Y_scores)==nrow(X_scores))
  Tn <- nrow(Y_scores); if (Tn <= lag) return(NA_real_)
  Yt <- Y_scores[(lag+1):Tn,, drop=FALSE]
  Yl <- Y_scores[1:(Tn-lag),  , drop=FALSE]
  Xl <- X_scores[1:(Tn-lag),  , drop=FALSE]
  D  <- cbind(Yl, Xl)
  B  <- ridge_solve(D, Yt)
  kY <- ncol(Yl); kX <- ncol(Xl)
  Bxy <- B[(kY+1):nrow(B), , drop=FALSE]
  sqrt(sum(Bxy^2)) / sqrt(kX * ncol(Yt))
}

#' Extract factor scores from results object
#'
#' Utility function to extract X and Y factor scores from a complete analysis
#' results object with validation.
#'
#' @param results_robust List. Output from \code{run_complete_factor_analysis_robust}.
#'
#' @return List with \code{scores_X} and \code{scores_Y} matrices.
#'
#' @keywords internal
#' @noRd

extract_scores_from_results <- function(results_robust) {
  if (!is.null(results_robust$factors) &&
      !is.null(results_robust$factors$scores_X) &&
      !is.null(results_robust$factors$scores_Y)) {
    return(list(
      scores_X = as.matrix(results_robust$factors$scores_X),
      scores_Y = as.matrix(results_robust$factors$scores_Y)
    ))
  }
  if (!is.null(results_robust$selection) &&
      !is.null(results_robust$selection$pls_model)) {
    pls_model <- results_robust$selection$pls_model
    k <- results_robust$selection$optimal_ncomp
    Sx <- pls_model$scores[,  1:k, drop=FALSE]
    Sy <- pls_model$Yscores[, 1:k, drop=FALSE]
    colnames(Sx) <- paste0("X", 1:ncol(Sx))
    colnames(Sy) <- paste0("Y", 1:ncol(Sy))
    return(list(scores_X = as.matrix(Sx), scores_Y = as.matrix(Sy)))
  }
  stop("No scores_X/scores_Y found in results_robust or pls_model.")
}

#' Null hypothesis test for X->Y factor coupling
#'
#' Tests whether the observed association between the lagged X factor space and
#' the Y factor space is stronger than expected when the X->Y temporal alignment
#' is broken.
#'
#' @section Why not a rotation null:
#' The coupling statistics used here (Procrustes nuclear-norm ratio, canonical
#' correlations, principal angles, dynamic-beta Frobenius norm) are
#' \emph{invariant} to orthogonal rotation of either factor space: rotating Y by
#' an orthogonal matrix leaves all of them unchanged, so a rotation-based null is
#' degenerate (its "distribution" is a point mass at the observed value and every
#' p-value is ~1). The default null (\code{null_method = "circular_shift"})
#' instead applies a random circular time shift to Y, which preserves each
#' series' own autocorrelation while destroying the cross-series alignment -- the
#' statistics then genuinely vary under the null. \code{null_method = "rotation"}
#' is retained only as an invariance diagnostic and will (correctly) yield
#' non-significant results.
#'
#' @param scores_X Factor scores from the first dataset (T x Kx).
#' @param scores_Y Factor scores from the second dataset (T x Ky).
#' @param lag Number of lags for the X->Y relationship (default 1).
#' @param B Number of null replicates (default 1000).
#' @param seed Random seed for reproducibility (default 123). The RNG state is
#'   restored on exit.
#' @param compute Statistics to compute: any of 'procrustes', 'cca', 'principal',
#'   'dynbeta'.
#' @param null_method Null-generating mechanism: "circular_shift" (default) or
#'   "rotation" (invariance diagnostic only).
#' @param progress Logical, show progress bar (default TRUE).
#'
#' @return List with components:
#'   \describe{
#'     \item{\code{observed}}{Observed statistics.}
#'     \item{\code{null_stats}}{List of length-\code{B} null vectors per statistic.}
#'     \item{\code{p_values}}{Monte Carlo one-sided p-values \code{(1 + #{null >=
#'       obs}) / (B + 1)}.}
#'     \item{\code{p_values_fdr}}{Benjamini-Hochberg FDR-adjusted p-values.}
#'     \item{\code{params}}{Test settings, including the seed and null method.}
#'   }
#' @export

rotation_null_test <- function(scores_X, scores_Y, lag = 1, B = 1000, seed = 123,
                               compute = c("procrustes", "cca", "principal", "dynbeta"),
                               null_method = c("circular_shift", "rotation"),
                               progress = TRUE) {
  null_method <- match.arg(null_method)

  if (!is.null(seed)) {
    if (exists(".Random.seed", envir = .GlobalEnv)) {
      old_seed <- get(".Random.seed", envir = .GlobalEnv)
      on.exit(assign(".Random.seed", old_seed, envir = .GlobalEnv), add = TRUE)
    }
    set.seed(seed)
  }

  X <- as.matrix(scores_X); Y <- as.matrix(scores_Y)
  stopifnot(nrow(X) == nrow(Y))
  Tn <- nrow(X); if (Tn <= lag + 1) stop("Series too short for the requested lag.")

  # All coupling statistics for a given (X, Y) pair, returned as a named vector.
  compute_stats <- function(Xfull, Yfull) {
    Xl <- Xfull[1:(Tn - lag), , drop = FALSE]
    Yt <- Yfull[(lag + 1):Tn, , drop = FALSE]
    out <- c(procrustes_R = NA_real_, cca_sumr2 = NA_real_, cca_maxr = NA_real_,
             princ_sum2 = NA_real_, princ_max = NA_real_, dynbeta = NA_real_)
    if ("procrustes" %in% compute) out["procrustes_R"] <- procrustes_metric(Xl, Yt)$R
    if ("cca" %in% compute) {
      cc <- cca_metric(Xl, Yt); out["cca_sumr2"] <- cc$sum_r2; out["cca_maxr"] <- cc$max_r
    }
    if ("principal" %in% compute) {
      pa <- principal_angles_metric(Xl, Yt)
      out["princ_sum2"] <- pa$sum_cos2; out["princ_max"] <- pa$max_cos
    }
    if ("dynbeta" %in% compute) out["dynbeta"] <- dyn_beta_norm(Yfull, Xfull, lag)
    out
  }

  obs_vec <- compute_stats(X, Y)
  stat_names <- names(obs_vec)
  null_mat <- matrix(NA_real_, B, length(stat_names),
                     dimnames = list(NULL, stat_names))

  if (progress) pb <- txtProgressBar(min = 0, max = B, style = 3)
  for (b in seq_len(B)) {
    if (null_method == "circular_shift") {
      tau <- sample.int(Tn - 1L, 1L)
      Yb <- circ_shift_rows(Y, tau)
    } else {                       # rotation: invariance diagnostic
      Yb <- Y %*% random_orthogonal(ncol(Y))
    }
    null_mat[b, ] <- compute_stats(X, Yb)
    if (progress) setTxtProgressBar(pb, b)
  }
  if (progress) close(pb)

  # Monte Carlo one-sided p-values with the (1 + .)/(B + 1) correction.
  pvals <- vapply(stat_names, function(nm) {
    if (!is.finite(obs_vec[[nm]])) return(NA_real_)
    nd <- null_mat[, nm]; nd <- nd[is.finite(nd)]
    if (!length(nd)) return(NA_real_)
    (1 + sum(nd >= obs_vec[[nm]])) / (length(nd) + 1)
  }, numeric(1))

  keep <- !is.na(pvals)
  pvals_fdr <- pvals
  if (any(keep)) pvals_fdr[keep] <- stats::p.adjust(pvals[keep], method = "BH")

  observed <- list(
    procrustes = list(R = obs_vec[["procrustes_R"]]),
    cca        = list(sum_r2 = obs_vec[["cca_sumr2"]], max_r = obs_vec[["cca_maxr"]]),
    principal  = list(sum_cos2 = obs_vec[["princ_sum2"]], max_cos = obs_vec[["princ_max"]]),
    dynbeta    = obs_vec[["dynbeta"]]
  )
  null_stats <- as.list(as.data.frame(null_mat))

  list(observed = observed, observed_vec = obs_vec,
       null_stats = null_stats, p_values = as.list(pvals),
       p_values_fdr = as.list(pvals_fdr),
       params = list(B = B, lag = lag, null_method = null_method, seed = seed))
}

#' Run the coupling null test on complete analysis results
#'
#' Wrapper that extracts factors from a complete analysis object and runs
#' \code{\link{rotation_null_test}} with the corrected time-shift null.
#'
#' @return List. Output from \code{\link{rotation_null_test}}.
#'
#' @param results_robust Output from \code{run_complete_factor_analysis_robust()}.
#' @param lag Number of lags for the X->Y relationship (default 1).
#' @param B Number of null replicates (default 1000).
#' @param seed Random seed for reproducibility (default 42).
#' @param null_method Null mechanism passed to \code{rotation_null_test()};
#'   "circular_shift" (default) or "rotation" (diagnostic).
#' @param compute Statistics to compute: 'procrustes', 'cca', 'principal',
#'   'dynbeta'.
#' @export

run_rotation_null_on_results <- function(results_robust, lag = 1, B = 1000, seed = 42,
                                         null_method = c("circular_shift", "rotation"),
                                         compute = c("procrustes", "cca", "principal", "dynbeta")) {
  null_method <- match.arg(null_method)
  sc <- extract_scores_from_results(results_robust)
  Tn <- min(nrow(sc$scores_X), nrow(sc$scores_Y))
  sc$scores_X <- sc$scores_X[seq_len(Tn), , drop = FALSE]
  sc$scores_Y <- sc$scores_Y[seq_len(Tn), , drop = FALSE]
  out <- rotation_null_test(sc$scores_X, sc$scores_Y, lag = lag, B = B, seed = seed,
                            compute = compute, null_method = null_method)
  results_robust$convergence_tests$rotation_null <- out
  out
}

#' Rescue short-run channel test
#'
#' Evaluates the short-run relationship from X to Y factors using three
#' complementary, valid pieces of evidence: (1) a coupling null obtained by
#' circularly time-shifting the X residual (breaking its alignment with Y while
#' preserving its own autocorrelation -- the rotation-based null used previously
#' was invariant and therefore vacuous); (2) an expanding-window out-of-sample
#' comparison of forecasts with and without lagged X, with a Clark-West test of
#' the improvement; and (3) a Granger causality test.
#'
#' @return List with components:
#'   \describe{
#'     \item{\code{p_values}, \code{p_values_fdr}}{Coupling p-values (time-shift
#'       null) and their Benjamini-Hochberg adjustment.}
#'     \item{\code{OOS}}{Out-of-sample RMSE/R2 comparison plus the Clark-West
#'       test (\code{cw_stat}, \code{cw_p}).}
#'     \item{\code{granger_p}}{Granger causality p-value.}
#'     \item{\code{observed}, \code{sim_quantiles}, \code{null_stats}}{Coupling
#'       statistics and their null distributions.}
#'   }
#'
#' @param results_robust Output from run_complete_factor_analysis_robust()
#' @param lag Number of lags for the model (default: 1)
#' @param B Number of null replicates (default: 1000)
#' @param seed Random seed for reproducibility (default: NULL)
#' @param ridge Ridge regularization parameter (default: 0.001)
#' @param oos_start Proportion of data to use for training (default: 0.6)
#' @param verbose Logical; print progress and diagnostic information. Default \code{TRUE}.
#' @export

rescue_short_run_channel <- function(results_robust, lag = 1, B = 1000, seed = NULL,
                                     ridge = 1e-3, oos_start = 0.6, verbose = TRUE) {
  if (!is.null(seed)) {
    if (exists(".Random.seed", envir = .GlobalEnv)) {
      old_seed <- get(".Random.seed", envir = .GlobalEnv)
      on.exit(assign(".Random.seed", old_seed, envir = .GlobalEnv), add = TRUE)
    }
    set.seed(seed)
  }

  fx <- tryCatch(as.matrix(results_robust$factors$scores_X), error = function(e) NULL)
  fy <- tryCatch(as.matrix(results_robust$factors$scores_Y), error = function(e) NULL)
  if (is.null(fx) || is.null(fy)) stop("scores_X / scores_Y not found in results_robust$factors.")
  
  Tn <- min(nrow(fx), nrow(fy))
  fx <- fx[1:Tn, , drop = FALSE]
  fy <- fy[1:Tn, , drop = FALSE]
  kx <- ncol(fx); ky <- ncol(fy)
  if (kx < 1 || ky < 1) stop("Invalid factor dimensions.")
  
  L <- as.integer(lag); if (L < 1) L <- 1L
  
  zscore <- function(M) { sc <- scale(M); sc[is.na(sc)] <- 0; sc }
  fx <- zscore(fx); fy <- zscore(fy)
  
  idx  <- (L + 1):Tn
  Yt   <- fy[idx, , drop = FALSE]
  Ylag <- fy[idx - L, , drop = FALSE]
  Xlag <- fx[idx - L, , drop = FALSE]
  
  ridge_fit <- function(Z, Y, lambda = ridge) {
    Z <- as.matrix(Z); Y <- as.matrix(Y)
    Zc <- cbind(1, Z)
    p <- ncol(Zc)
    Pen <- diag(p); Pen[1, 1] <- 0
    Bt <- tryCatch(
      solve(crossprod(Zc) + lambda * Pen, crossprod(Zc, Y)),
      error = function(e) solve(crossprod(Zc) + (lambda + 1e-6) * Pen, crossprod(Zc, Y))
    )
    list(B = Bt, fitted = Zc %*% Bt, resid = Y - (Zc %*% Bt))
  }
  
  resY  <- ridge_fit(Ylag, Yt, lambda = ridge);  Y_res <- resY$resid
  resX  <- ridge_fit(Ylag, Xlag, lambda = ridge); X_res <- resX$resid
  
  orth <- function(M) {
    M <- as.matrix(scale(M, center = TRUE, scale = FALSE))
    qrobj <- qr(M)
    q <- try(qr.Q(qrobj), silent = TRUE)
    if (inherits(q, "try-error")) {
      sv <- svd(M); q <- sv$u
    }
    q[is.na(q)] <- 0
    q
  }
  Qx <- orth(X_res); Qy <- orth(Y_res)
  rdim <- min(ncol(Qx), ncol(Qy)); if (rdim < 1) stop("Rango insuficiente tras residuar.")
  
  S <- svd(t(Qx[, 1:rdim, drop = FALSE]) %*% Qy[, 1:rdim, drop = FALSE])$d
  S <- pmin(pmax(S, 0), 1)
  princ_sum2_obs <- sum(S^2)
  princ_max_obs  <- max(S)
  
  cca <- tryCatch({
    cc <- cancor(X_res, Y_res)
    list(cor = cc$cor, sum_r2 = sum(cc$cor^2), max_r = max(cc$cor))
  }, error = function(e) list(cor = S, sum_r2 = princ_sum2_obs, max_r = princ_max_obs))
  
  cross  <- crossprod(scale(Y_res, center = TRUE, scale = FALSE),
                      scale(X_res, center = TRUE, scale = FALSE))
  svals  <- svd(cross)$d
  procr_R_obs <- sum(svals) / sqrt(sum(X_res^2) * sum(Y_res^2) + 1e-12)
  
  bigZ <- cbind(Ylag, Xlag)
  fit_dyn <- ridge_fit(bigZ, Yt, lambda = ridge)
  Bhat <- fit_dyn$B
  BX   <- Bhat[1 + ky + seq_len(kx), , drop = FALSE]
  dynbeta_obs <- sqrt(sum(BX^2))
  
  nboot <- as.integer(B)
  if (verbose) {
    message(sprintf("Generating null by circular time-shift of the X residual (B=%d)...", nboot))
  }

  procr_n     <- numeric(nboot)
  cca_sum2_n  <- numeric(nboot)
  cca_max_n   <- numeric(nboot)
  princ_sum2_n<- numeric(nboot)
  princ_max_n <- numeric(nboot)
  dynbeta_n   <- numeric(nboot)

  nX <- nrow(X_res)
  for (b in seq_len(nboot)) {
    # Break the X->Y alignment while preserving X's own autocorrelation.
    # The SAME shift is applied to the residual (subspace metrics) and to the
    # raw lagged X (dynamic-beta), so all statistics use a coherent null draw.
    tau_b <- sample.int(max(2L, nX) - 1L, 1L)
    X_rot   <- circ_shift_rows(X_res, tau_b)
    Xlag_b  <- circ_shift_rows(Xlag,  tau_b)

    Qxr <- orth(X_rot)
    Sr  <- svd(t(Qxr[, 1:rdim, drop = FALSE]) %*% Qy[, 1:rdim, drop = FALSE])$d
    Sr  <- pmin(pmax(Sr, 0), 1)
    
    princ_sum2_n[b] <- sum(Sr^2)
    princ_max_n[b]  <- max(Sr)
    
    cc_r <- tryCatch(cancor(X_rot, Y_res), error = function(e) NULL)
    if (is.null(cc_r)) {
      cca_sum2_n[b] <- princ_sum2_n[b]
      cca_max_n[b]  <- princ_max_n[b]
    } else {
      cca_sum2_n[b] <- sum(cc_r$cor^2)
      cca_max_n[b]  <- max(cc_r$cor)
    }
    
    cross_r <- crossprod(scale(Y_res, center = TRUE, scale = FALSE),
                         scale(X_rot, center = TRUE, scale = FALSE))
    svals_r <- svd(cross_r)$d
    procr_n[b] <- sum(svals_r) / sqrt(sum(X_rot^2) * sum(Y_res^2) + 1e-12)
    
    fit_r <- ridge_fit(cbind(Ylag, Xlag_b), Yt, lambda = ridge)
    Br    <- fit_r$B
    BXr   <- Br[1 + ky + seq_len(kx), , drop = FALSE]
    dynbeta_n[b] <- sqrt(sum(BXr^2))
  }

  mc_p <- function(null, obs) {
    null <- null[is.finite(null)]
    if (!length(null) || !is.finite(obs)) return(NA_real_)
    (1 + sum(null >= obs)) / (length(null) + 1)
  }
  pvals <- c(
    procrustes_R = mc_p(procr_n,     procr_R_obs),
    cca_sumr2    = mc_p(cca_sum2_n,  cca$sum_r2),
    cca_maxr     = mc_p(cca_max_n,   cca$max_r),
    princ_sum2   = mc_p(princ_sum2_n, princ_sum2_obs),
    princ_max    = mc_p(princ_max_n,  princ_max_obs),
    dynbeta      = mc_p(dynbeta_n,    dynbeta_obs)
  )
  pvals_fdr <- pvals
  keep <- is.finite(pvals)
  if (any(keep)) pvals_fdr[keep] <- stats::p.adjust(pvals[keep], method = "BH")
  
  observed <- list(
    procrustes = list(R = as.numeric(procr_R_obs)),
    principal  = list(sum_cos2 = as.numeric(princ_sum2_obs),
                      max_cos  = as.numeric(princ_max_obs),
                      cosines  = as.numeric(S)),
    cca        = list(cor   = as.numeric(cca$cor),
                      sum_r2= as.numeric(cca$sum_r2),
                      max_r = as.numeric(cca$max_r)),
    dynbeta    = as.numeric(dynbeta_obs)
  )
  
  sim_quant <- list(
    procrustes_R = quantile(procr_n,      c(.25, .5, .75, .95), na.rm = TRUE),
    cca_sumr2    = quantile(cca_sum2_n,   c(.25, .5, .75, .95), na.rm = TRUE),
    cca_maxr     = quantile(cca_max_n,    c(.25, .5, .75, .95), na.rm = TRUE),
    princ_sum2   = quantile(princ_sum2_n, c(.25, .5, .75, .95), na.rm = TRUE),
    princ_max    = quantile(princ_max_n,  c(.25, .5, .75, .95), na.rm = TRUE),
    dynbeta      = quantile(dynbeta_n,    c(.25, .5, .75, .95), na.rm = TRUE)
  )
  
  if (verbose) {
    message("OOS (expanding window) forecast with/without X_{t-1} ...")
  }
  Tobs <- nrow(Yt)
  start_idx <- max(5, floor(Tobs * oos_start))
  if (start_idx >= Tobs) start_idx <- floor(Tobs * 0.7)
  
  err_y <- c(); err_yx <- c(); cw_adj <- c()
  y_oos <- matrix(NA, Tobs - start_idx, ky)
  row <- 0
  for (i in (start_idx + 1):Tobs) {
    tr <- 1:(i - 1); te <- i

    fit_b <- ridge_fit(Ylag[tr, , drop = FALSE], Yt[tr, , drop = FALSE], lambda = ridge)
    pred_b <- cbind(1, Ylag[te, , drop = FALSE]) %*% fit_b$B

    fit_a <- ridge_fit(cbind(Ylag[tr, , drop = FALSE], Xlag[tr, , drop = FALSE]),
                       Yt[tr, , drop = FALSE], lambda = ridge)
    pred_a <- cbind(1, Ylag[te, , drop = FALSE], Xlag[te, , drop = FALSE]) %*% fit_a$B

    row <- row + 1
    y_oos[row, ] <- Yt[te, , drop = FALSE]
    err_y  <- c(err_y,  sum((Yt[te, ] - pred_b)^2))
    err_yx <- c(err_yx, sum((Yt[te, ] - pred_a)^2))
    cw_adj <- c(cw_adj, sum((pred_b - pred_a)^2))   # Clark-West adjustment
  }
  rmse_y  <- sqrt(mean(err_y))
  rmse_yx <- sqrt(mean(err_yx))
  sst <- sum((y_oos - matrix(colMeans(y_oos), nrow(y_oos), ky, byrow = TRUE))^2)
  r2_oos_y  <- 1 - sum(err_y)  / sst
  r2_oos_yx <- 1 - sum(err_yx) / sst

  # Clark-West test that adding lagged X significantly improves the (nested) OOS
  # forecast -- replaces the bare RMSE delta with an inferential statement.
  cw <- hac_mean_test(err_y - err_yx + cw_adj)

  OOS <- list(
    rmse_y = as.numeric(rmse_y),
    rmse_yx = as.numeric(rmse_yx),
    delta_rmse_pct = 100 * (rmse_y - rmse_yx) / rmse_y,
    r2_oos_y = as.numeric(r2_oos_y),
    r2_oos_yx = as.numeric(r2_oos_yx),
    delta_r2_oos = as.numeric(r2_oos_yx - r2_oos_y),
    cw_stat = cw$stat,
    cw_p = cw$p
  )
  
  gr_p <- NA
  if (requireNamespace("vars", quietly = TRUE)) {
    df <- as.data.frame(cbind(X = fx, Y = fy))
    names(df) <- c(paste0("X", seq_len(kx)), paste0("Y", seq_len(ky)))
    varfit <- tryCatch(vars::VAR(df, p = L, type = "const"), error = function(e) NULL)
    if (!is.null(varfit)) {
      cg <- tryCatch(vars::causality(varfit, cause = paste0("X", seq_len(kx))), error = function(e) NULL)
      if (!is.null(cg) && !is.null(cg$Granger)) gr_p <- cg$Granger$p.value
    }
  }
  
  list(
    p_values = pvals,
    p_values_fdr = pvals_fdr,
    observed = observed,
    sim_quantiles = sim_quant,
    null_stats = list(
      procrustes_R = procr_n,
      cca_sumr2 = cca_sum2_n,
      cca_maxr = cca_max_n,
      princ_sum2 = princ_sum2_n,
      princ_max = princ_max_n,
      dynbeta = dynbeta_n
    ),
    OOS = OOS,
    granger_p = matrix(gr_p, nrow = 1, dimnames = list(NULL, "p_value")),
    null_method = "circular_shift"
  )
}

#' Incremental R-squared from X in OU model
#'
#' Computes the incremental explanatory power (delta R-squared) contributed by
#' X factors in predicting Y factors, both in-sample and out-of-sample.
#'
#'
#' @return List with components:
#'   \describe{
#'     \item{\code{in_sample}}{Data frame with R2_full, R2_baseline, and deltaR2.}
#'     \item{\code{per_equation}}{Vector of deltaR2 for each Y factor equation.}
#'     \item{\code{OOS}}{List with out-of-sample RMSE and deltaR2 (if \code{oos = TRUE}).}
#'   }
#'
#'
#' @param results_robust Output from run_complete_factor_analysis_robust()
#' @param lag Number of lags for the model (default: 1)
#' @param oos Logical, perform out-of-sample validation (default: TRUE)
#' @param seed Random seed for reproducibility (default: 123)
#' @param ridge Ridge regularization parameter (default: 1e-08)
#' @param verbose Logical, print progress messages (default: TRUE)
#' @export

deltaR2_ou <- function(results_robust, lag = 1, oos = TRUE, seed = 123,
                       ridge = 1e-8, verbose = TRUE) {
  if (is.null(results_robust$factors) ||
      is.null(results_robust$factors$scores_X) ||
      is.null(results_robust$factors$scores_Y)) {
    stop("No encuentro scores_X/scores_Y en results_robust$factors.")
  }
  FX <- as.matrix(results_robust$factors$scores_X)
  FY <- as.matrix(results_robust$factors$scores_Y)
  
  Tn <- min(nrow(FX), nrow(FY))
  if (Tn <= lag + 2) stop("Serie demasiado corta para el lag solicitado.")
  FX <- FX[1:Tn, , drop = FALSE]
  FY <- FY[1:Tn, , drop = FALSE]
  
  nzv <- function(M) {
    v <- apply(M, 2, stats::var, na.rm = TRUE)
    keep <- which(is.finite(v) & v > 1e-10)
    if (!length(keep)) return(matrix(numeric(0), nrow(M), 0))
    M[, keep, drop = FALSE]
  }
  FX <- nzv(FX); FY <- nzv(FY)
  Kx <- ncol(FX); Ky <- ncol(FY)
  if (Kx == 0 || Ky == 0) stop("After dropping near-zero-variance columns (threshold = 1e-10), no usable factors remain in X or Y.")
  
  zscore <- function(M) {
    mu <- colMeans(M, na.rm = TRUE)
    sdv <- apply(M, 2, sd, na.rm = TRUE); sdv[sdv < 1e-8 | !is.finite(sdv)] <- 1
    sweep(sweep(M, 2, mu, "-"), 2, sdv, "/")
  }
  FXz <- zscore(FX)
  FYz <- zscore(FY)
  
  Y1 <- FYz[(1 + lag):Tn, , drop = FALSE]
  Y0 <- FYz[1:(Tn - lag), , drop = FALSE]
  X0 <- FXz[1:(Tn - lag), , drop = FALSE]
  T_eff <- nrow(Y1)
  
  solve_ridge <- function(A, b, lambda = ridge) {
    p <- ncol(A)
    XtX <- crossprod(A)
    B <- solve(XtX + diag(lambda, p), crossprod(A, b))
    return(B)
  }
  
  Br <- solve_ridge(Y0, Y1)
  E_r <- Y1 - Y0 %*% Br
  SSE_r <- sum(E_r^2)
  
  D  <- cbind(Y0, X0)
  Bf <- solve_ridge(D, Y1)
  E_f <- Y1 - D %*% Bf
  SSE_f <- sum(E_f^2)
  
  SST  <- sum(Y1^2)
  
  R2_r_global <- 1 - SSE_r / SST
  R2_f_global <- 1 - SSE_f / SST
  delta_R2_in <- R2_f_global - R2_r_global
  partial_R2  <- 1 - SSE_f / SSE_r
  
  per_eq <- data.frame(
    eq = paste0("Y", seq_len(Ky)),
    R2_restricted = NA_real_, R2_full = NA_real_,
    delta_R2 = NA_real_, partial_R2 = NA_real_
  )
  for (j in seq_len(Ky)) {
    sst_j <- sum(Y1[, j]^2)
    sse_rj <- sum(E_r[, j]^2)
    sse_fj <- sum(E_f[, j]^2)
    per_eq$R2_restricted[j] <- 1 - sse_rj / sst_j
    per_eq$R2_full[j]       <- 1 - sse_fj / sst_j
    per_eq$delta_R2[j]      <- per_eq$R2_full[j] - per_eq$R2_restricted[j]
    per_eq$partial_R2[j]    <- 1 - sse_fj / sse_rj
  }
  
  OOS <- NULL
  if (oos) {
    set.seed(seed)
    t0 <- max(20, floor(T_eff * 0.5))
    errs_r <- matrix(NA_real_, T_eff - t0, Ky)
    errs_f <- matrix(NA_real_, T_eff - t0, Ky)
    dpred  <- matrix(NA_real_, T_eff - t0, Ky)   # yhat_restricted - yhat_full
    y_act  <- matrix(NA_real_, T_eff - t0, Ky)

    for (t in t0:(T_eff - 1)) {
      Y1_tr <- Y1[1:t, , drop = FALSE]
      Y0_tr <- Y0[1:t, , drop = FALSE]
      X0_tr <- X0[1:t, , drop = FALSE]

      Br_t <- solve_ridge(Y0_tr, Y1_tr)
      Bf_t <- solve_ridge(cbind(Y0_tr, X0_tr), Y1_tr)

      y_next <- Y1[t + 1, , drop = FALSE]
      y0_next <- Y0[t + 1, , drop = FALSE]
      x0_next <- X0[t + 1, , drop = FALSE]

      yhat_r <- y0_next %*% Br_t
      yhat_f <- cbind(y0_next, x0_next) %*% Bf_t

      idx <- t - t0 + 1
      errs_r[idx, ] <- y_next - yhat_r
      errs_f[idx, ] <- y_next - yhat_f
      dpred[idx, ]  <- yhat_r - yhat_f
      y_act[idx, ]  <- y_next
    }

    valid <- complete.cases(errs_r) & complete.cases(errs_f) & complete.cases(y_act)
    if (any(valid)) {
      vr <- errs_r[valid, , drop = FALSE]
      vf <- errs_f[valid, , drop = FALSE]
      dp <- dpred[valid, , drop = FALSE]
      ya <- y_act[valid, , drop = FALSE]

      rmse_r <- sqrt(mean(vr^2))
      rmse_f <- sqrt(mean(vf^2))
      sse_r <- sum(vr^2); sse_f <- sum(vf^2)
      ya_center <- ya - matrix(colMeans(ya), nrow(ya), ncol(ya), byrow = TRUE)
      sst_oos <- sum(ya_center^2)

      R2_oos_r <- 1 - sse_r / sst_oos
      R2_oos_f <- 1 - sse_f / sst_oos

      # Clark-West test of the nested OOS improvement (per-time loss summed over
      # equations); gives an inferential statement, not just a RMSE delta.
      cw <- hac_mean_test(rowSums(vr^2) - rowSums(vf^2) + rowSums(dp^2))

      OOS <- list(
        rmse_y  = rmse_r,
        rmse_yx = rmse_f,
        delta_rmse_pct = 100 * (rmse_r - rmse_f) / rmse_r,
        r2_oos_y  = R2_oos_r,
        r2_oos_yx = R2_oos_f,
        delta_r2_oos = R2_oos_f - R2_oos_r,
        cw_stat = cw$stat,
        cw_p = cw$p
      )
    }
  }
  
  if (verbose) {
    message("\n--- R2 OU (in-sample) ---")
    message(sprintf("Restricted R^2 (Y_{t-1}): %.3f", R2_r_global))
    message(sprintf("Full R2 (Y_{t-1}, X_{t-1}): %.3f", R2_f_global))
    message(sprintf("Delta R^2 (in-sample) = %.3f | Partial R^2 (X|Y) = %.3f",
                    delta_R2_in, partial_R2))
    if (!is.null(OOS)) {
      message("\n---OU R2 (out-of-sample, 1 step) ---")
      message(sprintf("RMSE Y-only: %.3f | RMSE Y+X: %.3f | ?RMSE: %.1f%%",
                      OOS$rmse_y, OOS$rmse_yx, OOS$delta_rmse_pct))
      message(sprintf("R2_OOS (Y-only): %.3f | R2_OOS (Y+X): %.3f | Delta R2_OOS: %.3f",
                      OOS$r2_oos_y, OOS$r2_oos_yx, OOS$delta_r2_oos))
    }
  }
  
  list(
    in_sample = list(
      R2_restricted = R2_r_global,
      R2_full = R2_f_global,
      delta_R2 = delta_R2_in,
      partial_R2 = partial_R2
    ),
    per_equation = per_eq,
    OOS = OOS,
    details = list(T_eff = T_eff, Kx = Kx, Ky = Ky, lag = lag)
  )
}

#' Plot error correction panel
#'
#' Creates a two-panel plot showing (1) error correction terms over time and
#' (2) their estimated half-lives, providing visual assessment of convergence speed.
#'
#' @param results_robust List. Results from main analysis pipeline.
#' @param use_contemporaneous_X Logical. Use contemporaneous X (time t, the
#'   default) or lagged X (time t-1) for computing the equilibrium path.
#' @param main_left Expression/character. Title for the left panel.
#' @param main_right Character string. Title for the right panel.
#'
#' @return Invisibly returns a list with components:
#'   \describe{
#'     \item{\code{U}}{Matrix of error correction terms (T x Ky).}
#'     \item{\code{half_lives}}{Vector of estimated half-lives for each Y factor.}
#'   }
#'
#' @export

plot_error_correction_panel <- function(results_robust,
                                        use_contemporaneous_X = TRUE,
                                        main_left = expression(u[t] == Y[t] - hat(mu)[Y](X[t])),
                                        main_right = "Half-life de u[t]") {
  stopifnot(!is.null(results_robust$factor_ou), !is.null(results_robust$factors))
  
  ou <- results_robust$factor_ou
  X  <- as.matrix(results_robust$factors$scores_X)
  Y  <- as.matrix(results_robust$factors$scores_Y)
  
  Tn <- min(nrow(X), nrow(Y))
  X  <- X[1:Tn, , drop = FALSE]
  Y  <- Y[1:Tn, , drop = FALSE]
  
  Kx <- ncol(X); Ky <- ncol(Y)
  
  if (is.null(ou$beta) || any(dim(ou$beta) != c(Ky, Kx))) {
    stop("No encuentro una matriz beta Ky x Kx en results_robust$factor_ou$beta.")
  }
  
  mu_x <- if (!is.null(ou$mu_x)) as.numeric(ou$mu_x)[1:Kx] else rep(0, Kx)
  mu_y <- if (!is.null(ou$mu_y)) as.numeric(ou$mu_y)[1:Ky] else rep(0, Ky)
  B    <- ou$beta
  
  if (use_contemporaneous_X) {
    Y_star <- matrix(mu_y, Tn, Ky, byrow = TRUE) +
      (X - matrix(mu_x, Tn, Kx, byrow = TRUE)) %*% t(B)
  } else {
    Y_star <- matrix(mu_y, Tn, Ky, byrow = TRUE)
    if (Tn >= 2) {
      Y_star[2:Tn, ] <- Y_star[2:Tn, ] +
        (X[1:(Tn-1), , drop = FALSE] - matrix(mu_x, Tn-1, Kx, byrow = TRUE)) %*% t(B)
    }
  }
  
  U <- Y - Y_star
  
  hl_fun <- if (exists("safe_ar1_half_life", mode = "function")) {
    get("safe_ar1_half_life")
  } else {
    function(x) {
      x <- as.numeric(x); x <- x[is.finite(x)]
      if (length(x) < 3) return(NA_real_)
      x <- x - mean(x)
      z <- x[-length(x)]; y <- x[-1]
      phi <- sum(y * z) / sum(z * z)
      phi <- max(min(phi, 0.9999), 1e-6)
      if (phi <= 0) return(Inf)
      log(0.5) / log(phi)
    }
  }
  hl_u <- apply(U, 2, hl_fun)
  
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par), add = TRUE)
  par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))
  
  matplot(U, type = "l", lty = 1, lwd = 1.4,
          xlab = "Tiempo", ylab = expression(u[t]), main = main_left)
  abline(h = 0, lty = 2, col = "grey50")
  legend("topright",
         legend = paste0("u", seq_len(Ky)),
         lty = 1, col = seq_len(Ky), bty = "n", cex = 0.9)
  
  ymax <- max(hl_u[is.finite(hl_u)], na.rm = TRUE)
  bp <- barplot(hl_u, names.arg = paste0("u", seq_len(Ky)),
                ylim = c(0, ymax * 1.15), ylab = "Periods",
                main = main_right)
  abline(h = stats::median(hl_u[is.finite(hl_u)]), lty = 2, col = "red")
  text(bp, hl_u + 0.05 * ymax, labels = round(hl_u, 1), cex = 0.8)
  mtext(side = 3, line = 0.2,
        text = paste0("Mediana HL: ",
                      round(stats::median(hl_u[is.finite(hl_u)]), 2),
                      " periods"),
        col = "red", cex = 0.9)
  
  invisible(list(U = U, half_lives = hl_u))
}

#' Extract X innovations (reduced-form VAR residuals)
#'
#' Computes \emph{reduced-form} innovations for the X factors by fitting a VAR and
#' returning its residuals. These are not structural shocks: identifying
#' structural shocks would require an identification scheme (e.g. a Cholesky or
#' sign restriction); the residuals here are correlated across equations.
#'
#' @param results List. Output from main analysis.
#' @param p Integer. VAR lag order. If \code{NULL}, uses order from DFM estimation.
#'
#' @return List with components:
#'   \describe{
#'     \item{\code{shocks}}{Matrix (T x Kx) of reduced-form innovations with NA
#'       padding for the initial lags.}
#'     \item{\code{var_fit}}{Fitted VAR model object.}
#'     \item{\code{lag_used}}{Integer lag order used.}
#'   }
#'
#' @export

make_X_innovations <- function(results, p = NULL) {
  stopifnot(!is.null(results$factors$scores_X))
  FX <- as.matrix(results$factors$scores_X)
  colnames(FX) <- make.names(colnames(FX), unique = TRUE)
  
  if (is.null(p)) p <- if (!is.null(results$dfm$p_used)) results$dfm$p_used else 1
  
  FXc <- scale(FX, scale = FALSE)
  
  fitX <- vars::VAR(as.data.frame(FXc), p = p, type = "const")
  E <- residuals(fitX)
  E <- rbind(matrix(NA_real_, p, ncol(E)), E)
  colnames(E) <- colnames(FX)
  
  list(shocks = E, var_fit = fitX, lag_used = p)
}