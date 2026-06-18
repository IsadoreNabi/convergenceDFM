# =====================================================================
# Wedge diagnostics for the transformation problem
# (OU jerarquico de precios de produccion, decision D-S12.3 / D-S12.6).
#
# The "wedge" W_i = Phi_i - V_i = K_i * G' - p_i isolates the marxian
# redistribution of surplus (production price minus value), CANCELLING the
# shared cost price k = c + v. It sums to zero across sectors by construction
# (G' = sum p / sum K), so it is orthogonal to the level co-movement that a
# "spurious correlation" critique targets. The empirical content of gravitation
# lives in the DYNAMICS of the wedge (boundedness / mean reversion), not in the
# definitional level co-movement. These functions build the wedge, verify its
# cross-sectional zero-sum, test its per-sector and panel stationarity, and
# construct aggregate-preserving placebo values as negative controls.
# =====================================================================

#' Construct the transformation wedge W = Phi - V = K * G' - p
#'
#' Builds the sector-by-time wedge that isolates the marxian redistribution of
#' surplus from the shared cost price. The wedge cancels the cost price
#' \eqn{k = c + v} that both the production price \eqn{\Phi = k + K G'} and the
#' value \eqn{V = k + p} share, leaving \eqn{W_i = K_i G' - p_i}. By
#' construction of the general profit rate \eqn{G' = \sum_i p_i / \sum_i K_i},
#' the wedge sums to zero across sectors each period (redistribution, not
#' creation of value: the marxian invariance).
#'
#' @param K Numeric matrix \code{[T, S]} of total advanced capital per sector.
#' @param Gprime Numeric vector of length \code{T}, the general profit rate per
#'   period (typically \code{rowSums(p) / rowSums(K)}).
#' @param p Numeric matrix \code{[T, S]} of sectoral surplus (plusvalia, EBO).
#' @param tol Relative tolerance for the per-period zero-sum check. Default
#'   \code{1e-8}.
#' @param check_zero_sum Logical; if \code{TRUE} (default) verify
#'   \eqn{\sum_i W_{i,t} \approx 0} for every \code{t} and error out otherwise.
#'
#' @return A list with: \code{W} (the \code{[T, S]} wedge), \code{row_sums}
#'   (per-period cross-sectional sums, expected ~0), \code{max_abs_rowsum_rel}
#'   (the largest relative zero-sum violation), and \code{markup} (the profit
#'   markup \eqn{K_i G'}).
#'
#' @examples
#' set.seed(1)
#' Tn <- 20; S <- 5
#' K <- matrix(runif(Tn * S, 1, 5), Tn, S)
#' p <- matrix(runif(Tn * S, 0.1, 1), Tn, S)
#' Gp <- rowSums(p) / rowSums(K)
#' w <- compute_wedge(K, Gp, p)
#' max(abs(w$row_sums))  # ~ 0 by construction
#' @export
compute_wedge <- function(K, Gprime, p, tol = 1e-8, check_zero_sum = TRUE) {
  K <- as.matrix(K); p <- as.matrix(p)
  if (!all(dim(K) == dim(p)))
    stop("K and p must have the same dimensions [T, S].", call. = FALSE)
  if (length(Gprime) != nrow(K))
    stop("Gprime must have length T = nrow(K).", call. = FALSE)
  markup <- sweep(K, 1, Gprime, `*`)          # K_i * G'_t
  W <- markup - p
  row_sums <- rowSums(W)
  scale_t  <- rowSums(abs(markup)) + .Machine$double.eps
  max_abs_rowsum_rel <- max(abs(row_sums) / scale_t)
  if (isTRUE(check_zero_sum) && max_abs_rowsum_rel > tol)
    stop("Wedge does not sum to zero across sectors (max rel. violation ",
         signif(max_abs_rowsum_rel, 3), " > tol ", tol,
         "). Check that Gprime = rowSums(p)/rowSums(K).", call. = FALSE)
  list(W = W, row_sums = row_sums,
       max_abs_rowsum_rel = max_abs_rowsum_rel, markup = markup)
}

#' Per-sector and panel stationarity / mean-reversion of the wedge
#'
#' Tests whether the transformation wedge is a bounded, mean-reverting series
#' (evidence that the redistribution gravitates and does not diverge) rather
#' than a random walk. For each sector it fits a discrete AR(1)
#' \eqn{W_t = c + \rho W_{t-1} + e_t} (reporting \eqn{\rho}, the implied
#' mean-reversion speed \eqn{\kappa = -\log\rho} and half-life
#' \eqn{\log(0.5)/\log\rho}) and an augmented Dickey-Fuller test
#' (\code{urca::ur.df}, drift). The panel summary reports the fraction of
#' sectors rejecting the unit-root null and the distribution of reversion
#' speeds. Boundedness is reported as the range and the max-abs over standard
#' deviation. This is a descriptive diagnostic; no single pooled p-value is
#' claimed (anti-overreach: panel unit-root nulls require their own tabulated
#' distribution).
#'
#' @param W Numeric matrix \code{[T, S]}, the wedge (e.g. from
#'   \code{compute_wedge}).
#' @param adf_lags Maximum lag order for the ADF auxiliary regression
#'   (\code{urca::ur.df} \code{lags=}, with \code{selectlags = "AIC"}). Default
#'   \code{2}.
#' @param alpha Significance level for the ADF unit-root rejection. Default
#'   \code{0.05} (uses the 5pct critical value; \code{0.01}/\code{0.10}
#'   supported).
#' @param verbose Logical; print a short panel summary. Default \code{TRUE}.
#'
#' @return A list with \code{per_sector} (a data frame: \code{rho},
#'   \code{kappa}, \code{half_life}, \code{adf_stat}, \code{adf_crit},
#'   \code{reject_unit_root}, \code{range}, \code{maxabs_over_sd}) and
#'   \code{panel} (fraction rejecting, median/IQR of \code{rho} and
#'   \code{half_life}, and the boundedness summary).
#'
#' @examples
#' \donttest{
#' if (requireNamespace("urca", quietly = TRUE)) {
#'   set.seed(1); Tn <- 60; S <- 4
#'   W <- sapply(1:S, function(s) as.numeric(stats::arima.sim(list(ar = 0.6), Tn)))
#'   W <- W - rowMeans(W)  # impose zero-sum like a real wedge
#'   wedge_stationarity(W, verbose = FALSE)$panel$frac_reject
#' }
#' }
#' @export
wedge_stationarity <- function(W, adf_lags = 2, alpha = 0.05, verbose = TRUE) {
  if (!requireNamespace("urca", quietly = TRUE))
    stop("Package 'urca' is required for the ADF unit-root test.", call. = FALSE)
  W <- as.matrix(W); S <- ncol(W); Tn <- nrow(W)
  if (Tn < 10L) warning("Short series (T < 10): AR(1)/ADF estimates are noisy.", call. = FALSE)
  cl <- if (alpha <= 0.01) "1pct" else if (alpha >= 0.10) "10pct" else "5pct"
  ps <- lapply(seq_len(S), function(s) {
    w <- W[, s]
    ar1 <- stats::lm(w[-1] ~ w[-Tn])                 # W_t ~ c + rho W_{t-1}
    rho <- unname(stats::coef(ar1)[2])
    hl  <- if (!is.na(rho) && rho > 0 && rho < 1) log(0.5) / log(rho) else NA_real_
    kap <- if (!is.na(rho) && rho > 0) -log(rho) else NA_real_
    df  <- tryCatch(urca::ur.df(w, type = "drift", lags = adf_lags, selectlags = "AIC"),
                    error = function(e) NULL)
    if (is.null(df)) { stat <- NA_real_; crit <- NA_real_; rej <- NA }
    else {
      stat <- methods::slot(df, "teststat")[1]        # tau (unit-root) statistic
      cvm  <- methods::slot(df, "cval"); crit <- cvm[1, cl]
      rej  <- stat < crit                             # reject unit root if more negative
    }
    data.frame(sector = s, rho = rho, kappa = kap, half_life = hl,
               adf_stat = stat, adf_crit = crit, reject_unit_root = rej,
               range = diff(range(w)),
               maxabs_over_sd = max(abs(w - mean(w))) / (stats::sd(w) + .Machine$double.eps))
  })
  per_sector <- do.call(rbind, ps)
  if (!is.null(colnames(W))) per_sector$sector_name <- colnames(W)
  rej <- per_sector$reject_unit_root
  panel <- list(
    n_sectors = S, T = Tn,
    frac_reject = mean(rej, na.rm = TRUE),
    n_reject = sum(rej, na.rm = TRUE),
    rho_median = stats::median(per_sector$rho, na.rm = TRUE),
    half_life_median = stats::median(per_sector$half_life, na.rm = TRUE),
    half_life_iqr = stats::IQR(per_sector$half_life, na.rm = TRUE),
    maxabs_over_sd_median = stats::median(per_sector$maxabs_over_sd, na.rm = TRUE)
  )
  if (isTRUE(verbose)) {
    message(sprintf("Wedge stationarity: %d/%d sectors reject unit root at %.0f%% (%.0f%%);",
                    panel$n_reject, S, 100 * alpha, 100 * panel$frac_reject))
    message(sprintf("  median rho = %.3f, median half-life = %.1f periods.",
                    panel$rho_median, panel$half_life_median))
  }
  list(per_sector = per_sector, panel = panel)
}

#' Aggregate-preserving placebo values (negative controls for gravitation)
#'
#' Builds theoretically false "values" that respect the SAME aggregate marxian
#' constraints (per-period total surplus \eqn{\sum_i p_i} and the general profit
#' rate \eqn{G'}) but BREAK the marxian sectoral structure. Used as negative
#' controls: if the genuine sectoral value gravitates (wedge reverts / OOS)
#' better than placebos that share all aggregate accounting, the specific
#' marxian assignment carries information beyond the shared cost price.
#'
#' Schemes:
#' \describe{
#'   \item{\code{"uniform"}}{Every sector earns the average profit rate:
#'     \eqn{p^{plac}_i = G' K_i}, so \eqn{V^{plac} = k + G' K = \Phi} and the
#'     wedge is identically zero (the "no redistribution" null). Tests whether
#'     the value structure adds anything over the production price.}
#'   \item{\code{"permute"}}{A fixed (seeded) derangement \eqn{\sigma} reassigns
#'     each sector's surplus series to another sector:
#'     \eqn{p^{plac}_i = p_{\sigma(i)}}. Per-period totals are preserved
#'     (\eqn{\sum_i p_{\sigma(i)} = \sum_i p_i}) but the \eqn{K_i \leftrightarrow p_i}
#'     pairing is destroyed, so the wedge \eqn{K_i G' - p_{\sigma(i)}} mis-pairs
#'     capital and surplus.}
#' }
#'
#' @param k_cost Numeric matrix \code{[T, S]}, the cost price \eqn{k = c + v}.
#' @param K Numeric matrix \code{[T, S]}, total advanced capital.
#' @param Gprime Numeric vector length \code{T}, the general profit rate.
#' @param p Numeric matrix \code{[T, S]}, the true sectoral surplus.
#' @param scheme Either \code{"uniform"} or \code{"permute"}.
#' @param seed Integer seed for the permutation (reproducibility). Default
#'   \code{1L}. Only used by \code{"permute"}.
#'
#' @return A list with \code{p_placebo}, \code{V_placebo} (\eqn{= k + p^{plac}}),
#'   \code{W_placebo} (\eqn{= K G' - p^{plac}}), the \code{scheme}, and (for
#'   \code{"permute"}) the \code{perm} used.
#'
#' @examples
#' set.seed(1); Tn <- 12; S <- 6
#' k <- matrix(runif(Tn * S, 1, 3), Tn, S)
#' K <- k + matrix(runif(Tn * S, 1, 3), Tn, S)
#' p <- matrix(runif(Tn * S, 0.1, 1), Tn, S)
#' Gp <- rowSums(p) / rowSums(K)
#' pl <- placebo_values(k, K, Gp, p, scheme = "permute", seed = 42)
#' max(abs(rowSums(pl$p_placebo) - rowSums(p)))  # ~0: per-period total preserved
#' @export
placebo_values <- function(k_cost, K, Gprime, p,
                           scheme = c("uniform", "permute"), seed = 1L) {
  scheme <- match.arg(scheme)
  k_cost <- as.matrix(k_cost); K <- as.matrix(K); p <- as.matrix(p)
  S <- ncol(p)
  if (scheme == "uniform") {
    p_placebo <- sweep(K, 1, Gprime, `*`)            # G' * K_i (uniform rate)
    perm <- NULL
  } else {
    if (S < 2L) stop("Permutation placebo needs at least 2 sectors.", call. = FALSE)
    perm <- .derangement(S, seed)
    p_placebo <- p[, perm, drop = FALSE]             # fixed reassignment across time
  }
  V_placebo <- k_cost + p_placebo
  W_placebo <- sweep(K, 1, Gprime, `*`) - p_placebo
  if (!is.null(colnames(p))) {
    colnames(p_placebo) <- colnames(V_placebo) <- colnames(W_placebo) <- colnames(p)
  }
  list(p_placebo = p_placebo, V_placebo = V_placebo, W_placebo = W_placebo,
       scheme = scheme, perm = perm)
}

# Seeded derangement (no fixed points): each sector maps to a DIFFERENT sector.
.derangement <- function(n, seed) {
  old <- .Random.seed_backup()
  on.exit(.Random.seed_restore(old), add = TRUE)
  set.seed(seed)
  repeat {
    perm <- sample.int(n)
    if (all(perm != seq_len(n))) return(perm)        # reject any fixed point
  }
}
.Random.seed_backup <- function()
  if (exists(".Random.seed", envir = .GlobalEnv)) get(".Random.seed", envir = .GlobalEnv) else NULL
.Random.seed_restore <- function(old) {
  if (is.null(old)) {
    if (exists(".Random.seed", envir = .GlobalEnv))
      rm(".Random.seed", envir = .GlobalEnv)
  } else assign(".Random.seed", old, envir = .GlobalEnv)
}
