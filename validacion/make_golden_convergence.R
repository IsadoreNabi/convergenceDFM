# =============================================================================
# Regenerate the golden fixtures for the 0.3.0 additions:
#   * Leave-Cluster-Out (test_leave_cluster_out): in-memory, bit-exact golden.
#   * Reweighting robustness via the canonical conjugate disaggregation
#     (test_reweighting_robustness): xlsx fixtures + tolerant numeric golden
#     (the xlsx round-trip is a serialization layer, so the reweighting golden
#     is a numerical regression guard, not an algebraic identity -- the rigor
#     layers are kept separate).
#
# Run from the package root:
#   Rscript validacion/make_golden_convergence.R
# Requires openxlsx only here (a dev-time dependency, NOT a package dependency);
# the tests themselves read the xlsx fixtures with readxl (already in Imports).
# =============================================================================

suppressMessages(pkgload::load_all(".", quiet = TRUE))
stopifnot(requireNamespace("openxlsx", quietly = TRUE))

fix_dir <- file.path("tests", "testthat", "fixtures")
dir.create(fix_dir, recursive = TRUE, showWarnings = FALSE)

# ---- Leave-Cluster-Out golden (deterministic, no IO) ------------------------
make_lco_data <- function() {
  set.seed(123)
  T <- 30; K <- 6
  f <- cumsum(rnorm(T))
  Phi <- sapply(1:K, function(k) 100 + 5 * f + rnorm(T, 0, 1))
  phi <- sapply(1:K, function(k) Phi[, k] + rnorm(T, 0, 0.5))
  colnames(Phi) <- colnames(phi) <- paste0("sector_", 1:K)
  list(Phi = Phi, phi = phi)
}
d <- make_lco_data()
cluster_map <- c(1, 1, 2, 2, 3, 3)          # three balanced "value chains"
lco <- test_leave_cluster_out(d$Phi, d$phi, cluster_map = cluster_map,
                              seed = 7, verbose = FALSE)
golden_lco <- list(
  Phi = d$Phi, phi = d$phi, cluster_map = cluster_map, seed = 7L,
  ref = list(
    baseline          = lco$baseline,
    cluster_estimates = lco$cluster_estimates,
    influence         = lco$influence,
    retention         = lco$retention,
    bias              = lco$bias,
    se                = lco$se,
    robust            = lco$robust
  )
)
saveRDS(golden_lco, file.path(fix_dir, "golden_lco.rds"), version = 2)
cat("Wrote golden_lco.rds  (baseline =", round(lco$baseline, 6), ")\n")

# ---- Reweighting golden + xlsx fixtures -------------------------------------
make_reweight_data <- function() {
  set.seed(2)
  T <- 25; K <- 5
  years <- 1999:(1999 + T - 1)
  cpi <- 100 * cumprod(1 + rnorm(T, 0.02, 0.01)) + 50
  W <- matrix(runif(T * K), T, K); W <- W / rowSums(W)
  ind <- paste0("sector_", 1:K)
  Xm <- sapply(1:K, function(k) 100 + 3 * cumsum(rnorm(T)) + rnorm(T))
  colnames(Xm) <- ind
  list(years = years, cpi = cpi, W = W, ind = ind, Xm = Xm)
}
r <- make_reweight_data()

cpi_df <- data.frame(Date = r$years, CPI = r$cpi)
wdf <- data.frame(Industry = r$ind, check.names = FALSE)
for (i in seq_along(r$years)) wdf[[as.character(r$years[i])]] <- r$W[i, ]

f_cpi <- file.path(fix_dir, "cpi_demo.xlsx")
f_w   <- file.path(fix_dir, "weights_demo.xlsx")
openxlsx::write.xlsx(cpi_df, f_cpi)
openxlsx::write.xlsx(wdf, f_w)

# A sheet with no identifiable date/CPI columns, to exercise read_cpi's guard.
openxlsx::write.xlsx(data.frame(foo = 1:3, bar = 4:6),
                     file.path(fix_dir, "bad_cols.xlsx"))

rw <- test_reweighting_robustness(f_cpi, f_w, r$Xm, max_comp = 3, seed = 11,
                                  verbose = FALSE)
golden_reweight <- list(
  X_matrix = r$Xm, max_comp = 3L, seed = 11L,
  ref = list(
    couplings   = vapply(rw$results, function(x) x$coupling, numeric(1)),
    cv_coupling = rw$cv_coupling,
    robust      = rw$robust
  )
)
saveRDS(golden_reweight, file.path(fix_dir, "golden_reweight.rds"), version = 2)
cat("Wrote golden_reweight.rds + cpi_demo.xlsx + weights_demo.xlsx\n")
cat("  couplings =", paste(round(golden_reweight$ref$couplings, 6), collapse = " "), "\n")
