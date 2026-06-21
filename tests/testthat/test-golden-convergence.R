# Golden (regression) tests for the 0.3.0 additions.
#
# Layer separation (the three rigor layers are kept strictly apart):
#
#   * Leave-Cluster-Out is a *deterministic* in-memory computation: on a fixed
#     platform it is bit-reproducible. However, the shared coupling pipeline
#     goes through floating-point linear algebra (BLAS), whose last ~3
#     significant digits are not portable across CRAN's heterogeneous
#     platforms / BLAS builds (observed inter-platform drift ~1e-11). The
#     portable invariant is therefore *numerical* (a tight tolerance); the
#     stronger *bit-exact* invariant is retained only as a developer-machine
#     guard, skipped on CRAN.
#
#   * Reweighting via the conjugate engine additionally reads the xlsx
#     fixtures, so the serialization round-trip plus the conjugate smoother
#     widen the inter-platform drift (~1e-7). It is a numerical regression
#     guard compared with a tolerance, never an algebraic identity.
#
# A genuine algorithmic regression would move the estimates by O(1e-2) or more,
# so these tolerances stay far below the regression scale while absorbing the
# non-portable floating-point noise.
#
# Regenerate both with: Rscript validacion/make_golden_convergence.R

# Cross-platform tolerance for the leave-cluster-out golden (five-order margin
# over the observed ~1e-11 linear-algebra drift).
.lco_tol <- 1e-6

test_that("test_leave_cluster_out reproduces the frozen golden (numeric, portable)", {
  g <- readRDS(test_path("fixtures", "golden_lco.rds"))
  lco <- test_leave_cluster_out(g$Phi, g$phi, cluster_map = g$cluster_map,
                                seed = g$seed, verbose = FALSE)
  expect_equal(lco$baseline,          g$ref$baseline,          tolerance = .lco_tol)
  expect_equal(lco$cluster_estimates, g$ref$cluster_estimates, tolerance = .lco_tol)
  expect_equal(lco$influence,         g$ref$influence,         tolerance = .lco_tol)
  expect_equal(lco$retention,         g$ref$retention,         tolerance = .lco_tol)
  expect_equal(lco$bias,              g$ref$bias,              tolerance = .lco_tol)
  expect_equal(lco$se,                g$ref$se,                tolerance = .lco_tol)
  expect_identical(lco$robust,        g$ref$robust)  # logical: platform-invariant
})

test_that("test_leave_cluster_out reproduces the frozen golden bit-exactly (dev machine)", {
  skip_on_cran()
  g <- readRDS(test_path("fixtures", "golden_lco.rds"))
  lco <- test_leave_cluster_out(g$Phi, g$phi, cluster_map = g$cluster_map,
                                seed = g$seed, verbose = FALSE)
  expect_identical(lco$baseline,          g$ref$baseline)
  expect_identical(lco$cluster_estimates, g$ref$cluster_estimates)
  expect_identical(lco$influence,         g$ref$influence)
  expect_identical(lco$retention,         g$ref$retention)
  expect_identical(lco$bias,              g$ref$bias)
  expect_identical(lco$se,                g$ref$se)
  expect_identical(lco$robust,            g$ref$robust)
})

test_that("test_reweighting_robustness reproduces the frozen golden (numeric)", {
  fc <- test_path("fixtures", "cpi_demo.xlsx")
  fw <- test_path("fixtures", "weights_demo.xlsx")
  skip_if_not(file.exists(fc) && file.exists(fw))
  g <- readRDS(test_path("fixtures", "golden_reweight.rds"))

  rw <- test_reweighting_robustness(fc, fw, g$X_matrix,
                                    max_comp = g$max_comp, seed = g$seed,
                                    verbose = FALSE)
  coup <- vapply(rw$results, function(x) x$coupling, numeric(1))
  # xlsx round-trip + conjugate smoother widen the drift to ~1e-7; 1e-5 is
  # portable yet two orders below the regression scale.
  expect_equal(coup, g$ref$couplings, tolerance = 1e-5)
  expect_equal(rw$cv_coupling, g$ref$cv_coupling, tolerance = 1e-5)
  expect_identical(rw$robust, g$ref$robust)
})
