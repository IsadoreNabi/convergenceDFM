# Golden (regression) tests for the 0.3.0 additions.
#
# Layer separation (do not conflate):
#   * Leave-Cluster-Out: a deterministic, in-memory computation -> bit-exact
#     golden via expect_identical.
#   * Reweighting via the conjugate engine: reads the xlsx fixtures, so the
#     serialization round-trip makes it a *numerical* regression guard, compared
#     with a tight tolerance rather than an algebraic identity.
# Regenerate both with: Rscript validacion/make_golden_convergence.R

test_that("test_leave_cluster_out reproduces the frozen golden bit-exactly", {
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
  expect_equal(coup, g$ref$couplings, tolerance = 1e-8)
  expect_equal(rw$cv_coupling, g$ref$cv_coupling, tolerance = 1e-8)
  expect_identical(rw$robust, g$ref$robust)
})
