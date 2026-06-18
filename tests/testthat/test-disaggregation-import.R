# Import of the canonical disaggregation (BayesianDisaggregation) and removal of
# the local deterministic-blend duplicate. All Stan-free: the conjugate engine
# and the fallback factor-OU are pure R, so these run on any machine.

test_that("read_cpi (relocated, hardened) reads a clean CPI fixture", {
  f <- test_path("fixtures", "cpi_demo.xlsx")
  skip_if_not(file.exists(f))
  df <- read_cpi(f)
  expect_s3_class(df, "data.frame")
  expect_named(df, c("Year", "CPI"))
  expect_type(df$Year, "integer")
  expect_true(is.numeric(df$CPI))
  expect_false(is.unsorted(df$Year))
  expect_true(all(df$CPI > 0))
})

test_that("read_cpi errors cleanly on unidentifiable columns and missing files", {
  f <- test_path("fixtures", "bad_cols.xlsx")
  skip_if_not(file.exists(f))
  expect_error(read_cpi(f), "No date or CPI columns")
  expect_error(read_cpi(tempfile(fileext = ".xlsx")), "not found")
})

test_that("read_weights_matrix returns a row-normalized prior matrix", {
  f <- test_path("fixtures", "weights_demo.xlsx")
  skip_if_not(file.exists(f))
  w <- convergenceDFM:::read_weights_matrix(f)
  expect_true(is.matrix(w$P))
  expect_equal(unname(rowSums(w$P)), rep(1, nrow(w$P)), tolerance = 1e-10)
  expect_length(w$industries, ncol(w$P))
  expect_type(w$years, "integer")
})

test_that("test_reweighting_robustness runs via the canonical conjugate engine", {
  fc <- test_path("fixtures", "cpi_demo.xlsx")
  fw <- test_path("fixtures", "weights_demo.xlsx")
  skip_if_not(file.exists(fc) && file.exists(fw))
  g <- readRDS(test_path("fixtures", "golden_reweight.rds"))

  rw <- test_reweighting_robustness(fc, fw, g$X_matrix,
                                    max_comp = g$max_comp, seed = g$seed,
                                    verbose = FALSE)
  expect_type(rw, "list")
  expect_true(all(c("results", "cv_coupling", "robust") %in% names(rw)))
  coup <- vapply(rw$results, function(x) x$coupling, numeric(1))
  expect_true(all(is.finite(coup)))
  expect_true(is.logical(rw$robust))
})

test_that("the deterministic-blend duplicate has been removed", {
  ns <- asNamespace("convergenceDFM")
  for (fn in c("run_disaggregation_custom_prior", "compute_L_from_P",
               "spread_likelihood", "posterior_weighted")) {
    expect_false(exists(fn, envir = ns, inherits = FALSE),
                 info = paste("leftover duplicate:", fn))
  }
})
