# Recovery tests that lock in the methodological fixes of this release.
# They are deliberately Stan-free so they run on any machine.

make_coupled <- function(T = 200, k = 2, beta = 0.9, seed = 42) {
  set.seed(seed)
  X <- matrix(0, T, k)
  for (t in 2:T) X[t, ] <- 0.6 * X[t - 1, ] + rnorm(k)
  Y <- matrix(0, T, k)
  for (t in 2:T) Y[t, ] <- 0.5 * Y[t - 1, ] + beta * X[t - 1, ] + rnorm(k, 0, 0.3)
  list(X = X, Y = Y)
}

fb_coupling <- function(X, Y) {
  estimate_factor_OU_fallback(list(scores_X = X, scores_Y = Y),
                              verbose = FALSE)$coupling_strength
}

test_that("circ_shift_rows wraps rows correctly", {
  M <- matrix(1:6, ncol = 1)            # 6 x 1
  expect_equal(as.numeric(circ_shift_rows(M, 1)), c(2, 3, 4, 5, 6, 1))
  expect_equal(as.numeric(circ_shift_rows(M, 6)), c(2, 3, 4, 5, 6, 1)) # tau=0 -> 1
})

test_that("coupling strength is INVARIANT to a column permutation of Y (why the old null was vacuous)", {
  d <- make_coupled()
  base  <- fb_coupling(d$X, d$Y)
  cperm <- fb_coupling(d$X, d$Y[, c(2, 1)])
  expect_equal(cperm, base, tolerance = 1e-8)
})

test_that("coupling strength DOES change under a time shift of Y (why the new null is valid)", {
  d <- make_coupled()
  base   <- fb_coupling(d$X, d$Y)
  cshift <- fb_coupling(d$X, circ_shift_rows(d$Y, 53))
  expect_false(isTRUE(all.equal(cshift, base, tolerance = 1e-3)))
})

test_that("rotation_null_test: the rotation null is degenerate (invariance diagnostic)", {
  d <- make_coupled()
  out <- rotation_null_test(d$X, d$Y, lag = 1, B = 120, seed = 1,
                            null_method = "rotation", progress = FALSE)
  # All coupling statistics are rotation-invariant: the null collapses to a point
  # (range ~ machine epsilon), so the rotation null can never reject. The exact
  # Monte Carlo p-value is governed by floating-point tie-breaking, not by a real
  # distribution; the rigorous claim is degeneracy + non-significance.
  expect_lt(diff(range(out$null_stats$procrustes_R)), 1e-6)
  expect_lt(diff(range(out$null_stats$dynbeta)),      1e-6)
  expect_gt(out$p_values$procrustes_R, 0.1)   # not significant
})

test_that("rotation_null_test: the time-shift null is non-degenerate and detects coupling", {
  d <- make_coupled(beta = 0.9)
  out <- rotation_null_test(d$X, d$Y, lag = 1, B = 300, seed = 1,
                            null_method = "circular_shift", progress = FALSE)
  expect_gt(stats::sd(out$null_stats$dynbeta), 0)        # genuine variation
  expect_lt(out$p_values$dynbeta, 0.05)                  # detects real coupling
  expect_true(all(c("dynbeta") %in% names(out$p_values_fdr)))
})

test_that("rotation_null_test honours its seed (reproducibility)", {
  d <- make_coupled()
  a <- rotation_null_test(d$X, d$Y, B = 100, seed = 99,
                          null_method = "circular_shift", progress = FALSE)
  b <- rotation_null_test(d$X, d$Y, B = 100, seed = 99,
                          null_method = "circular_shift", progress = FALSE)
  expect_equal(a$p_values$dynbeta, b$p_values$dynbeta)
  expect_equal(a$observed_vec, b$observed_vec)
})

test_that("convergence index is NOT tautological: phi credible bands above 1 lower the score", {
  # Factor 2 of X has a credible interval reaching above 1 (non-convergent).
  post <- list(
    phi_x = rbind(median = c(0.5, 0.5),
                  lower  = c(0.3, 0.3),
                  upper  = c(0.9, 1.2)),
    phi_y = rbind(median = c(0.4),
                  lower  = c(0.2),
                  upper  = c(0.8))
  )
  idx <- calculate_continuous_indices(list(posteriors = post))
  # band_frac(X) = 0.5 (1 of 2 in band), band_frac(Y) = 1 -> mean = 0.75
  expect_equal(idx$convergence$scores$posteriors, 0.75)
})

test_that("dm/cw test flags a genuinely better forecast", {
  set.seed(8)
  n <- 200
  e_full <- rnorm(n, 0, 1)
  e_base <- e_full + rnorm(n, 0, 1.2)   # baseline strictly worse on average
  res <- dm_cw_test(e_base, e_full, nested = FALSE)
  expect_true(is.finite(res$stat))
  expect_lt(res$p, 0.05)
})
