# Helper: a stationary VAR(1) in 4 dimensions.
make_var_factors <- function(T = 150, seed = 11) {
  set.seed(seed)
  k <- 4
  A <- diag(c(0.5, 0.4, 0.6, 0.3))
  Z <- matrix(0, T, k)
  for (t in 2:T) Z[t, ] <- A %*% Z[t - 1, ] + rnorm(k, 0, 1)
  list(scores_X = Z[, 1:2, drop = FALSE],
       scores_Y = Z[, 3:4, drop = FALSE])
}

test_that("calculate_hc_manual returns finite HC and HC4 SEs", {
  set.seed(2)
  n <- 80
  x <- rnorm(n); y <- 1 + 0.5 * x + rnorm(n) * (1 + abs(x))  # heteroskedastic
  m <- lm(y ~ x)
  for (ty in c("HC0", "HC1", "HC2", "HC3", "HC4")) {
    res <- calculate_hc_manual(m, type = ty)
    expect_true(isTRUE(res$success))
    expect_equal(length(res$se), length(coef(m)))
    expect_true(all(is.finite(res$se)) && all(res$se > 0))
  }
})

test_that("estimate_DFM fits a stable VAR and returns robust inference", {
  fd <- make_var_factors()
  res <- suppressWarnings(estimate_DFM(fd, p = 1, compute_oos = FALSE,
                                       verbose = FALSE))
  expect_true(res$is_stable)
  expect_true(is.finite(res$r2_global))
  expect_true(is.finite(res$half_life_dominant) || is.infinite(res$half_life_dominant))
  expect_true(is.list(res$robust_inference))
  # robust_inference carries coef/se/t/p, not just SEs:
  first_eq <- res$robust_inference[[1]]
  expect_true(all(c("coef", "se", "t", "p") %in% names(first_eq)))
})

test_that("choose_var_lag returns a valid, stable model", {
  fd <- make_var_factors()
  Fc <- cbind(fd$scores_X, fd$scores_Y)
  sel <- suppressWarnings(choose_var_lag(Fc, lag.max = 3, oos_eval = FALSE,
                                         verbose = FALSE))
  expect_true(sel$p >= 1)
  expect_true(inherits(sel$fit, "varest"))
})
