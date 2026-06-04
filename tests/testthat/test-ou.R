test_that("ou_half_life is correct and flags non-convergence", {
  expect_equal(ou_half_life(0.5), log(0.5) / log(0.5))   # = 1
  expect_equal(ou_half_life(0.25), log(0.5) / log(0.25)) # = 0.5
  expect_identical(ou_half_life(1), Inf)                  # unit root -> Inf
  expect_identical(ou_half_life(1.2), Inf)               # explosive -> Inf
  expect_equal(ou_half_life(0), 0)
  expect_true(is.na(ou_half_life(NA_real_)))
})

test_that("safe_ar1_half_life recovers a persistent AR(1)", {
  set.seed(1)
  n <- 400; phi <- 0.7
  x <- numeric(n)
  for (t in 2:n) x[t] <- phi * x[t - 1] + rnorm(1)
  hl <- safe_ar1_half_life(x)
  expect_true(is.finite(hl) && hl > 0)
  expect_equal(hl, log(0.5) / log(phi), tolerance = 0.5)  # loose
})

test_that("fallback OU detects coupling when Y depends on lagged X", {
  set.seed(7)
  T <- 200; k <- 2
  X <- matrix(0, T, k)
  for (t in 2:T) X[t, ] <- 0.6 * X[t - 1, ] + rnorm(k)

  # Coupled Y
  Yc <- matrix(0, T, k)
  for (t in 2:T) Yc[t, ] <- 0.5 * Yc[t - 1, ] + 0.9 * X[t - 1, ] + rnorm(k, 0, 0.3)
  # Independent Y
  Yi <- matrix(0, T, k)
  for (t in 2:T) Yi[t, ] <- 0.5 * Yi[t - 1, ] + rnorm(k, 0, 0.3)

  cpl   <- estimate_factor_OU_fallback(list(scores_X = X, scores_Y = Yc),
                                       verbose = FALSE)$coupling_strength
  indep <- estimate_factor_OU_fallback(list(scores_X = X, scores_Y = Yi),
                                       verbose = FALSE)$coupling_strength
  expect_true(is.finite(cpl) && is.finite(indep))
  expect_gt(cpl, indep)
})

test_that("estimate_factor_OU falls back gracefully without Stan", {
  set.seed(3)
  fd <- list(scores_X = matrix(rnorm(60 * 2), 60, 2),
             scores_Y = matrix(rnorm(60 * 2), 60, 2))
  res <- suppressWarnings(estimate_factor_OU(fd, chains = 1, iter = 300,
                                             verbose = FALSE))
  expect_true(is.list(res))
  expect_true(!is.null(res$coupling_strength))
  expect_true(!is.null(res$half_lives_Y))
})
