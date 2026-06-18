# Wedge diagnostics (transformation problem, D-S12.3 / D-S12.6)

make_price_system <- function(Tn = 30, S = 6, seed = 1) {
  set.seed(seed)
  k  <- matrix(runif(Tn * S, 1, 3), Tn, S)
  K  <- k + matrix(runif(Tn * S, 1, 4), Tn, S)
  p  <- matrix(runif(Tn * S, 0.1, 1), Tn, S)
  Gp <- rowSums(p) / rowSums(K)                       # G' = sum p / sum K (canonical)
  list(k = k, K = K, p = p, Gp = Gp)
}

test_that("compute_wedge sums to zero across sectors by construction", {
  s <- make_price_system()
  w <- compute_wedge(s$K, s$Gp, s$p)
  expect_lt(w$max_abs_rowsum_rel, 1e-10)
  expect_equal(dim(w$W), dim(s$p))
  # wedge = markup - p, and markup row totals equal surplus row totals
  expect_equal(rowSums(w$markup), rowSums(s$p), tolerance = 1e-10)
})

test_that("compute_wedge rejects an inconsistent G' (non zero-sum)", {
  s <- make_price_system()
  expect_error(compute_wedge(s$K, s$Gp * 1.5, s$p), "sum to zero")
  expect_error(compute_wedge(s$K, s$Gp[-1], s$p), "length T")
  expect_error(compute_wedge(s$K, s$Gp, s$p[, -1]), "same dimensions")
})

test_that("wedge_stationarity recovers reversion on a stationary panel", {
  skip_if_not_installed("urca")
  set.seed(7); Tn <- 80; S <- 5
  W <- sapply(1:S, function(i) as.numeric(stats::arima.sim(list(ar = 0.5), Tn)))
  W <- W - rowMeans(W)                                # impose cross-sectional zero-sum
  out <- wedge_stationarity(W, verbose = FALSE)
  expect_s3_class(out$per_sector, "data.frame")
  expect_equal(nrow(out$per_sector), S)
  expect_true(out$panel$frac_reject >= 0.6)           # AR(0.5) is clearly stationary
  expect_true(all(out$per_sector$rho < 1, na.rm = TRUE))
})

test_that("wedge_stationarity flags a random walk as non-reverting", {
  skip_if_not_installed("urca")
  set.seed(11); Tn <- 80; S <- 4
  W <- sapply(1:S, function(i) cumsum(rnorm(Tn)))     # unit-root series
  W <- W - rowMeans(W)
  out <- wedge_stationarity(W, verbose = FALSE)
  expect_true(out$panel$frac_reject <= 0.5)           # should mostly NOT reject
})

test_that("placebo 'uniform' makes the wedge vanish and V equal Phi", {
  s <- make_price_system()
  Phi <- s$k + sweep(s$K, 1, s$Gp, `*`)               # production price
  pl <- placebo_values(s$k, s$K, s$Gp, s$p, scheme = "uniform")
  expect_lt(max(abs(pl$W_placebo)), 1e-10)            # wedge identically zero
  expect_equal(pl$V_placebo, Phi, tolerance = 1e-10)  # value == production price
  expect_equal(rowSums(pl$p_placebo), rowSums(s$p), tolerance = 1e-10) # total preserved
})

test_that("placebo 'permute' preserves per-period totals and is a derangement", {
  s <- make_price_system()
  pl <- placebo_values(s$k, s$K, s$Gp, s$p, scheme = "permute", seed = 42)
  expect_equal(rowSums(pl$p_placebo), rowSums(s$p), tolerance = 1e-10)
  expect_true(all(pl$perm != seq_along(pl$perm)))     # no fixed point
  # reproducible with the same seed, and does NOT touch the global RNG state
  set.seed(123); before <- runif(1)
  pl2 <- placebo_values(s$k, s$K, s$Gp, s$p, scheme = "permute", seed = 42)
  after <- runif(1)
  expect_identical(pl$perm, pl2$perm)
  set.seed(123); expect_equal(before, runif(1))       # global stream unaffected
})
