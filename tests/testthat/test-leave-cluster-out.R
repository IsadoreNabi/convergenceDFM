# Leave-Cluster-Out and its pluggable sector->cluster map. Stan-free.

make_nested <- function(T = 30, K = 6, seed = 123) {
  set.seed(seed)
  f <- cumsum(rnorm(T))
  Phi <- sapply(1:K, function(k) 100 + 5 * f + rnorm(T, 0, 1))   # production
  phi <- sapply(1:K, function(k) Phi[, k] + rnorm(T, 0, 0.5))    # market (couples)
  colnames(Phi) <- colnames(phi) <- paste0("sector_", 1:K)
  list(Phi = Phi, phi = phi)
}

test_that("build_cluster_map (correlation) returns one labelled cluster per sector", {
  d <- make_nested()
  cm <- build_cluster_map(d$phi, n_clusters = 3)
  expect_length(cm, ncol(d$phi))
  expect_named(cm, colnames(d$phi))
  expect_true(all(cm %in% 1:3))
  expect_true(length(unique(cm)) <= 3)
})

test_that("build_cluster_map collapses to one cluster when K < 2", {
  one <- matrix(rnorm(20), 20, 1, dimnames = list(NULL, "only"))
  cm <- build_cluster_map(one, n_clusters = 3)
  expect_equal(unname(cm), 1L)
})

test_that("build_cluster_map (com) bins by organic composition and validates input", {
  d <- make_nested()
  com <- c(0.1, 0.2, 0.5, 0.6, 0.9, 1.0)
  cm <- build_cluster_map(d$phi, n_clusters = 2, method = "com", com = com)
  expect_length(cm, ncol(d$phi))
  expect_true(all(cm %in% 1:2))
  expect_error(build_cluster_map(d$phi, method = "com"), "requires .com.")
  expect_error(build_cluster_map(d$phi, method = "com", com = com[1:2]),
               "one value per sector")
})

test_that("test_leave_cluster_out (fallback) returns the documented structure", {
  d <- make_nested()
  lco <- test_leave_cluster_out(d$Phi, d$phi, n_clusters = 3, seed = 7,
                                verbose = FALSE)
  expect_true(all(c("baseline", "cluster_estimates", "cluster_members",
                    "cluster_sizes", "bias", "se", "influence", "retention",
                    "influential_clusters", "n_clusters", "cluster_map",
                    "map_source", "robust") %in% names(lco)))
  expect_match(lco$map_source, "^fallback:")
  expect_true(is.finite(lco$baseline))
  expect_equal(sum(lco$cluster_sizes), ncol(d$phi))
  expect_true(is.logical(lco$robust))
})

test_that("LCO baseline matches the jackknife baseline (shared coupling pipeline)", {
  d <- make_nested()
  lco <- test_leave_cluster_out(d$Phi, d$phi, cluster_map = c(1,1,2,2,3,3),
                                seed = 7, verbose = FALSE)
  jk <- test_jackknife_sectors(d$Phi, d$phi, seed = 7, verbose = FALSE)
  expect_identical(lco$baseline, jk$baseline)
})

test_that("test_leave_cluster_out accepts a user vector map and a list map", {
  d <- make_nested()
  v <- test_leave_cluster_out(d$Phi, d$phi, cluster_map = c(1,1,2,2,3,3),
                              seed = 7, verbose = FALSE)
  expect_identical(v$map_source, "user")
  expect_equal(v$n_clusters, 3)

  lst <- list(chainA = c("sector_1", "sector_2"),
              chainB = c("sector_3", "sector_4"),
              chainC = c("sector_5", "sector_6"))
  l <- test_leave_cluster_out(d$Phi, d$phi, cluster_map = lst, seed = 7,
                              verbose = FALSE)
  expect_identical(l$map_source, "user")
  expect_setequal(names(l$cluster_estimates), names(lst))
})

test_that("a cluster map that misses a sector errors", {
  d <- make_nested()
  bad <- list(chainA = c("sector_1", "sector_2"),
              chainB = c("sector_3", "sector_4"))   # leaves 5,6 unassigned
  expect_error(
    test_leave_cluster_out(d$Phi, d$phi, cluster_map = bad, verbose = FALSE),
    "does not cover every sector"
  )
})
