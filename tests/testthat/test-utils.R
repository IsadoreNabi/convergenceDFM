test_that("to_num_commas handles European and Anglo-Saxon formats", {
  expect_equal(to_num_commas("1.234,56"), 1234.56)   # European
  expect_equal(to_num_commas("1,234.56"), 1234.56)   # US
  expect_equal(to_num_commas("1234.56"), 1234.56)    # plain decimal
  expect_equal(to_num_commas(1234.56), 1234.56)      # numeric passthrough
  expect_equal(to_num_commas("1.234.567"), 1234567)  # repeated dots = thousands
  expect_true(is.na(to_num_commas("abc")))           # unparseable -> NA
})

test_that("to_num_commas is vectorized", {
  expect_equal(to_num_commas(c("1.234,56", "1,234.56")), c(1234.56, 1234.56))
})

test_that("row_norm1 normalizes rows to sum 1", {
  M <- matrix(c(1, 2, 3, 4, 5, 6), nrow = 2)
  M_norm <- row_norm1(M)
  expect_equal(rowSums(M_norm), c(1, 1))
})

test_that("safe_div handles division by zero", {
  expect_true(is.finite(safe_div(1, 0)))
  expect_true(safe_div(0, 0) > 0)
})

test_that("geometric_mean_robust matches the closed form and guards empties", {
  x <- c(1, 4, 16)
  expect_equal(geometric_mean_robust(x), exp(mean(log(x))), tolerance = 1e-8)
  expect_true(is.na(geometric_mean_robust(numeric(0))))
})

test_that("logging threshold round-trips through .pkgenv", {
  old <- set_log_level("WARN")
  expect_equal(get_log_level(), "WARN")
  set_log_level(old)
  expect_equal(get_log_level(), old)
})
