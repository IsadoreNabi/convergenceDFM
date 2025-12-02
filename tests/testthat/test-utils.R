test_that("to_num_commas handles European number format", {
  expect_equal(to_num_commas("1.234,56"), 1234.56)
  expect_equal(to_num_commas("1.234,56"), 1234.56)
  expect_equal(to_num_commas(1234.56), 1234.56)
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
