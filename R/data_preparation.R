#' Diagnose and prepare data matrices
#'
#' Performs data validation, missing value imputation, and variance checks on
#' input matrices to prepare them for factor analysis. Handles dimension
#' compatibility, NA values, and zero-variance columns.
#'
#' @param X_matrix Numeric matrix or data frame of X variables (e.g., Marxist prices).
#' @param Y_matrix Numeric matrix or data frame of Y variables (e.g., market prices/CPI).
#' @param verbose Logical; print diagnostic information. Default \code{TRUE}.
#'
#' @return List with components:
#'   \describe{
#'     \item{\code{X_matrix}}{Cleaned and prepared X matrix.}
#'     \item{\code{Y_matrix}}{Cleaned and prepared Y matrix.}
#'   }
#'
#' @details The function:
#'   \itemize{
#'     \item Converts to matrix format if needed
#'     \item Validates dimensional compatibility
#'     \item Imputes missing values via interpolation (using zoo if available)
#'     \item Adds minimal noise to zero-variance columns
#'     \item Reports diagnostic information
#'   }
#'
#' @export

diagnose_data <- function(X_matrix, Y_matrix, verbose = TRUE) {
  if (verbose) {
    message("========================================")
    message("DATA DIAGNOSIS")
    message("========================================\n")
  }
  
  na_approx_vec <- function(v) {
    if (all(!is.finite(v))) return(rep(NA_real_, length(v)))
    
    if (requireNamespace("zoo", quietly = TRUE)) {
      out <- zoo::na.approx(v, na.rm = FALSE)
      if (anyNA(out)) out[is.na(out)] <- mean(out, na.rm = TRUE)
      return(out)
    } else {
      idx <- which(is.finite(v))
      if (length(idx) == 0) return(rep(NA_real_, length(v)))
      approx(x = idx, y = v[idx], xout = seq_along(v), rule = 2)$y
    }
  }
  
  if (verbose) {
    message("DATA STRUCTURE:")
    message("-------------------")
    message("X_matrix:")
    message("  - Class: ", paste(class(X_matrix), collapse = " "))
    message("  - Dimensions: ", paste(dim(X_matrix), collapse = " "))
    
    message("\nY_matrix:")
    message("  - Class: ", paste(class(Y_matrix), collapse = " "))
    message("  - Dimensions: ", paste(dim(Y_matrix), collapse = " "))
  }
  
  if (!is.matrix(X_matrix)) {
    if (verbose) message("\nWARNING: Converting X_matrix to matrix...")
    X_matrix <- as.matrix(X_matrix)
  }
  
  if (!is.matrix(Y_matrix)) {
    if (verbose) message("\nWARNING: Converting Y_matrix to matrix...")
    Y_matrix <- as.matrix(Y_matrix)
  }
  
  if (nrow(X_matrix) != nrow(Y_matrix)) {
    stop("ERROR: X_matrix and Y_matrix have different number of rows")
  }
  if (verbose) message("\nOK: Compatible dimensions: ", nrow(X_matrix), " observations")
  
  na_x <- sum(!is.finite(X_matrix))
  na_y <- sum(!is.finite(Y_matrix))
  
  if (verbose) {
    message("\nMISSING/INFINITE VALUES:")
    message("  - In X_matrix: ", na_x, " (", round(100*na_x/length(X_matrix), 2), "%)")
    message("  - In Y_matrix: ", na_y, " (", round(100*na_y/length(Y_matrix), 2), "%)")
  }
  
  if (na_x > 0) {
    if (verbose) message("WARNING: Imputing X by interpolation...")
    for (j in seq_len(ncol(X_matrix))) {
      v <- X_matrix[, j]
      v[!is.finite(v)] <- NA_real_
      if (anyNA(v)) v <- na_approx_vec(v)
      if (anyNA(v)) v[is.na(v)] <- mean(v, na.rm = TRUE)
      X_matrix[, j] <- v
    }
  }
  
  if (na_y > 0) {
    if (verbose) message("WARNING: Imputing Y by interpolation...")
    for (j in seq_len(ncol(Y_matrix))) {
      v <- Y_matrix[, j]
      v[!is.finite(v)] <- NA_real_
      if (anyNA(v)) v <- na_approx_vec(v)
      if (anyNA(v)) v[is.na(v)] <- mean(v, na.rm = TRUE)
      Y_matrix[, j] <- v
    }
  }
  
  vx <- apply(X_matrix, 2, stats::var)
  vy <- apply(Y_matrix, 2, stats::var)
  
  zero_var_x <- which(vx < 1e-10)
  zero_var_y <- which(vy < 1e-10)
  
  if (length(zero_var_x) > 0) {
    if (verbose) message("WARNING: Adding noise to columns with variance ~0 in X")
    X_matrix[, zero_var_x] <- X_matrix[, zero_var_x] + 
      matrix(rnorm(nrow(X_matrix) * length(zero_var_x), 0, 1e-6), 
             nrow(X_matrix))
  }
  
  if (length(zero_var_y) > 0) {
    if (verbose) message("WARNING: Adding noise to columns with variance ~0 in Y")
    Y_matrix[, zero_var_y] <- Y_matrix[, zero_var_y] + 
      matrix(rnorm(nrow(Y_matrix) * length(zero_var_y), 0, 1e-6), 
             nrow(Y_matrix))
  }
  
  if (verbose) message("\nOK: Data ready for analysis\n")
  
  list(X_matrix = X_matrix, Y_matrix = Y_matrix)
}