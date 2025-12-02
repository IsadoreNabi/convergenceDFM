#' Log message with timestamp and level
#'
#' Internal logging function that prints messages with timestamp and severity level.
#' Respects the package-level log threshold set in \code{.pkgenv$log_level}.
#'
#' @param level Character string indicating message severity. One of "DEBUG", "INFO",
#'   "WARN", or "ERROR". Defaults to "INFO".
#' @param ... Message components to be concatenated and logged.
#'
#' @return Invisibly returns \code{NULL}. Called for side effect of printing message.
#'
#' @keywords internal
#' @noRd

log_level <- "INFO"  
.level_order <- c("DEBUG"=1, "INFO"=2, "WARN"=3, "ERROR"=4)

log_msg <- function(level="INFO", ...) {
  if (.level_order[[level]] >= .level_order[[log_level]]) {
    message(sprintf("[%s][%s] %s",
                    format(Sys.time(), "%Y-%m-%d %H:%M:%S"), 
                    level,
                    paste(..., collapse=" ")))
  }
}

#' Convert European number format to numeric
#'
#' Converts character strings in European number format (using comma as decimal
#' separator and period as thousands separator) to numeric values. Already numeric
#' inputs are returned unchanged.
#'
#' @param x Numeric value or character string in European format (e.g., "1.234,56").
#'
#' @return Numeric value. Returns \code{NA} if conversion fails.
#'
#' @examples
#' \donttest{
#' to_num_commas("1.234,56")
#' to_num_commas("1234.56")
#' to_num_commas(1234.56)
#' }
#'
#' @export

to_num_commas <- function(x) {
  if (is.numeric(x)) return(as.numeric(x))
  x <- stringr::str_replace_all(as.character(x), "\\.", "")
  x <- stringr::str_replace_all(x, ",", ".")
  suppressWarnings(as.numeric(x))
}

#' Normalize matrix rows to sum to one
#'
#' Scales each row of a matrix so that its elements sum to 1. Handles zero-sum
#' rows by leaving them unchanged (avoiding division by zero).
#'
#' @param M Numeric matrix to be row-normalized.
#'
#' @return Matrix of the same dimensions as \code{M} with each row summing to 1
#'   (except rows that originally summed to zero).
#'
#' @examples
#' \donttest{
#' M <- matrix(c(1, 2, 3, 4, 5, 6), nrow = 2)
#' M_norm <- row_norm1(M)
#' rowSums(M_norm)  # Should be c(1, 1)
#' }
#'
#' @export

row_norm1 <- function(M) {
  rs <- rowSums(M)
  rs[rs == 0] <- 1
  sweep(M, 1, rs, "/")
}

#' Safe division with floor protection
#'
#' Performs division while protecting against division by zero and ensuring both
#' numerator and denominator are above machine epsilon.
#'
#' @param p Numeric vector, numerator.
#' @param q Numeric vector, denominator.
#'
#' @return Numeric vector of length equal to \code{max(length(p), length(q))}.
#'
#' @keywords internal
#' @noRd

safe_div <- function(p, q) { 
  p <- pmax(p, .Machine$double.eps)
  q <- pmax(q, .Machine$double.eps)
  p/q 
}

#' Robust correlation coefficient
#'
#' Computes correlation between two vectors using both Pearson and Spearman methods,
#' returning the one with larger absolute value. Handles missing data and edge cases.
#'
#' @param x Numeric vector.
#' @param y Numeric vector of the same length as \code{x}.
#'
#' @return Numeric scalar. Returns 0 if fewer than 3 valid pairs exist or if both
#'   correlation methods fail.
#'
#' @keywords internal
#' @noRd

robust_cor <- function(x, y) {
  if (length(x) != length(y) || length(x) < 3) return(0)
  
  valid <- complete.cases(x, y)
  if (sum(valid) < 3) return(0)
  
  x <- x[valid]
  y <- y[valid]
  
  cp <- suppressWarnings(cor(x, y, method = "pearson"))
  cs <- suppressWarnings(cor(x, y, method = "spearman"))
  
  if (is.na(cp) && is.na(cs)) return(0)
  if (is.na(cp)) return(cs)
  if (is.na(cs)) return(cp)
  if (abs(cs) > abs(cp)) cs else cp
}

#' Geometric mean
#'
#' Computes the geometric mean of a numeric vector by exponentiating the mean
#' of logarithms.
#'
#' @param x Numeric vector with positive values.
#'
#' @return Numeric scalar, the geometric mean.
#'
#' @keywords internal
#' @noRd

geometric_mean_robust <- function(x, epsilon = 1e-10) {
  x <- pmax(x, epsilon)
  exp(mean(log(x)))
}