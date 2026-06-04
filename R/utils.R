#' Get or set the package logging threshold
#'
#' Internal helpers to read and update the package-level logging threshold stored
#' in the private environment \code{.pkgenv}. Messages with a severity below the
#' threshold are suppressed by \code{log_msg()}.
#'
#' @param level Character string, one of "DEBUG", "INFO", "WARN", "ERROR".
#'
#' @return \code{get_log_level()} returns the current threshold (character);
#'   \code{set_log_level()} returns the previous threshold invisibly.
#'
#' @keywords internal
#' @noRd

.level_order <- c("DEBUG" = 1L, "INFO" = 2L, "WARN" = 3L, "ERROR" = 4L)

get_log_level <- function() {
  lv <- tryCatch(get("log_level", envir = .pkgenv), error = function(e) NULL)
  if (is.null(lv) || !lv %in% names(.level_order)) "INFO" else lv
}

set_log_level <- function(level = c("DEBUG", "INFO", "WARN", "ERROR")) {
  level <- match.arg(level)
  old <- get_log_level()
  assign("log_level", level, envir = .pkgenv)
  invisible(old)
}

#' Log message with timestamp and level
#'
#' Internal logging function that prints messages with timestamp and severity
#' level, respecting the package-level threshold in \code{.pkgenv$log_level}
#' (see \code{set_log_level()}). Messages are emitted through \code{message()} so
#' they go to \code{stderr()} and can be suppressed with
#' \code{suppressMessages()}.
#'
#' @param level Character string indicating message severity. One of "DEBUG",
#'   "INFO", "WARN", or "ERROR". Defaults to "INFO".
#' @param ... Message components to be concatenated and logged.
#'
#' @return Invisibly returns \code{NULL}. Called for side effect of printing.
#'
#' @keywords internal
#' @noRd

log_msg <- function(level = "INFO", ...) {
  if (!level %in% names(.level_order)) level <- "INFO"
  if (.level_order[[level]] >= .level_order[[get_log_level()]]) {
    message(sprintf("[%s][%s] %s",
                    format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
                    level,
                    paste(..., collapse = " ")))
  }
  invisible(NULL)
}

#' Convert localized number strings to numeric
#'
#' Converts character strings that may use either the European convention
#' (comma decimal separator, period thousands separator, e.g. "1.234,56") or the
#' Anglo-Saxon convention (period decimal separator, comma thousands separator,
#' e.g. "1,234.56") to numeric. Plain numbers and already-numeric inputs are
#' returned unchanged. The function is vectorized.
#'
#' @details Disambiguation rules, applied per element after stripping whitespace:
#'   \itemize{
#'     \item If both \code{.} and \code{,} are present, the \emph{rightmost} of
#'       the two is treated as the decimal mark and the other as the thousands
#'       separator.
#'     \item If only commas are present: a single comma is treated as a decimal
#'       mark; multiple commas are treated as thousands separators and removed.
#'     \item If only periods are present: a single period is treated as a decimal
#'       mark (Anglo-Saxon); multiple periods are treated as thousands separators
#'       and removed.
#'   }
#'   The single ambiguous case "1.234" (one period, no comma) is interpreted as
#'   the decimal 1.234, not as 1234.
#'
#' @param x Numeric value or character vector with localized number strings.
#'
#' @return Numeric vector the same length as \code{x}. Entries that cannot be
#'   parsed become \code{NA}.
#'
#' @examples
#' to_num_commas("1.234,56")   # European  -> 1234.56
#' to_num_commas("1,234.56")   # US        -> 1234.56
#' to_num_commas("1234.56")    # plain     -> 1234.56
#' to_num_commas(1234.56)      # numeric   -> 1234.56
#'
#' @export

to_num_commas <- function(x) {
  if (is.numeric(x)) return(as.numeric(x))

  parse_one <- function(s) {
    if (is.na(s)) return(NA_real_)
    s <- gsub("[[:space:]]", "", as.character(s))
    if (!nzchar(s)) return(NA_real_)

    pos_dot   <- gregexpr(".", s, fixed = TRUE)[[1]]
    pos_comma <- gregexpr(",", s, fixed = TRUE)[[1]]
    n_dot   <- if (pos_dot[1]   == -1L) 0L else length(pos_dot)
    n_comma <- if (pos_comma[1] == -1L) 0L else length(pos_comma)

    if (n_dot > 0L && n_comma > 0L) {
      if (max(pos_dot) > max(pos_comma)) {
        s <- gsub(",", "", s, fixed = TRUE)       # comma = thousands
      } else {
        s <- gsub(".", "", s, fixed = TRUE)       # dot   = thousands
        s <- sub(",", ".", s, fixed = TRUE)       # comma = decimal
      }
    } else if (n_comma > 0L) {
      if (n_comma == 1L) {
        s <- sub(",", ".", s, fixed = TRUE)       # single comma = decimal
      } else {
        s <- gsub(",", "", s, fixed = TRUE)       # repeated commas = thousands
      }
    } else if (n_dot > 1L) {
      s <- gsub(".", "", s, fixed = TRUE)         # repeated dots = thousands
    }

    suppressWarnings(as.numeric(s))
  }

  vapply(x, parse_one, numeric(1), USE.NAMES = FALSE)
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
#' M <- matrix(c(1, 2, 3, 4, 5, 6), nrow = 2)
#' M_norm <- row_norm1(M)
#' rowSums(M_norm)  # Should be c(1, 1)
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
#' numerator and denominator are above machine epsilon. Intended for ratios of
#' non-negative quantities (e.g. variance shares); not a general-purpose divide.
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
  p / q
}

#' Geometric mean (numerically guarded)
#'
#' Computes the geometric mean of a numeric vector by exponentiating the mean of
#' logarithms, flooring values at \code{epsilon} to avoid \code{log(0)}.
#'
#' @param x Numeric vector with non-negative values.
#' @param epsilon Positive floor applied before taking logs. Default \code{1e-10}.
#'
#' @return Numeric scalar, the geometric mean. Returns \code{NA_real_} for an
#'   empty input.
#'
#' @keywords internal
#' @noRd

geometric_mean_robust <- function(x, epsilon = 1e-10) {
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  x <- pmax(x, epsilon)
  exp(mean(log(x)))
}
