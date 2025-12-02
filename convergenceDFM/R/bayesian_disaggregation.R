#' Read CPI data from Excel file
#'
#' Reads Consumer Price Index data from an Excel file with robust error handling
#' and data validation.
#'
#'
#' @return Data frame with CPI data. Returns \code{NULL} if file not found or
#'   read operation fails.
#'
#'
#' @param path_cpi Path to the CPI data file
#' @export

read_cpi <- function(path_cpi) {
  log_msg("INFO", "Reading CPI:", path_cpi)
  
  if (!file.exists(path_cpi)) {
    stop("CPI file not found: ", path_cpi)
  }
  
  df <- readxl::read_excel(path_cpi)
  cn <- tolower(names(df))
  
  col_date <- which.max(stringr::str_detect(cn, "date|year|anio|period"))
  col_cpi <- which.max(stringr::str_detect(cn, "cpi|indice|price|ipc"))
  
  if (length(col_date) == 0 || length(col_cpi) == 0) {
    stop("No date or CPI columns found")
  }
  
  df <- df %>%
    transmute(
      Year = as.integer(substr(as.character(.[[col_date]]), 1, 4)),
      CPI = to_num_commas(.[[col_cpi]])
    ) %>% 
    filter(!is.na(Year), !is.na(CPI)) %>% 
    arrange(Year)
  
  if (anyDuplicated(df$Year)) {
    log_msg("WARN", "CPI with duplicate years; will be added by average.")
    df <- df %>% 
      group_by(Year) %>% 
      summarise(CPI = mean(CPI), .groups = "drop")
  }
  
  df
}

#' Read and normalize industry weights
#'
#' Loads an industry-by-year weight table from Excel, reshapes to long/wide, and
#' normalizes weights to sum to one by year.
#'
#' @param path_weights Path to an Excel file with weights.
#' @return A list with matrix \code{P} (T x K), \code{industries}, and \code{years}.
#' @keywords internal
#' @noRd

read_weights_matrix <- function(path_weights) {
  log_msg("INFO", "Reading weight matrix:", path_weights)
  
  if (!file.exists(path_weights)) {
    stop("Weight file not found: ", path_weights)
  }
  
  raw <- readxl::read_excel(path_weights)
  names(raw)[1] <- "Industry"
  
  long <- raw %>%
    tidyr::pivot_longer(
      cols = -Industry,
      names_to = "Year",
      values_to = "Weight_raw"
    ) %>%
    dplyr::mutate(
      Year = as.integer(stringr::str_extract(Year, "\\d{4}")),
      Weight = to_num_commas(Weight_raw)
    ) %>%
    dplyr::filter(!is.na(Year), !is.na(Weight))
  
  long <- long %>% 
    dplyr::group_by(Year) %>% 
    dplyr::mutate(Weight = Weight / sum(Weight, na.rm = TRUE)) %>% 
    dplyr::ungroup()
  
  wide <- long %>%
    dplyr::select(Industry, Year, Weight) %>%
    tidyr::pivot_wider(
      names_from = Industry,
      values_from = Weight,
      values_fill = 0
    ) %>%
    dplyr::arrange(Year)
  
  P <- as.matrix(wide[, -1, drop = FALSE])
  rownames(P) <- wide$Year
  P <- row_norm1(P)
  
  list(
    P = P,
    industries = colnames(P),
    years = as.integer(rownames(P))
  )
}

#' Likelihood weights from a weight matrix via SVD
#'
#' Uses the first right singular vector (absolute value) as cross-sectional
#' weights and normalizes to sum to one.
#'
#' @param P Numeric matrix (T x K).
#' @return Numeric vector of length K.
#' @keywords internal
#' @noRd

compute_L_from_P <- function(P) {
  X <- scale(P, center = TRUE, scale = FALSE)
  sv <- svd(X)
  l <- abs(sv$v[, 1])
  l / sum(l)
}

#' Likelihood weights from a weight matrix via SVD
#'
#' Uses the first right singular vector (absolute value) as cross-sectional
#' weights and normalizes to sum to one.
#'
#' @param P Numeric matrix (T x K).
#' @return Numeric vector of length K.
#' @keywords internal
#' @noRd

spread_likelihood <- function(L, T_periods, pattern = c("recent", "constant", "linear", "bell")) {
  pattern <- match.arg(pattern)
  L <- L / sum(L)
  t <- seq_len(T_periods)
  
  w <- switch(pattern,
              recent = 0.5 + 0.5 * (t / T_periods),
              constant = rep(1, T_periods),
              linear = 0.3 + 1.4 * ((t - 1) / max(1, T_periods - 1)),
              bell = exp(-(t - T_periods/2)^2 / (2 * (T_periods/4)^2))
  )
  
  LT <- matrix(rep(L, each = T_periods), nrow = T_periods) * w
  row_norm1(LT)
}


#' Weighted posterior update of weights
#'
#' Blends prior \code{P} with time-distributed likelihood \code{LT} via convex
#' combination.
#'
#' @param P Prior T x K matrix.
#' @param LT T x K likelihood matrix.
#' @param lambda Mixing parameter between 0 and 1.
#' @return Updated T x K matrix with row sums equal to 1.
#' @keywords internal
#' @noRd

posterior_weighted <- function(P, LT, lambda = 0.7) {
  row_norm1(lambda * P + (1 - lambda) * LT)
}

#' Disaggregate CPI using custom prior and Bayesian updating
#'
#' Reads CPI and industry weights, optionally replaces weights with a custom prior,
#' and derives a time-consistent disaggregation matrix.
#'
#' @param path_cpi Path to CPI Excel.
#' @param path_weights Path to weights Excel.
#' @param custom_prior Optional length-K numeric prior (sums to 1).
#' @param method Character; reserved for future methods (\code{"weighted"}).
#' @param lambda Mixing weight for posterior update.
#' @return A list with \code{Y_matrix}, \code{industries}, \code{years}.
#' @examples
#' # out <- run_disaggregation_custom_prior("cpi.xlsx","weights.xlsx")
#' @keywords internal
#' @noRd

# Función corregida con el parámetro verbose agregado
# Reemplazá tu función actual con esta versión

run_disaggregation_custom_prior <- function(path_cpi, path_weights,
                                            custom_prior = NULL,
                                            method = "weighted",
                                            lambda = 0.7,
                                            verbose = FALSE) {  # PARÁMETRO AGREGADO
  
  cpi <- read_cpi(path_cpi)
  Wraw <- read_weights_matrix(path_weights)
  Pfull <- Wraw$P
  
  if (!is.null(custom_prior)) {
    for (i in 1:nrow(Pfull)) {
      Pfull[i, ] <- custom_prior
    }
  }
  
  years_common <- intersect(as.integer(rownames(Pfull)), cpi$Year)
  if (length(years_common) < 2) {
    stop("Very few years in common between CPI and WEIGHTS.")
  }
  
  P <- Pfull[as.character(sort(years_common)), , drop = FALSE]
  Tn <- nrow(P)
  
  L <- compute_L_from_P(P)
  LT <- spread_likelihood(L, Tn, pattern = "recent")
  
  W <- posterior_weighted(P, LT, lambda)
  
  list(
    Y_matrix = W,
    industries = colnames(P),
    years = as.integer(rownames(P))
  )
}