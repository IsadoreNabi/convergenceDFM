# =============================================================================
# Excel input/output for the disaggregation/convergence workflow.
#
# These are generic readers (parse an aggregate CPI series and a sector-by-year
# weight table from Excel). They are intentionally kept in convergenceDFM rather
# than imported from BayesianDisaggregation: the canonical asset reused across
# the two packages is the disaggregation ENGINE (BayesianDisaggregation::
# disaggregate_conjugate), which operates on a numeric CPI vector and a weight
# matrix, not on file paths. The readers here are the hardened versions (clean
# error on an unidentifiable column; see NEWS 0.2.0) and are the single readers
# convergenceDFM relies on.
# =============================================================================

#' Read CPI data from an Excel file
#'
#' Reads an aggregate Consumer Price Index series from an Excel file, with robust
#' column detection and data validation. The date/year and CPI columns are
#' located by pattern-matching on lower-cased header names; localized numerics
#' (comma decimals) are parsed with [to_num_commas()]; duplicate years are
#' collapsed by averaging.
#'
#' @param path_cpi Path to the CPI Excel file.
#'
#' @return A `data.frame` with two columns: `Year` (integer) and `CPI`
#'   (numeric), sorted ascending by `Year`. The function errors (it does not
#'   return `NULL`) when the file is missing or no date/CPI column can be
#'   identified.
#'
#' @seealso [to_num_commas()]
#' @export
read_cpi <- function(path_cpi) {
  log_msg("INFO", "Reading CPI:", path_cpi)

  if (!file.exists(path_cpi)) {
    stop("CPI file not found: ", path_cpi)
  }

  df <- readxl::read_excel(path_cpi)
  cn <- tolower(names(df))

  # which.max() on a logical returns 1 even when nothing matches, which silently
  # selected the first column. Use which() and fail cleanly if no column matches.
  date_hits <- which(stringr::str_detect(cn, "date|year|anio|period"))
  cpi_hits  <- which(stringr::str_detect(cn, "cpi|indice|price|ipc"))

  if (length(date_hits) == 0 || length(cpi_hits) == 0) {
    stop("No date or CPI columns found (looked for date/year/anio/period and ",
         "cpi/indice/price/ipc in the column names).", call. = FALSE)
  }
  col_date <- date_hits[1]
  col_cpi  <- cpi_hits[1]

  df <- df %>%
    transmute(
      Year = as.integer(stringr::str_extract(as.character(.[[col_date]]), "\\d{4}")),
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
#' normalizes weights to sum to one by year. The first column is assumed to
#' contain sector names (renamed to `Industry`); the remaining columns are years.
#'
#' @param path_weights Path to an Excel file with weights.
#' @return A list with matrix `P` (T x K, rows sum to 1), `industries`, and
#'   `years`.
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
