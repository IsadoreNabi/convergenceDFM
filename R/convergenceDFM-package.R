#' convergenceDFM: Dynamic Factor Models for Economic Convergence Analysis
#'
#' Tools to analyze economic convergence using Dynamic Factor Models (DFM) and
#' Factor Ornstein-Uhlenbeck (OU) processes. Includes data preparation, PLS-based
#' factor selection, DFM estimation with robust inference, OU estimation (Stan
#' or fallback), formal convergence tests, robustness checks, rotation/coupling
#' tests, and visualization utilities.
#'
#' @section Main functions:
#' \itemize{
#'   \item \code{\link{run_complete_factor_analysis_robust}}: End-to-end pipeline.
#'   \item \code{\link{estimate_DFM}}: Estimate a VAR-based Dynamic Factor Model.
#'   \item \code{\link{estimate_factor_OU}}: Estimate Factor OU model (Stan or fallback).
#'   \item \code{\link{run_convergence_robustness_tests}}: Robustness test suite.
#'   \item \code{\link{visualize_factor_dynamics}}: Visualization helpers.
#' }
#'
#' @importFrom grDevices dev.new dev.cur dev.off pdf heat.colors
#' @importFrom graphics abline barplot box grid layout legend matplot mtext par plot.new text title
#' @importFrom stats approx ar cancor coef complete.cases cor df.residual fitted median model.matrix na.omit p.adjust pnorm predict pt quantile residuals rnorm runif sd var
#' @importFrom parallel detectCores
#' @importFrom utils install.packages setTxtProgressBar txtProgressBar
#' @importFrom magrittr %>%
#' @importFrom dplyr filter arrange group_by summarise transmute
#'
#' @seealso \code{vignette("convergence-analysis")}
#' @keywords internal
"_PACKAGE"

# NOTE: the mutable package environment `.pkgenv` is created and initialised in
# zzz.R (.onLoad). It must remain a real environment; do not re-bind it here.
