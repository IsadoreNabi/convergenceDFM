#' Select optimal number of PLS components with cross-validation
#'
#' Determines the optimal number of Partial Least Squares components using
#' leave-one-out or k-fold cross-validation, with robust error handling.
#'
#'   default) or "CV" (k-fold cross-validation).
#'   validation = "LOO".
#'
#' @return List with components:
#'   \describe{
#'     \item{\code{ncomp_optimal}}{Optimal number of components selected.}
#'     \item{\code{RMSEP}}{Root Mean Squared Error of Prediction for each component.}
#'     \item{\code{R2}}{R-squared values for each component.}
#'     \item{\code{validation_method}}{Validation method used.}
#'     \item{\code{pls_fit}}{Fitted PLS model object with optimal components.}
#'   }
#'
#' @details The function fits PLS models with 1 to \code{max_comp} components,
#'   evaluates each via cross-validation, and selects the number that minimizes
#'   prediction error while avoiding overfitting.
#'
#'
#' @param X First data matrix
#' @param Y Second data matrix
#' @param max_comp Maximum number of components to test (default: 15)
#' @param verbose Logical; print progress and diagnostic information. Default \code{TRUE}.
#' @export
select_optimal_components_safe <- function(X, Y, max_comp = 15, verbose = TRUE) {
  if (verbose) {
    message("========================================")
    message("OPTIMAL COMPONENT SELECTION")
    message("========================================\n")
  }
  
  if (is.null(dim(Y))) Y <- as.matrix(Y)
  
  if (verbose) {
    message("Dimensions for PLS: X =", paste(dim(X), collapse = "x"),
            ", Y =", paste(dim(Y), collapse = "x"))
  }
  
  max_comp <- max(1L, min(max_comp, nrow(X) - 1L, ncol(X)))
  if (verbose) {
    message("Maximum fitted components:", max_comp, "\n")
  }
  
  if (verbose) {
    message("Running PLS with cross-validation...")
  }
  pls_cv <- pls::plsr(Y ~ X, ncomp = max_comp, validation = "CV",
                      segments = min(10L, nrow(X)))
  
  r2_arr <- pls::R2(pls_cv, estimate = "CV")$val
  r2_cv_vec <- as.numeric(r2_arr[1, 1, -1, drop = TRUE])
  
  rmsep_obj <- pls::RMSEP(pls_cv, estimate = "CV")
  rm_cv_vec <- as.numeric(rmsep_obj$val[1, 1, -1, drop = TRUE])
  se_cv_vec <- if (!is.null(rmsep_obj$se)) {
    as.numeric(rmsep_obj$se[1, 1, -1, drop = TRUE])
  } else {
    rep(0, length(rm_cv_vec))
  }
  
  best_rm <- which.min(rm_cv_vec)
  one_se_thr <- rm_cv_vec[best_rm] + se_cv_vec[best_rm]
  one_se_opt <- which(rm_cv_vec <= one_se_thr)[1]
  
  improvements <- c(r2_cv_vec[1], diff(r2_cv_vec))
  elbow <- which(improvements < 0.01)[1]
  if (is.na(elbow)) elbow <- which.max(r2_cv_vec)
  
  optimal_ncomp <- one_se_opt
  if (!is.finite(optimal_ncomp)) optimal_ncomp <- elbow
  if (!is.finite(optimal_ncomp)) optimal_ncomp <- best_rm
  optimal_ncomp <- max(1L, min(optimal_ncomp, max_comp))
  
  varY_explained_pct <- 100 * r2_cv_vec[optimal_ncomp]
  
  if (verbose) {
    message("\nRESULTS:")
    message("- Selected components (one-SE):", optimal_ncomp)
    message("- R2 (CV):", round(r2_cv_vec[optimal_ncomp], 4))
    message("- Variance Y explained (CV):", round(varY_explained_pct, 1), "%\n")
  }
  
  list(
    optimal_ncomp = optimal_ncomp,
    pls_model = pls_cv,
    R2_cv = r2_cv_vec[optimal_ncomp],
    var_explained = varY_explained_pct,
    rmsep_cv = rm_cv_vec
  )
}