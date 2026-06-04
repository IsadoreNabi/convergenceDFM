#' Select optimal number of PLS components with cross-validation
#'
#' Determines the optimal number of Partial Least Squares (PLS) components using
#' k-fold or leave-one-out cross-validation, with robust error handling. The
#' number of components is chosen with the one-standard-error rule (the most
#' parsimonious model whose cross-validated RMSEP is within one standard error of
#' the minimum), falling back to an R2 elbow and then to the RMSEP minimizer.
#'
#' @param X First data matrix (predictors), T x p.
#' @param Y Second data matrix (responses), T x q (or a vector).
#' @param max_comp Maximum number of components to test. Default 15. It is
#'   internally capped at \code{min(nrow(X) - 1, ncol(X))}.
#' @param validation Character; cross-validation type, one of \code{"CV"} (k-fold,
#'   default) or \code{"LOO"} (leave-one-out, deterministic).
#' @param segments Integer; number of folds when \code{validation = "CV"}.
#'   Default \code{min(10, nrow(X))}.
#' @param seed Optional integer. With \code{validation = "CV"} the folds are
#'   random; supplying a seed makes the selection reproducible. The RNG state is
#'   restored on exit. Default \code{NULL}.
#' @param verbose Logical; print progress and diagnostic information. Default
#'   \code{TRUE}.
#'
#' @return A list with components:
#'   \describe{
#'     \item{\code{optimal_ncomp}}{Integer, number of components selected.}
#'     \item{\code{pls_model}}{Fitted \code{pls::plsr} model (all \code{max_comp}
#'       components, validated).}
#'     \item{\code{R2_cv}}{Cross-validated R2 at \code{optimal_ncomp}.}
#'     \item{\code{var_explained}}{Percent of Y variance explained (CV) at
#'       \code{optimal_ncomp}.}
#'     \item{\code{rmsep_cv}}{Numeric vector of cross-validated RMSEP for
#'       1..\code{max_comp} components.}
#'     \item{\code{validation}}{The validation method actually used.}
#'   }
#'
#' @details The function fits PLS models with 1 to \code{max_comp} components,
#'   evaluates each via cross-validation, and selects the number that minimizes
#'   prediction error while avoiding overfitting (one-SE rule).
#'
#' @export
select_optimal_components_safe <- function(X, Y, max_comp = 15,
                                           validation = c("CV", "LOO"),
                                           segments = NULL, seed = NULL,
                                           verbose = TRUE) {
  validation <- match.arg(validation)

  if (verbose) {
    message("========================================")
    message("OPTIMAL COMPONENT SELECTION")
    message("========================================\n")
  }

  if (is.null(dim(Y))) Y <- as.matrix(Y)

  if (verbose) {
    message("Dimensions for PLS: X = ", paste(dim(X), collapse = "x"),
            ", Y = ", paste(dim(Y), collapse = "x"))
  }

  max_comp <- max(1L, min(max_comp, nrow(X) - 1L, ncol(X)))
  if (verbose) message("Maximum fitted components: ", max_comp, "\n")

  if (validation == "CV" && !is.null(seed)) {
    if (exists(".Random.seed", envir = .GlobalEnv)) {
      old_seed <- get(".Random.seed", envir = .GlobalEnv)
      on.exit(assign(".Random.seed", old_seed, envir = .GlobalEnv), add = TRUE)
    }
    set.seed(seed)
  }

  if (verbose) message("Running PLS with ", validation, " validation...")

  if (validation == "LOO") {
    pls_cv <- pls::plsr(Y ~ X, ncomp = max_comp, validation = "LOO")
  } else {
    if (is.null(segments)) segments <- min(10L, nrow(X))
    pls_cv <- pls::plsr(Y ~ X, ncomp = max_comp, validation = "CV",
                        segments = segments)
  }

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
    message("- Selected components (one-SE): ", optimal_ncomp)
    message("- R2 (CV): ", round(r2_cv_vec[optimal_ncomp], 4))
    message("- Variance Y explained (CV): ", round(varY_explained_pct, 1), "%\n")
  }

  list(
    optimal_ncomp = optimal_ncomp,
    pls_model = pls_cv,
    R2_cv = r2_cv_vec[optimal_ncomp],
    var_explained = varY_explained_pct,
    rmsep_cv = rm_cv_vec,
    validation = validation
  )
}
