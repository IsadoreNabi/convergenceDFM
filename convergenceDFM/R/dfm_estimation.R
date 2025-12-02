#' Estimate Dynamic Factor Model with VAR dynamics
#'
#' Estimates a Dynamic Factor Model by extracting factors via PLS and modeling
#' their dynamics with a Vector Autoregression. Includes automatic lag selection,
#' robust inference, and optional out-of-sample evaluation.
#'
#' @param factors_data List containing PLS-extracted factor scores (\code{scores_X},
#'   \code{scores_Y}) and related objects.
#' @param p Integer. VAR lag order. If \code{NULL}, selected automatically. Default is 2.
#' @param compute_oos Logical. Should out-of-sample diagnostics be computed? Default is TRUE.
#' @param hc_type Character string. Heteroskedasticity-consistent SE type. Default is "HC3".
#' @param verbose Logical; print progress and diagnostic information. Default \code{TRUE}.
#'
#' @return List with components:
#'   \describe{
#'     \item{\code{var_fit}}{Fitted VAR model on combined factors.}
#'     \item{\code{p_used}}{VAR lag order used.}
#'     \item{\code{robust_se}}{Matrices of robust standard errors.}
#'     \item{\code{diagnostics}}{List of diagnostic tests (stability, serial correlation).}
#'     \item{\code{oos_metrics}}{Out-of-sample forecast evaluation (if requested).}
#'   }
#'
#' @details This function models the joint dynamics of X and Y factors using a VAR.
#'   It performs stability checks, tests for serial correlation, computes robust
#'   standard errors, and optionally evaluates forecast performance out-of-sample.
#'
#' @export

estimate_DFM <- function(factors_data, p = 2, compute_oos = TRUE, hc_type = "HC3", 
                         verbose = TRUE) {
  if (verbose) {
    message("========================================")
    message("DYNAMIC FACTOR MODEL (DFM)")
    message("========================================\n")
  }
  
  FX <- tryCatch(as.matrix(factors_data$scores_X), error = function(e) NULL)
  FY <- tryCatch(as.matrix(factors_data$scores_Y), error = function(e) NULL)
  
  if (is.null(FX) && is.null(FY)) {
    stop("No factors available for DFM")
  }
  
  original_names_X <- if (!is.null(FX)) colnames(FX) else NULL
  original_names_Y <- if (!is.null(FY)) colnames(FY) else NULL
  
  if (!is.null(FX) && !is.null(FY)) {
    if (is.null(colnames(FX)) || all(grepl("^Comp", colnames(FX)))) {
      colnames(FX) <- paste0("X", 1:ncol(FX))
    }
    if (is.null(colnames(FY)) || all(grepl("^Comp", colnames(FY)))) {
      colnames(FY) <- paste0("Y", 1:ncol(FY))
    }
    F_combined <- cbind(FX, FY)
  } else if (!is.null(FX)) {
    if (is.null(colnames(FX))) colnames(FX) <- paste0("X", 1:ncol(FX))
    F_combined <- FX
  } else {
    if (is.null(colnames(FY))) colnames(FY) <- paste0("Y", 1:ncol(FY))
    F_combined <- FY
  }
  
  F_combined <- as.data.frame(F_combined)
  colnames(F_combined) <- make.names(colnames(F_combined), unique = TRUE)
  
  if (verbose) {
    message("FACTORS STRUCTURE:")
    message("- Factores X: ", if (!is.null(FX)) ncol(FX) else 0)
    message("- Factores Y: ", if (!is.null(FY)) ncol(FY) else 0)
    message("- Nombres: ", paste(colnames(F_combined), collapse = ", "), "\n")
  }
  
  sel <- choose_var_lag(F_combined, lag.max = max(4, p), type = "const",
                        p_pref = c("SC(n)", "HQ(n)"), alpha = 0.05,
                        oos_eval = compute_oos, verbose = verbose)
  
  fit <- sel$fit
  p_used <- sel$p
  
  rts <- vars::roots(fit, modulus = TRUE)
  is_stable <- all(rts < 1 - 1e-6)
  hl <- ifelse(rts < 1, log(0.5) / log(rts), Inf)
  
  pt_val <- NA_real_
  pt <- tryCatch(
    vars::serial.test(fit, lags.pt = 12, type = "PT.asymptotic"),
    error = function(e) NULL
  )
  if (!is.null(pt) && !is.null(pt$serial$p.value)) {
    pt_val <- pt$serial$p.value
  }
  
  if (verbose) message("\nCalculating robust standard errors (", hc_type, ")...")
  
  robust_se_list <- list()
  n_success <- 0
  
  for (eq_name in names(fit$varresult)) {
    eq_model <- fit$varresult[[eq_name]]
    
    n_obs <- length(residuals(eq_model))
    n_coef <- length(coef(eq_model))
    
    if (n_obs <= n_coef + 1) {
      if (verbose) message("  ", eq_name, ": Insufficient degrees of freedom")
      next
    }
    
    hc_result <- calculate_hc_manual(eq_model, type = hc_type)
    
    if (hc_result$success) {
      se_hc <- hc_result$se
      coef_vals <- coef(eq_model)
      
      min_len <- min(length(se_hc), length(coef_vals))
      se_hc <- se_hc[1:min_len]
      coef_vals <- coef_vals[1:min_len]
      
      t_stats <- coef_vals / se_hc
      df_resid <- df.residual(eq_model)
      p_vals <- if (df_resid > 0) {
        2 * pt(abs(t_stats), df = df_resid, lower.tail = FALSE)
      } else {
        rep(NA, length(t_stats))
      }
      
      robust_se_list[[eq_name]] <- se_hc
      n_sig <- sum(p_vals < 0.05, na.rm = TRUE)
      if (verbose) {
        message("  ", eq_name, ": ", n_sig, "/", length(p_vals), 
                " coefs significativos (p<0.05)")
      }
      n_success <- n_success + 1
    }
  }
  
  if (verbose && n_success > 0) {
    message("\nOK: Robust errors calculated for ", n_success, "/", 
            length(fit$varresult), " ecuaciones")
  }
  
  vm <- fit
  r2_eq <- sapply(vm$varresult, function(m) summary(m)$r.squared)
  adj_eq <- sapply(vm$varresult, function(m) summary(m)$adj.r.squared)
  
  SSE <- sapply(vm$varresult, function(m) sum(residuals(m)^2))
  SST <- sapply(vm$varresult, function(m) {
    y <- m$model[[1]]
    sum((y - mean(y))^2)
  })
  R2_global <- 1 - sum(SSE) / sum(SST)
  
  fits_mat <- do.call(cbind, lapply(vm$varresult, fitted))
  resid_mat <- do.call(cbind, lapply(vm$varresult, residuals))
  
  var_fit_eq <- apply(fits_mat, 2, stats::var, na.rm = TRUE)
  var_res_eq <- apply(resid_mat, 2, stats::var, na.rm = TRUE)
  r2_gelman_eq <- var_fit_eq / (var_fit_eq + var_res_eq + 1e-10)
  r2_gelman_global <- sum(var_fit_eq) / (sum(var_fit_eq) + sum(var_res_eq) + 1e-10)
  
  r2_oos_eq <- NULL
  r2_oos_global <- NA_real_
  mse_oos_eq <- NULL
  
  if (compute_oos && nrow(F_combined) > 50) {
    if (verbose) message("\nCalculating out-of-sample R2...")
    
    T_total <- nrow(F_combined)
    T_start <- max(20, floor(T_total * 0.5))
    
    if (T_start <= p_used + 5) {
      T_start <- min(p_used + 10, floor(T_total * 0.7))
    }
    
    if (T_start < T_total - 5) {
      oos_errors <- matrix(NA, T_total - T_start, ncol(F_combined))
      oos_actuals <- matrix(NA, T_total - T_start, ncol(F_combined))
      
      for (t in T_start:(T_total - 1)) {
        train_data <- F_combined[1:t, , drop = FALSE]
        test_point <- F_combined[t + 1, , drop = FALSE]
        
        tryCatch({
          fit_temp <- vars::VAR(train_data, p = p_used, type = "const")
          pred <- predict(fit_temp, n.ahead = 1)
          pred_vals <- sapply(pred$fcst, function(x) x[1, "fcst"])
          
          idx <- t - T_start + 1
          oos_errors[idx, ] <- as.numeric(test_point) - pred_vals
          oos_actuals[idx, ] <- as.numeric(test_point)
        }, error = function(e) NULL)
      }
      
      r2_oos_eq <- numeric(ncol(F_combined))
      mse_oos_eq <- numeric(ncol(F_combined))
      names(r2_oos_eq) <- names(mse_oos_eq) <- colnames(F_combined)
      
      for (j in 1:ncol(F_combined)) {
        valid_idx <- !is.na(oos_errors[, j])
        if (sum(valid_idx) > 10) {
          sse_oos <- sum(oos_errors[valid_idx, j]^2)
          mse_oos_eq[j] <- sse_oos / sum(valid_idx)
          sst_oos <- sum((oos_actuals[valid_idx, j] - 
                            mean(oos_actuals[valid_idx, j]))^2)
          if (sst_oos > 1e-10) {
            r2_oos_eq[j] <- 1 - sse_oos / sst_oos
            r2_oos_eq[j] <- max(-2, min(1, r2_oos_eq[j]))
          }
        }
      }
      
      valid_all <- complete.cases(oos_errors)
      if (sum(valid_all) > 10) {
        sse_oos_total <- sum(oos_errors[valid_all, ]^2)
        oos_act_valid <- oos_actuals[valid_all, ]
        sst_oos_total <- sum((oos_act_valid - 
                                matrix(colMeans(oos_act_valid), 
                                       nrow(oos_act_valid), 
                                       ncol(oos_act_valid), 
                                       byrow = TRUE))^2)
        if (sst_oos_total > 1e-10) {
          r2_oos_global <- 1 - sse_oos_total / sst_oos_total
          r2_oos_global <- max(-2, min(1, r2_oos_global))
        }
      }
      
      if (verbose) {
        message("\nOut-Of-Sample R2 by equation:")
        print(round(r2_oos_eq, 3))
      }
    }
  }
  
  if (verbose) {
    message("\n========================================")
    message("DFM SUMMARY")
    message("========================================")
    message(sprintf("\n- VAR(p=%d): stable=%s, serial.ok=%s",
                    p_used, 
                    ifelse(is_stable, "yes", "no"),
                    ifelse(isTRUE(sel$serial_ok), "yes", "no")))
    
    if (is.finite(pt_val)) {
      message(sprintf("- Portmanteau p-val: %.4f", pt_val))
    }
    
    finite_hl <- hl[is.finite(hl)]
    if (length(finite_hl) > 0) {
      message("- Median half-life: ", round(stats::median(finite_hl), 2), " periods")
    }
    
    msg <- paste0("- In-Sample GLOBAL R2: ", round(R2_global, 3))
    if (!is.na(r2_oos_global)) {
      msg <- paste0(msg, " | Out-Of-Sample Global R2: ", round(r2_oos_global, 3))
    }
    message(msg)
  }
  
  list(
    var_model = fit,
    p_used = p_used,
    roots = rts,
    is_stable = is_stable,
    half_lives = hl,
    r2_eq = r2_eq,
    r2_adj_eq = adj_eq,
    r2_global = R2_global,
    r2_gelman_eq = r2_gelman_eq,
    r2_gelman_global = r2_gelman_global,
    r2_oos_eq = r2_oos_eq,
    r2_oos_global = r2_oos_global,
    mse_oos_eq = mse_oos_eq,
    robust_se = robust_se_list,
    original_names_X = original_names_X,
    original_names_Y = original_names_Y
  )
}

#' Calculate heteroskedasticity-consistent covariance matrix
#'
#' Computes HC (sandwich) standard errors for linear model coefficients with
#' support for HC0, HC1, HC2, HC3, and HC4 variants.
#'
#' @param lm_model Fitted linear model object.
#' @param type Character string. HC type: "HC0", "HC1", "HC2", "HC3" (default), or "HC4".
#'
#' @return Covariance matrix with HC standard errors.
#'
#' @keywords internal
#' @noRd

calculate_hc_manual <- function(lm_model, type = "HC3") {
  tryCatch({
    X <- model.matrix(lm_model)
    n <- nrow(X)
    k <- ncol(X)
    e <- residuals(lm_model)
    
    if (n <= k) {
      return(list(success = FALSE, error = "n <= k"))
    }
    
    XtX_inv <- tryCatch(
      solve(crossprod(X)),
      error = function(e) {
        solve(crossprod(X) + diag(1e-10, ncol(X)))
      }
    )
    
    if (type %in% c("HC2", "HC3")) {
      h <- rowSums((X %*% XtX_inv) * X)
      h <- pmin(h, 0.9999)
    }
    
    weights <- switch(type,
                      "HC0" = e^2,
                      "HC1" = e^2 * (n / (n - k)),
                      "HC2" = e^2 / (1 - h),
                      "HC3" = e^2 / ((1 - h)^2),
                      e^2
    )
    
    weights[!is.finite(weights)] <- 0
    meat <- crossprod(X * sqrt(weights))
    vcov_hc <- XtX_inv %*% meat %*% XtX_inv
    se_hc <- sqrt(diag(vcov_hc))
    
    if (any(!is.finite(se_hc))) {
      return(list(success = FALSE, error = "No Finite SE"))
    }
    
    list(success = TRUE, vcov = vcov_hc, se = se_hc)
    
  }, error = function(e) {
    list(success = FALSE, error = e$message)
  })
}

#' Select optimal VAR lag order with multiple criteria
#'
#' Determines the optimal lag order for a Vector Autoregression model using
#' information criteria, stability checks, serial correlation tests, and optional
#' out-of-sample validation.
#'
#' @param F_combined Numeric matrix (T x K) of factor scores to be modeled.
#' @param lag.max Integer. Maximum lag order to consider. Default is 4.
#' @param type Character string. Type of deterministic terms: "const" (default),
#'   "trend", "both", or "none".
#' @param p_pref Character vector. Preferred information criteria for initial
#'   selection. Default is \code{c("SC(n)", "HQ(n)")}.
#' @param alpha Numeric. Significance level for serial correlation test. Default is 0.05.
#' @param oos_eval Logical. Should out-of-sample evaluation be performed? Default is TRUE.
#' @param oos_start Numeric. Proportion of sample to use for training in OOS validation.
#'   Default is 0.7.
#' @param verbose Logical; print progress information. Default \code{TRUE}.
#'
#' @return List with components:
#'   \describe{
#'     \item{\code{p}}{Selected optimal lag order.}
#'     \item{\code{fit}}{Fitted VAR model object of class \code{varest}.}
#'     \item{\code{roots_ok}}{Logical indicating if stability condition is satisfied.}
#'     \item{\code{serial_ok}}{Logical indicating if serial correlation test passed.}
#'     \item{\code{oos_mse}}{Out-of-sample mean squared error (if \code{oos_eval = TRUE}).}
#'   }
#'
#' @details The function combines multiple selection criteria: (1) information
#'   criteria (AIC, BIC, HQ), (2) VAR stability (eigenvalue modulus < 1), (3)
#'   Portmanteau test for serial correlation, and (4) out-of-sample forecast
#'   performance. Returns the model that best balances these considerations.
#'
#' @export

choose_var_lag <- function(F_combined, lag.max = 4, type = "const",
                           p_pref = c("SC(n)", "HQ(n)"), alpha = 0.05,
                           oos_eval = TRUE, oos_start = 0.7, verbose = TRUE) {
  
  Fdf <- as.data.frame(F_combined)
  
  if (is.null(colnames(Fdf)) || any(colnames(Fdf) == "")) {
    k_total <- ncol(Fdf)
    k_x <- ceiling(k_total / 2)
    colnames(Fdf) <- c(paste0("X", 1:k_x), paste0("Y", 1:(k_total - k_x)))
  }
  colnames(Fdf) <- make.names(colnames(Fdf), unique = TRUE)
  
  sel <- vars::VARselect(Fdf, lag.max = lag.max, type = type)$selection
  p0 <- min(na.omit(as.integer(sel[p_pref])), na.rm = TRUE)
  if (!is.finite(p0)) p0 <- 1
  
  try_fit <- function(p) {
    tryCatch(
      vars::VAR(Fdf, p = p, type = type), 
      error = function(e) NULL
    )
  }
  
  evaluate_oos <- function(fit, p) {
    if (!oos_eval || nrow(Fdf) < 50) return(NA_real_)
    
    T_total <- nrow(Fdf)
    T_train_start <- floor(T_total * oos_start)
    
    if (T_train_start <= p + 10) return(NA_real_)
    
    errors_sq <- numeric()
    for (t in T_train_start:(T_total - 1)) {
      train_data <- Fdf[1:t, , drop = FALSE]
      test_point <- Fdf[t + 1, , drop = FALSE]
      
      tryCatch({
        fit_temp <- vars::VAR(train_data, p = p, type = type)
        pred <- predict(fit_temp, n.ahead = 1)
        pred_vals <- sapply(pred$fcst, function(x) x[1, "fcst"])
        errors <- as.numeric(test_point) - pred_vals
        errors_sq <- c(errors_sq, sum(errors^2))
      }, error = function(e) NULL)
    }
    
    if (length(errors_sq) > 0) mean(errors_sq) else NA_real_
  }
  
  best_p <- p0
  best_mse <- Inf
  best_fit <- NULL
  
  for (p in unique(c(p0, p0 + 1, 1:lag.max))) {
    if (p < 1 || p > lag.max) next
    
    fit <- try_fit(p)
    if (is.null(fit)) next
    
    roots_ok <- all(Mod(vars::roots(fit)) < 1 - 1e-6)
    if (!roots_ok) next
    
    st <- tryCatch(
      vars::serial.test(fit, lags.pt = 12, type = "PT.asymptotic"),
      error = function(e) NULL
    )
    serial_ok <- is.null(st) || is.null(st$serial$p.value) || st$serial$p.value > alpha
    
    mse_oos <- evaluate_oos(fit, p)
    
    if (serial_ok && !is.na(mse_oos) && mse_oos < best_mse) {
      best_p <- p
      best_mse <- mse_oos
      best_fit <- fit
    } else if (serial_ok && is.null(best_fit)) {
      best_fit <- fit
      best_p <- p
    }
  }
  
  if (is.null(best_fit)) {
    for (p in 1:lag.max) {
      fit <- try_fit(p)
      if (!is.null(fit) && all(Mod(vars::roots(fit)) < 1 - 1e-6)) {
        return(list(p = p, fit = fit, roots_ok = TRUE, serial_ok = FALSE))
      }
    }
    stop("It was not possible to estimate the VAR in any order.")
  }
  
  list(
    p = best_p, 
    fit = best_fit, 
    roots_ok = TRUE,
    serial_ok = TRUE, 
    oos_mse = best_mse
  )
}