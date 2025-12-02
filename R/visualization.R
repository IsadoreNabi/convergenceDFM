#' Simple factor dynamics visualization
#'
#' Creates a simplified multi-panel plot of factor scores over time with
#' minimal dependencies. Useful for quick diagnostics.
#'
#' @param factors_data List containing factor scores.
#' @param output_file Character string. Optional file path for saving. Default is \code{NULL}.
#' @param verbose Logical; print progress and diagnostic information. Default \code{TRUE}.
#'
#' @return Invisibly returns \code{NULL}.
#'
#' @examples
#' \donttest{
#' data <- list(scores_X = matrix(rnorm(100), ncol = 2),
#'              scores_Y = matrix(rnorm(100), ncol = 2))
#' tmp_file <- file.path(tempdir(), "test_plot.pdf")
#' visualize_factor_dynamics_simple(data, output_file = tmp_file)
#' }
#' @export

visualize_factor_dynamics_simple <- function(factors_data, output_file = NULL, verbose = TRUE) {
  if (verbose) {
    message("\nCreating alternative simple visualization...")
  }
  
  if (!is.null(output_file)) {
    pdf(output_file, width = 10, height = 6)
    on.exit(dev.off())
  }
  
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar), add = TRUE)
  
  par(mfrow = c(1, 2))
  
  if (!is.null(factors_data$scores_X)) {
    X_data <- as.matrix(factors_data$scores_X)
    if (ncol(X_data) > 0 && nrow(X_data) > 1) {
      plot(X_data[, 1], type = "l", 
           main = "First X Factor", 
           ylab = "Score", xlab = "Time",
           col = "blue", lwd = 2)
      grid()
    }
  }
  
  if (!is.null(factors_data$scores_Y)) {
    Y_data <- as.matrix(factors_data$scores_Y)
    if (ncol(Y_data) > 0 && nrow(Y_data) > 1) {
      plot(Y_data[, 1], type = "l", 
           main = "First Y Factor", 
           ylab = "Score", xlab = "Time",
           col = "red", lwd = 2)
      grid()
    }
  }
  
  if (verbose) {
    message("  OK: Simple visualization completed")
  }
  return(invisible(TRUE))
}

#' Visualize factor dynamics comprehensively
#'
#' Creates a multi-panel visualization summarizing all aspects of the factor
#' model: scores over time, loadings, correlations, OU dynamics, and convergence
#' patterns.
#'
#'   Default is \code{NULL} (display only).
#'
#' @return Invisibly returns \code{NULL}. Called for side effect of creating plots.
#'
#'
#' @param dfm_result Result object from DFM analysis
#' @param ou_result Result object from OU estimation
#' @param factors_data Data frame with factor information
#' @param save_plot Logical, save plot to file (default: FALSE)
#' @param plot_file File name for saved plot (default: 'factor_dynamics.pdf')
#' @param use_device Graphics device to use: 'default', 'pdf', 'png' (default: 'default')
#' @param verbose Logical; print progress and diagnostic information. Default \code{TRUE}.
#' @export

visualize_factor_dynamics <- function(dfm_result, ou_result, factors_data, 
                                      save_plot = FALSE, 
                                      plot_file = NULL,
                                      use_device = "default",
                                      verbose = TRUE) {
  
  if (verbose) {
    message("Visualizing factor dynamics...")
  }
  
  open_graphics_device <- function(use_device) {
    device_opened <- FALSE
    
    if (use_device == "pdf" || save_plot) {
      tryCatch({
        pdf(plot_file, width = 10, height = 8)
        device_opened <- TRUE
        if (verbose) {
          message("  Dispositivo PDF abierto:", plot_file)
        }
      }, error = function(e) {
        if (verbose) {
          message("  WARNING: Could not open PDF:", e$message)
        }
      })
    }
    
    if (!device_opened && use_device != "none") {
      in_rstudio <- Sys.getenv("RSTUDIO") == "1"
      
      if (in_rstudio) {
        tryCatch({
          if (dev.cur() == 1) {
            dev.new(noRStudioGD = FALSE)
          }
          device_opened <- TRUE
        }, error = function(e) {
          if (verbose) {
            message("  WARNING: Using existing RStudio device")
          }
          device_opened <- TRUE
        })
      } else {
        tryCatch({
          grDevices::dev.new(width = 10, height = 8)
          device_opened <- TRUE
          if (verbose) {
            message("  OK: Graphics device opened successfully")
          }
        }, error = function(e) {
          if (verbose) {
            message("  WARNING: Could not open graphics device:", e$message)
          }
        })
      }
    }
    
    return(device_opened)
  }
  
  setup_plot_layout <- function() {
    old_par <- NULL
    layout_set <- FALSE
    
    tryCatch({
      old_par <- par(no.readonly = TRUE)
      on.exit(par(old_par))
      
      par(mfrow = c(2, 2))
      par(mar = c(4, 4, 3, 1))
      par(oma = c(0, 0, 2, 0))
      layout_set <- TRUE
      
      if (verbose) {
        message("  OK: Configured layout (simple method)")
      }
    }, error = function(e) {
      if (verbose) {
        message("  WARNING: Error in simple configuration:", e$message)
      }
    })
    
    if (!layout_set) {
      tryCatch({
        layout(matrix(1:4, 2, 2))
        layout_set <- TRUE
        if (verbose) {
          message("  OK: Configured layout (layout method)")
        }
      }, error = function(e) {
        if (verbose) {
          message("  WARNING: Layout error:", e$message)
        }
      })
    }
    
    if (!layout_set) {
      tryCatch({
        old_par <- par(no.readonly = TRUE)
        on.exit(par(old_par))
        par(mfrow = c(2, 2))
        layout_set <- TRUE
        if (verbose) {
          message("  OK: Configured layout (minimal method)")
        }
      }, error = function(e) {
        if (verbose) {
          message("  OK: Layout could not be configured")
        }
      })
    }
    
    return(list(old_par = old_par, success = layout_set))
  }
  
  safe_plot <- function(plot_fn, main_title = "", fallback_msg = "Without data") {
    plot_created <- FALSE
    
    tryCatch({
      plot_fn()
      plot_created <- TRUE
    }, error = function(e) {
      tryCatch({
        plot(1, type = "n", axes = FALSE, xlab = "", ylab = "", 
             xlim = c(0, 1), ylim = c(0, 1))
        text(0.5, 0.5, fallback_msg, cex = 1.2, col = "gray50")
        box()
        title(main = main_title)
        plot_created <- TRUE
      }, error = function(e2) {
        tryCatch({
          plot.new()
          text(0.5, 0.5, fallback_msg)
          title(main = main_title)
          plot_created <- TRUE
        }, error = function(e3) {
          if (verbose) {
            message("    ERROR: Could not create panel:", main_title)
          }
        })
      })
    })
    
    return(plot_created)
  }
  
  extract_valid_data <- function(data, max_cols = 3) {
    if (is.null(data)) return(NULL)
    
    if (!is.matrix(data)) {
      data <- as.matrix(data)
    }
    
    if (nrow(data) < 2 || ncol(data) < 1) return(NULL)
    
    n_cols <- min(max_cols, ncol(data))
    data <- data[, 1:n_cols, drop = FALSE]
    
    if (!any(is.finite(data))) return(NULL)
    
    return(data)
  }
  
  device_opened <- open_graphics_device(use_device)
  
  if (!device_opened && use_device != "none") {
    if (verbose) {
      message("WARNING: Could not open graphics device, trying current device...")
    }
  }
  
  layout_config <- setup_plot_layout()
  
  if (!layout_config$success) {
    if (verbose) {
      message("ERROR: The display layout could not be configured.")
    }
    if (save_plot && dev.cur() > 1) dev.off()
    return(invisible(FALSE))
  }
  
  on.exit({
    if (save_plot && dev.cur() > 1) {
      dev.off()
      if (verbose) {
        message("  OK: Plot saved in:", plot_file)
      }
    }
    
    if (!save_plot && !is.null(layout_config$old_par)) {
      tryCatch({
        suppressWarnings(par(layout_config$old_par))
      }, error = function(e) {})
    }
  })
  
  plots_created <- 0
  
  if (verbose) {
    message("  Creating Panel 1: X Factors...")
  }
  
  X_data <- extract_valid_data(factors_data$scores_X)
  
  plot1_ok <- safe_plot(
    function() {
      if (!is.null(X_data)) {
        cols <- c("blue", "red", "darkgreen")[1:ncol(X_data)]
        
        matplot(X_data, type = "l", lty = 1, lwd = 2, col = cols,
                main = "Dynamics of X Factors", 
                xlab = "Time", ylab = "Score",
                las = 1)
        
        grid(col = "gray90", lty = 2)
        
        if (ncol(X_data) > 1) {
          legend("topright", 
                 legend = paste("X", 1:ncol(X_data)), 
                 col = cols, lty = 1, lwd = 2, 
                 bty = "n", cex = 0.9)
        }
        
        abline(h = 0, col = "gray30", lty = 2)
      } else {
        plot(1, type = "n", axes = FALSE, xlab = "", ylab = "")
        text(1, 1, "Sin datos X", cex = 1.2, col = "gray50")
        box()
        title("Factores X")
      }
    },
    "X Factors",
    "No X scores available"
  )
  if (plot1_ok) plots_created <- plots_created + 1
  
  if (verbose) {
    message("  Creating Panel 2: Y Factors...")
  }
  
  Y_data <- extract_valid_data(factors_data$scores_Y)
  
  plot2_ok <- safe_plot(
    function() {
      if (!is.null(Y_data)) {
        cols <- c("darkblue", "darkred", "darkgreen")[1:ncol(Y_data)]
        
        matplot(Y_data, type = "l", lty = 1, lwd = 2, col = cols,
                main = "Y Factors Dynamics", 
                xlab = "Time", ylab = "Score",
                las = 1)
        
        grid(col = "gray90", lty = 2)
        
        if (ncol(Y_data) > 1) {
          legend("topright", 
                 legend = paste("Y", 1:ncol(Y_data)), 
                 col = cols, lty = 1, lwd = 2, 
                 bty = "n", cex = 0.9)
        }
        
        abline(h = 0, col = "gray30", lty = 2)
      } else {
        plot(1, type = "n", axes = FALSE, xlab = "", ylab = "")
        text(1, 1, "Sin datos Y", cex = 1.2, col = "gray50")
        box()
        title("Factores Y")
      }
    },
    "Y Factors",
    "No scores available"
  )
  if (plot2_ok) plots_created <- plots_created + 1
  
  if (verbose) {
    message("  Creating panel 3: DFM half-lives ...")
  }
  
  plot3_ok <- safe_plot(
    function() {
      hl_dfm <- NULL
      
      if (!is.null(dfm_result) && !is.null(dfm_result$half_lives)) {
        hl_dfm <- dfm_result$half_lives[is.finite(dfm_result$half_lives)]
        hl_dfm <- hl_dfm[hl_dfm > 0 & hl_dfm < 1000]
      }
      
      if (!is.null(hl_dfm) && length(hl_dfm) > 0) {
        bar_names <- if (length(hl_dfm) <= 10) {
          paste0("VAR", 1:length(hl_dfm))
        } else {
          rep("", length(hl_dfm))
        }
        
        cols <- heat.colors(length(hl_dfm))
        
        bp <- barplot(hl_dfm, 
                      names.arg = bar_names,
                      main = "Half-lives (DFM/VAR)", 
                      ylab = "Periods", 
                      col = cols,
                      las = 1,
                      ylim = c(0, max(hl_dfm) * 1.1))
        
        if (length(hl_dfm) <= 10) {
          text(bp, hl_dfm + max(hl_dfm) * 0.02, 
               round(hl_dfm, 1), 
               cex = 0.8, pos = 3)
        }
        
        abline(h = median(hl_dfm), col = "red", lty = 2, lwd = 2)
        
        mtext(paste("Median:", round(median(hl_dfm), 1), "periods"), 
              side = 3, line = 0, cex = 0.8, col = "red")
      } else {
        plot(1, type = "n", axes = FALSE, xlab = "", ylab = "")
        text(1, 1, "Sin half-lives DFM", cex = 1.2, col = "gray50")
        box()
        title("DFM half-lives")
      }
    },
    "Half-lives (DFM/VAR)",
    "No available DFM half-lives"
  )
  if (plot3_ok) plots_created <- plots_created + 1
  
  if (verbose) {
    message("  Creating panel 4: Half-lives OU...")
  }
  
  plot4_ok <- safe_plot(
    function() {
      hl_ou <- c()
      
      if (!is.null(ou_result)) {
        if (!is.null(ou_result$half_lives_X)) {
          hl_x <- ou_result$half_lives_X[is.finite(ou_result$half_lives_X)]
          hl_x <- hl_x[hl_x > 0 & hl_x < 1000]
          hl_ou <- c(hl_ou, hl_x)
        }
        if (!is.null(ou_result$half_lives_Y)) {
          hl_y <- ou_result$half_lives_Y[is.finite(ou_result$half_lives_Y)]
          hl_y <- hl_y[hl_y > 0 & hl_y < 1000]
          hl_ou <- c(hl_ou, hl_y)
        }
      }
      
      if (length(hl_ou) > 0) {
        n_x <- sum(!is.null(ou_result$half_lives_X))
        n_y <- sum(!is.null(ou_result$half_lives_Y))
        
        bar_names <- if (length(hl_ou) <= 10) {
          c(if(n_x > 0) paste0("X", 1:length(ou_result$half_lives_X[is.finite(ou_result$half_lives_X)])),
            if(n_y > 0) paste0("Y", 1:length(ou_result$half_lives_Y[is.finite(ou_result$half_lives_Y)])))
        } else {
          rep("", length(hl_ou))
        }
        
        cols <- c(rep("lightblue", length(ou_result$half_lives_X[is.finite(ou_result$half_lives_X)])),
                  rep("lightcoral", length(ou_result$half_lives_Y[is.finite(ou_result$half_lives_Y)])))
        
        bp <- barplot(hl_ou, 
                      names.arg = bar_names[1:length(hl_ou)],
                      main = "Half-lives (Factor-OU)", 
                      ylab = "Periods", 
                      col = cols[1:length(hl_ou)],
                      las = 1,
                      ylim = c(0, max(hl_ou) * 1.1))
        
        if (length(hl_ou) <= 10) {
          text(bp, hl_ou + max(hl_ou) * 0.02, 
               round(hl_ou, 1), 
               cex = 0.8, pos = 3)
        }
        
        abline(h = median(hl_ou), col = "red", lty = 2, lwd = 2)
        
        mtext(paste("Median:", round(median(hl_ou), 1), "periods"), 
              side = 3, line = 0, cex = 0.8, col = "red")
        
        if (n_x > 0 && n_y > 0) {
          legend("topright", 
                 legend = c("Factores X", "Factores Y"),
                 fill = c("lightblue", "lightcoral"),
                 bty = "n", cex = 0.9)
        }
      } else {
        plot(1, type = "n", axes = FALSE, xlab = "", ylab = "")
        text(1, 1, "Sin half-lives OU", cex = 1.2, col = "gray50")
        box()
        title("Half-lives OU")
      }
    },
    "Half-lives (OU Factor Model)",
    "Half-lives OU not available"
  )
  if (plot4_ok) plots_created <- plots_created + 1
  
  tryCatch({
    if (!save_plot) {
      mtext("Factorial Dynamics Analysis", 
            side = 3, line = -1, outer = TRUE, cex = 1.2, font = 2)
    }
  }, error = function(e) {})
  
  if (verbose) {
    message(sprintf("  OK: Full display: %d of/4 panels created", plots_created))
  }
  
  return(invisible(plots_created == 4))
}