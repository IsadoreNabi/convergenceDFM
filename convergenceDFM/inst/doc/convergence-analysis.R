## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(convergenceDFM)

## ----example, eval=FALSE------------------------------------------------------
# # Load example data
# data("example_marxist_data")
# 
# # Run complete analysis
# results <- run_complete_factor_analysis_robust(
#   X_matrix = marxist_prices[, -1],
#   Y_matrix = bayesian_cpi[, -1],
#   max_comp = 3,
#   dfm_lags = 1,
#   ou_chains = 4,
#   ou_iter = 2000
# )
# 
# # View results
# summary(results)

## ----viz, eval=FALSE----------------------------------------------------------
# # Visualize factor dynamics
# visualize_factor_dynamics(
#   dfm_result = results$dfm,
#   ou_result = results$factor_ou,
#   factors_data = results$factors
# )

