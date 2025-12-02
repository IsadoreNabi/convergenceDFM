
# Declaración de variables globales para evitar NOTEs de R CMD check
# Estas variables son usadas dentro de pipes de dplyr/tidyverse
utils::globalVariables(c(".", "Year", "CPI", "Industry", "Weight_raw", "Weight"))
