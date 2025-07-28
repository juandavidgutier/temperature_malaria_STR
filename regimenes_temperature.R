library(data.table)
library(haldensify)
library(sl3)
library(tmle3)
library(tmle3shift)
library(dplyr)
library(tidyr)
library(ggplot2)

# --- Data Loading ---
data_all <- read.csv("D:/data07-23.csv")

# Apply temperature filter
dataset_raw <- data_all[data_all$Temperature >= 15 & data_all$Temperature <= 30, ]

# Select relevant columns (assuming columns 7 to 28 are the ones needed)
# *** CRITICAL: Verify column indices/names ***
# selected_columns <- names(data_all)[7:28] # Get names first
dataset_selected <- dataset_raw[, 7:28] # Then select

# --- Explicitly Handle Missing Values First ---
# Check for NAs before processing
cat("Checking for missing values in selected columns before processing:\n")
print(colSums(is.na(dataset_selected)))

# Drop rows with ANY NAs in the selected columns upfront
dataset_clean_na <- dataset_selected %>% drop_na()
cat("Rows before NA removal:", nrow(dataset_selected), "\n")
cat("Rows after NA removal:", nrow(dataset_clean_na), "\n")

# --- Convert all columns to numeric explicitly ---
# This is a robust way to ensure all columns are numeric
# It will coerce factors/characters to numeric (NAs if coercion fails)
# and warn if there are issues.
cat("Converting all selected columns to numeric...\n")
dataset_numeric <- dataset_clean_na %>%
  mutate(across(everything(), as.numeric))

# Check if conversion introduced any new NAs (shouldn't happen if data was clean numeric/char)
conversion_na_check <- colSums(is.na(dataset_numeric))
if (any(conversion_na_check > 0)) {
  cat("!!! WARNING: Converting to numeric introduced NAs in the following columns: !!!\n")
  print(conversion_na_check[conversion_na_check > 0])
  # You might want to stop here if this happens
  # stop("Data type conversion failed, introducing NAs.")
} else {
  cat("All columns successfully converted to numeric.\n")
}

# --- Binarize Confounders (Columns 3 to 13 of the SELECTED dataset) ---
# *** CRITICAL: Ensure indices 3:13 are correct for YOUR selected columns ***
# Let's check the names first to be sure
confounder_names <- names(dataset_numeric)[3:14]
cat("Binarizing confounders:", paste(confounder_names, collapse = ", "), "\n")

for (col_name in confounder_names) {
  # Check if column exists and is numeric
  if (!(col_name %in% names(dataset_numeric))) {
    stop(paste("Column", col_name, "not found in dataset_numeric. Check indices."))
  }
  if (!is.numeric(dataset_numeric[[col_name]])) {
    stop(paste("Column", col_name, "is not numeric. Conversion failed."))
  }
  
  mediana <- median(dataset_numeric[[col_name]], na.rm = TRUE)
  # Explicitly convert result to numeric (shouldn't be needed, but defensive)
  dataset_numeric[[col_name]] <- as.numeric(ifelse(dataset_numeric[[col_name]] > mediana, 1, 0))
}

# Final dataset after all processing
dataset <- dataset_numeric

# --- Final Data Verification ---
cat("\n--- Final Dataset Structure ---\n")
str(dataset)
cat("\n--- Summary of Key Variables ---\n")
print(summary(dataset[c("Temperature", "excess", confounder_names[1])])) # Check a couple

# learners used for conditional mean of the outcome
bart_lrnr <- Lrnr_dbarts$new()
earth_lrnr <- Lrnr_earth$new()
rf_lrnr <- Lrnr_ranger$new()
hal_lrnr <- Lrnr_hal9001$new(max_degree = 3)

# SL for the outcome regression
sl_reg_lrnr <- Lrnr_sl$new(
  learners = list(bart_lrnr, earth_lrnr, rf_lrnr, hal_lrnr),
  metalearner = Lrnr_nnls$new()
)

sl3_list_learners("density")

# learners used for conditional densities for (g_n)
haldensify_lrnr <- Lrnr_haldensify$new(
  n_bins = c(5, 10, 20),
  lambda_seq = exp(seq(-1, -10, length = 200))
)
# semiparametric density estimator with homoscedastic errors (HOSE)
hose_hal_lrnr <- make_learner(Lrnr_density_semiparametric,
                              mean_learner = hal_lrnr
)
# semiparametric density estimator with heteroscedastic errors (HESE)
hese_rf_glm_lrnr <- make_learner(Lrnr_density_semiparametric,
                                 mean_learner = rf_lrnr,
                                 var_learner = earth_lrnr
)

# SL for the conditional treatment density
sl_dens_lrnr <- Lrnr_sl$new(
  learners = list(hose_hal_lrnr, hese_rf_glm_lrnr),
  metalearner = Lrnr_solnp_density$new()
)

learner_list <- list(Y = sl_reg_lrnr, A = sl_dens_lrnr)

# --- Node list ---
cat("\n--- Final Column Names in 'dataset' ---\n")
print(names(dataset))

node_list <- list(
  W = names(dataset)[1:14], # Explicitly use column names
  A = "Temperature",        # Must match exactly
  Y = "excess"              # Must match exactly
)

# --- Ensure Node Variables are Present ---
missing_nodes <- setdiff(c(node_list$A, node_list$Y, node_list$W), names(dataset))
if (length(missing_nodes) > 0) {
  stop(paste("The following nodes are missing from dataset:", paste(missing_nodes, collapse = ", ")))
}

# --- Function to fit TMLE and store results ---
fit_tmle_and_store <- function(shift_val, index) {
  cat("  -> Fitting TMLE for shift =", shift_val, "...\n")
  tmle_spec <- tmle_shift(shift_val = shift_val, shift_fxn = shift_additive, shift_fxn_inv = shift_additive_inv)

  tmle_fit <- tryCatch({
    tmle3(tmle_spec, dataset, node_list, learner_list)
  }, error = function(e) {
    cat("  -> ERROR fitting shift =", shift_val, ":", conditionMessage(e), "\n")
    # Print problematic columns if error is related to data types
    if (grepl("numeric|matrix", conditionMessage(e))) {
      cat("     -> Suspect data type issue. Checking column classes:\n")
      problematic_cols <- c(node_list$W[1:min(5, length(node_list$W))], node_list$A, node_list$Y)
      for(col in problematic_cols) {
        cat("        ", col, ":", class(dataset[[col]]), "\n")
      }
    }
    return(NULL)
  })
  
  if (!is.null(tmle_fit)) {
    main_results[index, "EY"] <- tmle_fit$summary$tmle_est
    main_results[index, "Lower"] <- tmle_fit$summary$lower
    main_results[index, "Upper"] <- tmle_fit$summary$upper
  } else {
    main_results[index, "EY"] <- NA
    main_results[index, "Lower"] <- NA
    main_results[index, "Upper"] <- NA
  }
  return(tmle_fit)
}

# --- Dataframe for saving results ---
main_results <- data.frame(
  `labs` = c("Observed temperature", "Temperature + 0.5 °C", "Temperature + 1.0 °C",
             "Temperature + 1.5 °C", "Temperature + 2.0 °C"),
  EY = numeric(5),
  Lower = numeric(5),
  Upper = numeric(5),
  ATE = numeric(5),
  ATE_Lower = numeric(5),
  ATE_Upper = numeric(5)
)
main_results[1, "ATE"] <- 0
main_results[1, "ATE_Lower"] <- 0
main_results[1, "ATE_Upper"] <- 0

# --- Fit models ---
cat("\n--- Starting TMLE Fits ---\n")
fit_observed <- fit_tmle_and_store(0, 1)
fit_05 <- fit_tmle_and_store(0.5, 2)
fit_10 <- fit_tmle_and_store(1.0, 3)
fit_15 <- fit_tmle_and_store(1.5, 4)
fit_20 <- fit_tmle_and_store(2.0, 5)


# --- Dataframe for saving results ---
main_results <- data.frame(
  `labs` = c("Observed temperature", "Temperature + 0.5 °C", "Temperature + 1.0 °C",
             "Temperature + 1.5 °C", "Temperature + 2.0 °C"),
  EY = c(
    if (!is.null(fit_observed)) fit_observed$summary$tmle_est else NA,
    if (!is.null(fit_05)) fit_05$summary$tmle_est else NA,
    if (!is.null(fit_10)) fit_10$summary$tmle_est else NA,
    if (!is.null(fit_15)) fit_15$summary$tmle_est else NA,
    if (!is.null(fit_20)) fit_20$summary$tmle_est else NA
  ),
  Lower = c(
    if (!is.null(fit_observed)) fit_observed$summary$lower else NA,
    if (!is.null(fit_05)) fit_05$summary$lower else NA,
    if (!is.null(fit_10)) fit_10$summary$lower else NA,
    if (!is.null(fit_15)) fit_15$summary$lower else NA,
    if (!is.null(fit_20)) fit_20$summary$lower else NA
  ),
  Upper = c(
    if (!is.null(fit_observed)) fit_observed$summary$upper else NA,
    if (!is.null(fit_05)) fit_05$summary$upper else NA,
    if (!is.null(fit_10)) fit_10$summary$upper else NA,
    if (!is.null(fit_15)) fit_15$summary$upper else NA,
    if (!is.null(fit_20)) fit_20$summary$upper else NA
  ),
  ATE = numeric(5),
  ATE_Lower = numeric(5),
  ATE_Upper = numeric(5)
)

# ATE vs itself is 0
main_results[1, "ATE"] <- 0
main_results[1, "ATE_Lower"] <- 0
main_results[1, "ATE_Upper"] <- 0

# ATE 
calculate_ate_manual <- function(index_shift, index_obs = 1) {
  if (is.na(main_results[index_shift, "EY"]) || is.na(main_results[index_obs, "EY"])) {
    warning("can not estimate ATE: Fault EY values.")
    return(c(ATE = NA, Lower = NA, Upper = NA))
  }
  
  psi_obs <- main_results[index_obs, "EY"]
  psi_shift <- main_results[index_shift, "EY"]
  
  # ATE as difference
  ate_est <- psi_shift - psi_obs
  
  # Standard errors
  se_obs <- (main_results[index_obs, "Upper"] - main_results[index_obs, "Lower"]) / (2 * 1.96)
  se_shift <- (main_results[index_shift, "Upper"] - main_results[index_shift, "Lower"]) / (2 * 1.96)
  
  var_ate <- se_shift^2 + se_obs^2
  se_ate <- sqrt(var_ate)
  
  # Confidence interval for ATE
  margin_of_error <- 1.96 * se_ate
  ate_ci_lower <- ate_est - margin_of_error
  ate_ci_upper <- ate_est + margin_of_error
  
  return(c(ATE = ate_est, Lower = ate_ci_lower, Upper = ate_ci_upper))
}

# --- Calculate and Store ATEs  ---
manual_ate_05 <- calculate_ate_manual(2, 1) # +0.5 vs Observed
main_results[2, c("ATE", "ATE_Lower", "ATE_Upper")] <- manual_ate_05

manual_ate_10 <- calculate_ate_manual(3, 1) # +1.0 vs Observed
main_results[3, c("ATE", "ATE_Lower", "ATE_Upper")] <- manual_ate_10

manual_ate_15 <- calculate_ate_manual(4, 1) # +1.5 vs Observed
main_results[4, c("ATE", "ATE_Lower", "ATE_Upper")] <- manual_ate_15

manual_ate_20 <- calculate_ate_manual(5, 1) # +2.0 vs Observed
main_results[5, c("ATE", "ATE_Lower", "ATE_Upper")] <- manual_ate_20

# --- Final results ---
print(main_results)


# --- Figure 3: Plot ATEs ---
ate_plot_data <- main_results[-1, , drop = FALSE] 
ate_plot_data <- ate_plot_data[!is.na(ate_plot_data$ATE), , drop = FALSE]

if (nrow(ate_plot_data) > 0) {
  ate_plot_data$labs <- factor(ate_plot_data$labs, levels = ate_plot_data$labs)
  
  p_ate <- ggplot(ate_plot_data, aes(x = labs, y = ATE)) +
    geom_point(size = 4, color = "darkblue") +
    geom_errorbar(aes(ymin = ATE_Lower, ymax = ATE_Upper), width = 0.2, color = "darkblue") +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red", linewidth = 0.8) +
    labs(
      x = "Interventions",
      y = "ATE"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 12)
    )
  
  # Print
  print(p_ate)
  
  
} else {
  cat("\n!!!No ATEs valid. !!!\n")
}
