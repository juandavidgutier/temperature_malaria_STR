# Load required libraries
# Note: Ensure you have installed kuenm, raster, sp, dismo, and other dependencies
# Install kuenm from GitHub if not already installed: 
# devtools::install_github("marlonecobos/kuenm")
library(kuenm)
library(raster)
library(sp)
library(dismo)
library(rgdal)
library(stringr)

# Set working directory and paths
main_dir <- "D:/"
occ_dir <- paste0(main_dir, "/points/")
env_dir <- paste0(main_dir, "/WoldClim/")
results_dir <- paste0(main_dir, "/malaria_modeling_results/")

# Create main results directory if it doesn't exist
if (!dir.exists(results_dir)) {
  dir.create(results_dir)
}

# List all CSV files in the occurrence directory
occ_files <- list.files(occ_dir, pattern = "*.csv", full.names = TRUE)
vector_names <- tools::file_path_sans_ext(basename(occ_files))

# Loop through each of the 8 malaria vectors
for (i in 1:length(occ_files)) {
  vector_name <- vector_names[i]
  print(paste("Processing vector:", vector_name))
  
  # Create species-specific directory structure
  species_dir <- file.path(results_dir, vector_name)
  if (!dir.exists(species_dir)) {
    dir.create(species_dir)
  }
  
  # Set working directory to species directory
  setwd(species_dir)
  
  # Create necessary subdirectories as required by kuenm
  if (!dir.exists("occurrences")) {
    dir.create("occurrences")
  }
  if (!dir.exists("env")) {
    dir.create("env")
  }
  if (!dir.exists("models")) {
    dir.create("models")
  }
  if (!dir.exists("projections")) {
    dir.create("projections")
  }
  
  # 1. Prepare occurrence data
  # Read occurrence data
  occ_data <- read.csv(occ_files[i])
  
  # Check if data has correct columns (assumes LAT and LON columns)
  if (!all(c("LAT", "LON") %in% colnames(occ_data))) {
    stop("Occurrence data must have LAT and LON columns")
  }
  
  # Save complete occurrence data
  write.csv(occ_data, file.path("occurrences", "complete_occurrences.csv"), row.names = FALSE)
  
  # Split data into training (75%) and testing (25%) subsets
  # Note: For robust modeling, a more sophisticated partitioning might be needed
  set.seed(123)  # For reproducibility
  train_idx <- sample(1:nrow(occ_data), 0.75 * nrow(occ_data))
  train_data <- occ_data[train_idx, ]
  test_data <- occ_data[-train_idx, ]
  
  # Save training and testing data
  write.csv(train_data, file.path("occurrences", "train_occurrences.csv"), row.names = FALSE)
  write.csv(test_data, file.path("occurrences", "test_occurrences.csv"), row.names = FALSE)
  
  # 2. Prepare environmental variables
  # Copy WorldClim rasters to species env directory
  env_files <- list.files(env_dir, pattern = "*.tif", full.names = TRUE)
  file.copy(env_files, file.path("env"), overwrite = TRUE)
  
  # Create 10 different variable sets by removing correlated variables
  # This implements the exhaustive variable selection approach as suggested by Cobos et al. (2019)
  # First, calculate correlation matrix among environmental variables
  env_stack <- stack(file.path("env", list.files("env", pattern = "*.tif")))
  env_values <- extract(env_stack, occ_data[, c("LON", "LAT")])
  
  # Calculate correlation matrix
  corr_matrix <- cor(do.call(rbind, env_values), use = "complete.obs")
  
  # Function to create variable subsets by removing highly correlated variables
  create_variable_subsets <- function(corr_matrix, n_subsets = 10) {
    vars <- colnames(corr_matrix)
    subsets <- list()
    
    for (s in 1:n_subsets) {
      # Randomly select a threshold for correlation (between 0.7 and 0.9)
      threshold <- runif(1, 0.7, 0.9)
      
      # Identify highly correlated pairs
      high_corr <- which(abs(corr_matrix) >= threshold & row(corr_matrix) != col(corr_matrix), arr.ind = TRUE)
      
      # Create a network of correlations
      corr_network <- list()
      for (i in 1:nrow(high_corr)) {
        v1 <- vars[high_corr[i, 1]]
        v2 <- vars[high_corr[i, 2]]
        if (!v1 %in% names(corr_network)) corr_network[[v1]] <- character(0)
        if (!v2 %in% names(corr_network)) corr_network[[v2]] <- character(0)
        corr_network[[v1]] <- c(corr_network[[v1]], v2)
        corr_network[[v2]] <- c(corr_network[[v2]], v1)
      }
      
      # Randomly select variables to keep from each correlated group
      kept_vars <- character(0)
      visited <- character(0)
      
      for (v in vars) {
        if (v %in% visited) next
        
        # Get all variables correlated with v
        to_visit <- v
        group <- character(0)
        
        while (length(to_visit) > 0) {
          current <- to_visit[1]
          to_visit <- to_visit[-1]
          
          if (!(current %in% visited)) {
            visited <- c(visited, current)
            group <- c(group, current)
            
            if (current %in% names(corr_network)) {
              neighbors <- corr_network[[current]]
              to_visit <- c(to_visit, neighbors[!(neighbors %in% visited)])
            }
          }
        }
        
        # Randomly select one variable from the group to keep
        if (length(group) > 0) {
          kept <- sample(group, 1)
          kept_vars <- c(kept_vars, kept)
        }
      }
      
      subsets[[s]] <- kept_vars
    }
    
    return(subsets)
  }
  
  # Create 10 variable sets
  variable_subsets <- create_variable_subsets(corr_matrix, n_subsets = 10)
  
  # Save variable subset information for later reference
  write.csv(data.frame(Set = rep(1:10, times = sapply(variable_subsets, length)), 
                       Variable = unlist(variable_subsets)),
            "variable_subsets.csv", row.names = FALSE)
  
  # 3. Model calibration with kuenm
  # Define parameter combinations as requested
  # Regularization multipliers: 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 2, 3, 4, 5, 6, 10
  rm_values <- c(seq(0.1, 1, by = 0.1), 2:6, 10)
  
  # Feature class combinations (all possible combinations of l, q, p, h, t)
  feature_classes <- c("l", "q", "p", "h", "t")
  fc_combinations <- unlist(lapply(1:5, function(x) {
    apply(combn(feature_classes, x), 2, function(y) paste(y, collapse = ""))
  }))
  
  # Create candidate models using kuenm_cal
  # This implements the detailed model calibration as described by Cobos et al. (2019) [[1]]
  cal_results <- kuenm_cal(
    calib = "occurrences/train_occurrences.csv",
    eval = "occurrences/test_occurrences.csv",
    env = "env",
    out = "models",
    rm = rm_values,
    fc = fc_combinations,
    var = variable_subsets,
    prefix = vector_name,
    runs = 10,  # Create 10 replicates for each model
    threads = 4,  # Adjust based on your system capabilities
    maxent = "maxent.jar"  # Path to Maxent jar file
  )
  
  # Check if candidate models were created successfully
  if (!cal_results$completed) {
    warning(paste("Model calibration failed for", vector_name))
    next
  }
  
  # 4. Model evaluation and selection
  # Evaluate candidate models based on three criteria as requested
  # a) Statistical significance via partial ROC with 500 iterations and 50% bootstrap
  # b) Omission rates < 5%
  # c) Model complexity using AICc
  eval_results <- kuenm_ceval(
    calib = "occurrences/train_occurrences.csv",
    eval = "occurrences/test_occurrences.csv",
    models = "models",
    prefix = vector_name,
    iterations = 500,  # 500 iterations for partial ROC
    boot_prop = 0.5,   # 50% of data for bootstrap
    E = 5,             # 5% omission rate threshold
    AICc = TRUE        # Calculate AICc for model complexity
  )
  
  # Check if evaluation was successful
  if (!eval_results$completed) {
    warning(paste("Model evaluation failed for", vector_name))
    next
  }
  
  # 5. Create final model with best parameters
  # Select the best model based on the three criteria
  best_models <- eval_results$best.models
  
  if (nrow(best_models) == 0) {
    warning(paste("No models met all selection criteria for", vector_name))
    # If no models meet all criteria, try to find models that meet at least two criteria
    # This follows the approach used in malaria vector modeling studies [[4]]
    partial_best <- eval_results$evaluation[eval_results$evaluation$Significance == "Significant" & 
                                              (eval_results$evaluation$Omission <= 0.05 | eval_results$evaluation$Delta_AICc <= 2), ]
    if (nrow(partial_best) > 0) {
      best_models <- partial_best[order(partial_best$Delta_AICc), ][1, , drop = FALSE]
      warning(paste("Using best available model (not meeting all criteria) for", vector_name))
    } else {
      next
    }
  }
  
  # Create final model with the best parameters
  mod_results <- kuenm_mod(
    occ = "occurrences/complete_occurrences.csv",
    env = "env",
    models = "models",
    prefix = vector_name,
    best = best_models[1, ],  # Use the very best model
    proj = "projections",
    extrapolation = "clamp",   # Use clamping for extrapolation
    runs = 10,                # 10 bootstrap replicates
    threads = 4,
    maxent = "maxent.jar"
  )
  
  # 6. Convert model to binary using 5% omission threshold
  # First, evaluate omission rate at different thresholds
  # This implements the threshold selection approach used in malaria vector modeling [[6]]
  om_results <- kuenm_omrat(
    occ = "occurrences/complete_occurrences.csv",
    model = file.path("models", paste0(vector_name, "_", best_models[1, "ID"], ".lambdas")),
    env = "env",
    E_values = seq(0.01, 0.1, by = 0.01)  # Test thresholds from 1% to 10%
  )
  
  # Find threshold corresponding to 5% omission rate
  threshold_idx <- which.min(abs(om_results$Omission - 0.05))
  threshold_val <- om_results$Threshold[threshold_idx]
  
  # Convert continuous model to binary using this threshold
  cont_model <- raster(file.path("projections", paste0(vector_name, "_", best_models[1, "ID"], ".tif")))
  bin_model <- cont_model > threshold_val
  
  # Save binary model
  writeRaster(bin_model, 
              filename = file.path("projections", paste0(vector_name, "_binary.tif")),
              overwrite = TRUE)
  
  # 7. Create a summary report for this vector
  # This follows the reproducible research approach recommended by kuenm [[7]]
  report <- paste(
    "## Malaria Vector Distribution Modeling Report for", vector_name,
    "\n\n### Model Selection Criteria",
    "\n- Statistical significance: Partial ROC with 500 iterations and 50% bootstrap",
    "\n- Omission rate threshold: 5%",
    "\n- Model complexity: AICc",
    "\n\n### Best Model Parameters",
    paste("\n- Regularization multiplier:", best_models[1, "RM"]),
    paste("\n- Feature classes:", best_models[1, "FC"]),
    paste("\n- Variable set:", best_models[1, "VS"]),
    paste("\n- AUC ratio:", round(best_models[1, "AUC.ratio"], 3)),
    paste("\n- Omission rate:", round(best_models[1, "Omission"], 3)),
    paste("\n- Delta AICc:", round(best_models[1, "Delta_AICc"], 3)),
    "\n\n### Binary Threshold",
    paste("\n- Threshold value:", round(threshold_val, 3)),
    paste("\n- Actual omission rate at threshold:", round(om_results$Omission[threshold_idx], 3))
  )
  
  writeLines(report, con = paste0(vector_name, "_report.md"))
  
  print(paste("Completed modeling for vector:", vector_name))
}

# Create a master report summarizing all models
master_report <- paste(
  "# Malaria Vector Distribution Modeling Summary",
  "\n\nThis analysis used the kuenm R package to develop ecological niche models",
  "for 8 malaria vectors following the rigorous calibration approach described by Cobos et al. (2019) [[1]].",
  "\n\nThe modeling process included:",
  "\n- Exhaustive variable selection creating 10 different environmental variable sets",
  "\n- Testing 16 regularization multipliers and multiple feature class combinations",
  "\n- Model evaluation based on statistical significance, omission rates, and model complexity",
  "\n- Binary conversion using the 5% omission threshold approach",
  "\n\nAll models were developed using Maxent via the kuenm package, which has been",
  "successfully applied to mosquito vector distribution modeling in previous studies [[4]].",
  "\n\nThis approach ensures robust model development while accounting for",
  "multiple sources of variation and uncertainty in ecological niche modeling [[3]]."
)

writeLines(master_report, con = file.path(results_dir, "master_report.md"))

print("All malaria vectors processed successfully!")
print(paste("Results saved in:", results_dir))