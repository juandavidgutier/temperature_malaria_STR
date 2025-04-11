library(data.table)
library(haldensify)
library(sl3)
library(tmle3)
library(tmle3shift)
library(dplyr)
library(tidyr)
library(ggplot2)


#Dataframe for saving results
main_results <- data.frame(`labs` = c("Observed temperature", "Temperature + 0.5 °C", "Temperature + 1.0 °C", 
                                      "Temperature + 1.5 °C", "Temperature + 2.0 °C"),
                           ATE = numeric(5),  
                           Lower = numeric(5),    
                           Upper = numeric(5))

data_all <- read.csv("D:/data07-23.csv")
str(data_all)

# Temperature between 15-30°C
dataset <- data_all[data_all$Temperature >= 15 & data_all$Temperature <= 30, ]

dataset <- data_all[,7:28]

# Confounders as binary 
for (i in 3:13) {
  mediana <- median(dataset[[i]], na.rm = TRUE)
  dataset[[i]] <- ifelse(dataset[[i]] > mediana, 1, 0)
}
str(dataset)

dataset <- dataset %>% drop_na()

# learners used for conditional mean of the outcome
bart_lrnr <- Lrnr_dbarts$new()
earth_lrnr <- Lrnr_earth$new()
rf_lrnr <- Lrnr_ranger$new()
hal_lrnr <- Lrnr_hal9001$new(max_degree = 3, n_folds = 3)

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


#############################################################################################
#node list
node_list <- list(
  W = names(dataset[,1:14]), #covariates
  A = "Temperature", #exposure
  Y = "excess" #outcome
)

# initialize a tmle specification for currently observed temperature
tmle_spec0 <- tmle_shift(
  shift_val = 0.0001, #We specify the argument shift_val = 0.0001°
  shift_fxn = shift_additive,
  shift_fxn_inv = shift_additive_inv
)

# Population current effect:
tmle_fit0 <- tmle3(tmle_spec0, dataset, node_list, learner_list)
print(tmle_fit0)

main_results[1,2] <- tmle_fit0$summary$tmle_est
main_results[1,3] <- tmle_fit0$summary$lower
main_results[1,4] <- tmle_fit0$summary$upper

# initialize a tmle specification for 0.5°
tmle_spec05 <- tmle_shift(
  shift_val = 0.5, #We specify the argument shift_val = 0.5°
  shift_fxn = shift_additive,
  shift_fxn_inv = shift_additive_inv
)

#Population effect for 0.5°:
tmle_fit05 <- tmle3(tmle_spec05, dataset, node_list, learner_list)
print(tmle_fit05)

main_results[2,2] <- tmle_fit05$summary$tmle_est
main_results[2,3] <- tmle_fit05$summary$lower
main_results[2,4] <- tmle_fit05$summary$upper

# initialize a tmle specification for 1.0°
tmle_spec1 <- tmle_shift(
  shift_val = 1.0, #We specify the argument shift_val = 1.0°
  shift_fxn = shift_additive,
  shift_fxn_inv = shift_additive_inv
)

#Population effect for 1°:
tmle_fit1 <- tmle3(tmle_spec1, dataset, node_list, learner_list)
print(tmle_fit1)

main_results[3,2] <- tmle_fit1$summary$tmle_est
main_results[3,3] <- tmle_fit1$summary$lower
main_results[3,4] <- tmle_fit1$summary$upper

# initialize a tmle specification for 1.5°
tmle_spec15 <- tmle_shift(
  shift_val = 1.5, #We specify the argument shift_val = 1.5° 
  shift_fxn = shift_additive,
  shift_fxn_inv = shift_additive_inv
)

#Population effect for 1.5°:
tmle_fit15 <- tmle3(tmle_spec15, dataset, node_list, learner_list)
print(tmle_fit15)

main_results[4,2] <- tmle_fit15$summary$tmle_est
main_results[4,3] <- tmle_fit15$summary$lower
main_results[4,4] <- tmle_fit15$summary$upper

# initialize a tmle specification for 2.0°
tmle_spec2 <- tmle_shift(
  shift_val = 2.0, #We specify the argument shift_val = 2.0°
  shift_fxn = shift_additive,
  shift_fxn_inv = shift_additive_inv
)

#Population effect for 2.0°:
tmle_fit2 <- tmle3(tmle_spec2, dataset, node_list, learner_list)
print(tmle_fit2)

main_results[5,2] <- tmle_fit2$summary$tmle_est
main_results[5,3] <- tmle_fit2$summary$lower
main_results[5,4] <- tmle_fit2$summary$upper


# Define the desired order for the labs field
desired_order <- c("Temperature + 2.0 °C", "Temperature + 1.5 °C", "Temperature + 1.0 °C", "Temperature + 0.5 °C", "Observed temperature")

# Convert labs to a factor with the specified levels to control ordering
main_results$labs <- factor(main_results$labs, levels = desired_order)

# Fig. 3
ggplot(main_results, aes(y = labs, x = ATE)) +
  geom_point(color = "red", size = 2.5, shape = 15) +
  geom_errorbarh(aes(xmin = Lower, xmax = Upper), height = 0.0, lineend = "butt") +
  geom_text(aes(x = 0.95, label = sprintf("%.3f", ATE)), hjust = 0) +
  geom_text(aes(x = 1.05, label = sprintf("(%.3f, %.3f)", Lower, Upper)), hjust = 0) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  coord_cartesian(xlim = c(0.6, 1.1)) +
  scale_x_continuous(expand = c(0.5, 0.05)) +
  annotate("text", x = 1.08, y = max(as.numeric(main_results$labs)) + 0.5, 
           label = "ATE (95% CI)", hjust = 0.5, fontface = "bold") +
  
  # Clean up the theme
  theme_minimal() +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title = element_blank(),
    axis.text.y = element_text(hjust = 0, size = 12), 
    axis.text.x = element_text(hjust = 0, size = 12),
    plot.margin = margin(r = 3, l = 5, t = 40, unit = "pt")  
  )


