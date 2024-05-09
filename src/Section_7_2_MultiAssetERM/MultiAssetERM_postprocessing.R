################################################################################
# This code is for postprocessing the results of multi-asset ERM experiments   #
################################################################################

rm(list = ls()) # clear all variables

#-- load and combine the results
# some codes to create a list of file names, e.g., filename_list, to load
load_filename_list <- list.files(pattern = "MultiAssetERM_10Stocks1000Outer10Macro")
num_files <- length(load_filename_list) # total number of files to load
load(file = load_filename_list[1]) # load the first file to get num_macro

#-- placeholders for combined results. All combined results have _all suffix
num_macro_all <- num_macro * num_files # total number of macro scenarios

# all outer scenarios
s_tau_mac_all <- array(0, dim = c(num_outer, num_stocks, num_macro_all))

# closed-form portfolio values
val_true_all <- matrix(0, nrow = num_outer, ncol = num_macro_all)

# standard nested simulation with n_inner in each outer scenario
val_nested_plus_all <- matrix(0, nrow = num_outer, ncol = num_macro_all)
time_nested_plus_all <- rep(0, 2) # total time for nested plus simulation

# optimal LR design
val_opt_nested_all <- matrix(0, nrow = num_outer, ncol = num_macro_all)
sim_budgets_all <- matrix(0, nrow = num_stocks, ncol = num_macro_all)
total_num_samp_all <- matrix(0, nrow = num_stocks, ncol = num_macro_all)
time_opt_nested_opt_all <- rep(0, 2) # runtime for optimization
time_opt_nested_sim_all <- rep(0, 2) # runtime for simulation
time_opt_nested_tot_all <- rep(0, 2) # total runtime

# regression-based nested simulation
val_reg_nested_all <- matrix(0, nrow = num_outer, ncol = num_macro_all)
time_reg_nested_sim_all <- rep(0, 2) # runtime for simulation
time_reg_nested_fit_all <- rep(0, 2) # runtime for fitting
time_reg_nested_tot_all <- rep(0, 2) # total runtime
s_tau_reg_nested_all <- vector("list", num_macro_all) # outer scens

# standard nested Sim with Gamma^(1/3) inner reps & Gamma^(2/3) outer scen
val_true_nested_all <- vector("list", num_macro_all) # true value for MSEs
val_std_nested_all <- vector("list", num_macro_all) # standard nested estimates
time_std_nested_all <- rep(0, 2) # runtime for standard nested simulation
s_tau_std_nested_all <- vector("list", num_macro_all) # outer scens
sim_budget_std_nested_all <- matrix(0, nrow = 3, ncol = num_macro_all)

#-- load and combine the results
for (fileID in 1:num_files) {
  # load the results
  filename <- load_filename_list[[fileID]]
  load(filename)

  # calculate the indices for the current file
  macro_start <- (fileID - 1) * num_macro + 1 # start index
  macro_end <- fileID * num_macro # end index
  macro_index <- c(macro_start:macro_end) # indices for the current file

  s_tau_mac_all[, , macro_index] <- s_tau_mac # outer scenarios
  val_true_all[, macro_index] <- val_true # true values

  # standard nested simulation with n_inner in each outer scenario
  val_nested_plus_all[, macro_index] <- val_nested_plus
  time_nested_plus_all <- time_nested_plus_all + time_nested_plus

  # optimal LR design
  val_opt_nested_all[, macro_index] <- val_opt_nested
  sim_budgets_all[, macro_index] <- sim_budgets
  total_num_samp_all[, macro_index] <- total_num_samp
  time_opt_nested_opt_all <- time_opt_nested_opt_all + time_opt_nested_opt
  time_opt_nested_sim_all <- time_opt_nested_sim_all + time_opt_nested_sim
  time_opt_nested_tot_all <- time_opt_nested_tot_all + time_opt_nested_tot

  # regression-based nested simulation
  val_reg_nested_all[, macro_index] <- val_reg_nested
  time_reg_nested_sim_all <- time_reg_nested_sim_all + time_reg_nested_sim
  time_reg_nested_fit_all <- time_reg_nested_fit_all + time_reg_nested_fit
  time_reg_nested_tot_all <- time_reg_nested_tot_all + time_reg_nested_tot
  s_tau_reg_nested_all[macro_index] <- s_tau_reg_nested

  # standard nested Sim with Gamma^(1/3) inner reps & Gamma^(2/3) outer scen
  val_true_nested_all[macro_index] <- val_true_nested
  val_std_nested_all[macro_index] <- val_std_nested
  time_std_nested_all <- time_std_nested_all + time_std_nested
  s_tau_std_nested_all[macro_index] <- s_tau_std_nested
  sim_budget_std_nested_all[, macro_index] <- sim_budget_std_nested

  print(sprintf(
    "Processing %d-th result file, max value is %6.2f",
    fileID, max(val_nested_plus_all)
  ))
}

#-- calculate simulation budgets (first row of Table 2)
budget_nested_plus <- num_inner * num_outer
budget_opt_nested <- mean(sim_budgets_all)
budget_reg_nested <- mean(ceiling(colMeans(sim_budgets_all) / num_inner_reg))
budget_std_nested <- mean(sim_budget_std_nested_all[3, ])
budget <- c(
  budget_nested_plus, budget_opt_nested,
  budget_std_nested, budget_reg_nested
)

#-- calculate runtime (second row of Table 2)
runtime <- cbind(
  time_nested_plus_all[1], time_opt_nested_tot_all[1],
  time_std_nested_all[1], time_reg_nested_tot_all[1]
) / num_macro_all

#-- calculate number of outer scenarios (third row of Table 2)
num_outer_scen_nested_plus <- num_outer
num_outer_scen_opt_nested <- mean(total_num_samp_all)
num_outer_scen_reg_nested <- mean(ceiling(colMeans(sim_budgets_all)))
num_outer_scen_std_nested <- mean(sim_budget_std_nested_all[1, ])
num_outer_scen <- c(
  num_outer_scen_nested_plus, num_outer_scen_opt_nested,
  num_outer_scen_std_nested, num_outer_scen_reg_nested
)

table2 <- rbind(budget, runtime, num_outer_scen)
colnames(table2) <- c("SNS_plus", "Optimal Design", "SNS", "Regression")
rownames(table2) <- c("Simulation Budget", "CPU time", "Num. Outer Scen.")
table2
table2_reduction <- table2[, 1] / table2
table2_reduction

#-- Calculate different error measures
# average MSE (AMSE in paper) for different methods (first row in Table 3)
amse_nested_plus <- mean((val_true_all - val_nested_plus_all)^2)
amse_opt_nested <- mean((val_true_all - val_opt_nested_all)^2)
amse_reg_nested <- mean((val_true_all - val_reg_nested_all)^2)
amse_std_nested <- mean(mapply(
  function(x, y) mean((x - y)^2),
  val_true_nested_all, val_std_nested_all
))
amse <- c(
  amse_nested_plus, amse_opt_nested,
  amse_std_nested, amse_reg_nested
)
amse # print the AMSEs to the console, first row in Table 3


erm_anova_calc <- function(risks_true, val_est, alpha, threshold) {
  # Calculate the bias, bias sqaures, variance, and MSE of the estimates
  # Input:
  #   risks_true: a vector of true risk measures
  #   val_est: a matrix or a list of estimates
  #   alpha: VaR level
  #   threshold: threshold for indicator, hockey stick, and square risk measures
  # Output:
  #   a data frame of bias, bias squares, variance, and MSE

  if (is.matrix(val_est)) {
    # split into list of column vectors for convenience
    val_est <- asplit(val_est, MARGIN = 2)
  }

  mean_est <- sapply(val_est, function(x) mean(x))
  VaR_est <- sapply(val_est, function(x) quantile(x, alpha))
  indicator_est <- sapply(val_est, function(x) mean(x > threshold))
  hockey_stick_est <- sapply(val_est, function(x) mean(pmax(x - threshold, 0)))
  square_est <- sapply(val_est, function(x) mean((pmax(x - threshold, 0))^2))
  risks_est <- cbind(
    mean_est, VaR_est,
    indicator_est, hockey_stick_est, square_est
  )

  bias <- colMeans(sweep(risks_est, MARGIN = 2, risks_true, "-"))
  bias2 <- bias^2
  var <- colMeans((sweep(risks_est, MARGIN = 2, colMeans(risks_est), "-"))^2)
  mse <- colMeans((sweep(risks_est, MARGIN = 2, risks_true, "-"))^2)

  return(data.frame(Bias = bias, Bias2 = bias2, Var = var, MSE = mse))
}

alpha <- 0.995 # VaR level
# threshold for indicator, hockey stick, and square risk measures
threshold <- 2800

# calculate true risk measures
val_true <- c(val_true_all) # convert to a vector for convenience
mean_true <- mean(val_true) # mean of the true values
VaR_true <- quantile(val_true, alpha) # VaR of the true values
indicator_true <- mean(val_true > threshold)
hockey_stick_true <- mean(pmax(val_true - threshold, 0))
square_true <- mean((pmax(val_true - threshold, 0))^2)
risks_true <- cbind(
  mean_true, VaR_true,
  indicator_true, hockey_stick_true, square_true
)
risks_true # print the true risk measures to the console

ANOVA_nested_plus <- erm_anova_calc(risks_true, val_nested_plus_all, alpha, threshold)
ANOVA_opt_nested <- erm_anova_calc(risks_true, val_opt_nested_all, alpha, threshold)
ANOVA_reg_nested <- erm_anova_calc(risks_true, val_reg_nested_all, alpha, threshold)
ANOVA_std_nested <- erm_anova_calc(risks_true, val_std_nested_all, alpha, threshold)

table3 <- rbind(
  amse,
  cbind(
    ANOVA_nested_plus[2:5, 4], ANOVA_opt_nested[2:5, 4],
    ANOVA_std_nested[2:5, 4], ANOVA_reg_nested[2:5, 4]
  )
)
rownames(table3) <- c(
  "AMSE", "MSE_Quantile", "MSE_Indicator",
  "MSE_HockeyStick", "MSE_Square"
)
table3

table3_reduction <- table3 / table3[, 1]
table3_reduction

#-- save the results
filename <- sprintf(
  "MultiAssetERM_%dStocks%dOuter%dMacro.RData",
  num_stocks, num_outer, num_macro_all
)
save(
  list = c(
    "num_outer", "num_inner", "num_macro", # simulation design (per seed)
    "num_stocks", "S0", "mu", "rf", # multivariate Black-Scholes parameters
    "vol", "Sigma", "T", "tau", "h",
    "K_euro", "type_euro", "position_euro", "num_euro", # European options
    "K_asian", "type_asian", "position_asian", "num_asian", # Asian options
    "K_barrier", "H_barrier", "U_barrier", "type_barrier", # Barrier options
    "position_barrier", "num_barrier",
    "num_macro_all", "s_tau_mac_all", "val_true_all", # combined results
    "val_nested_plus_all", "time_nested_plus_all", # nested plus simulation
    "val_opt_nested_all", "sim_budgets_all", # optimal LR design
    "total_num_samp_all", "time_opt_nested_opt_all",
    "time_opt_nested_sim_all", "time_opt_nested_tot_all",
    "val_reg_nested_all", "time_reg_nested_sim_all", # regression-based nested
    "time_reg_nested_fit_all", "time_reg_nested_tot_all",
    "s_tau_reg_nested_all",
    "val_true_nested_all", "val_std_nested_all", # standard nested simulation
    "time_std_nested_all", "s_tau_std_nested_all", "sim_budget_std_nested_all",
    "table2", "table2_reduction", "table3", "table3_reduction"
  ),
  file = filename
)
