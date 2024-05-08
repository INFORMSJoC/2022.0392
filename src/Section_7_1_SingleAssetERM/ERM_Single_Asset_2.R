################################################################################
# ERM_Single_Asset2.R                                                          #
# Single-Asset ERM example in Section 7.1 in the paper                         #
# Straddle option portfolio that consists of a long call and a long put        #
# This script produces Table 1 in Section 7.1 in the paper                     #
################################################################################

#-- Clear the workspace and load the necessary libraries
rm(list = ls())
library(lpSolve) # for solving LP (e.g., Eq. (11) in the paper)
library(latex2exp) # for TeX expressions
library(ggplot2) # for plotting
library(dplyr) # for data manipulations (for plotting)
library(tidyverse) # for data manipulations (for plotting)
source("NestedSim_ERM_Helpers.R") # helper functions

#-- Black-Scholes model parameters
s0 <- 100 # initial asset price
mu <- 0.05 # annual expected return (real world drift)
rf <- 0.02 # annual risk-free rate (risk neutral drift)
div <- 0 # annual diviDend rate
vol <- 0.3 # annual volatility
str_k <- 110 # strike price (for both options)
big_t <- 2 # time-to-maturity
tau <- 1 / 4 # risk horizon

#-- Derived parameters, mean and standard deviation of the normal distribution
#   for the outer and inner simulation
meanlog_outer <- (mu - div - vol^2 / 2) * tau # outer sim. normal mean
sdlog_outer <- vol * sqrt(tau) # outer sim. normal sd
meanlog_inner <- (rf - div - vol^2 / 2) * (big_t - tau) # inner sim. normal mean
sdlog_inner <- vol * sqrt(big_t - tau) # inner sim. normal sd

#-- Calculate "true" 99%-quantile using 10^8 conditional means
n_true <- 1e8 # number of outer scenarios in "true" distribution
qts_true <- (1:n_true) / (n_true + 1)
s_tau_true <- qlnorm(qts_true, log(s0) + meanlog_outer, sdlog_outer)
mu_s_tau_true <- sapply(
  s_tau_true,
  function(x) bs_price_straddle(x, str_k, div, rf, vol, (big_t - tau))
) # true conditional mean
rm(list = c("qts_true", "s_tau_true")) # remove intermediate variables

#-- Calculate different tail risk measures
alpha <- 0.99
VaR_alpha <- quantile(mu_s_tau_true, alpha) # 99%-quantile
print(sprintf(
  "The 99%%-quantile based on %0.2g scenario is %4.2f.",
  n_true, VaR_alpha
))

threshold <- 49
indicator_true <- mean(mu_s_tau_true > threshold) # indicator function loss
hockeystick_true <- mean(pmax(mu_s_tau_true - threshold, 0)) # hockey-stick loss
square_true <- mean((pmax(mu_s_tau_true - threshold, 0))^2) # square loss
risks_true <- cbind(VaR_alpha, indicator_true, hockeystick_true, square_true)

#-- Run nested simulation with different number of outer scenarios
n_outer_set <- 2^(7:12) # a set of numbers of outer scenarios
n_cases <- length(n_outer_set) # number of cases

#-- Placeholders of results
#-- List of outer scenarios in different cases
s_tau_list <- vector(mode = "list", length = n_cases)
sim_budget <- rep(0, n_cases) # total simulation budgets in different cases

#-- list of estimated conditional means for different methods
val_opt_nested_list <- vector(mode = "list", length = n_cases)
val_reg_nested_list <- vector(mode = "list", length = n_cases)
val_std_nested_list <- vector(mode = "list", length = n_cases)
val_std_plus_list <- vector(mode = "list", length = n_cases)

n_macro <- 1e4 # nubmer of macro reps
set.seed(2024)
for (case in 1:n_cases) { # loop over different n_outer cases

  n_outer <- n_outer_set[case] # number of outer scenarios
  n_inner <- n_outer_set[case] # number of inner reps, set the same as n_outer

  #-- Outer simulation (real world projection from time 0 to tau)
  #-- deterministic outer scenarios
  qts <- (1:n_outer) / (n_outer + 1) # equally-spaced quantiles
  s_tau <- qlnorm(qts, log(s0) + meanlog_outer, sdlog_outer)
  # mean(s_tau) - s0*exp(mu*tau) # unit test, result should be close to 0
  #-- stochastic outer scenarios (not used in the paper)
  # s_tau<-s0*exp(meanlog_Outer + rnorm(n_outer)*sdlog_Outer)
  # s_tau<-rlnorm(n_outer,log(s0)+meanlog_Outer,sdlog_Outer)
  s_tau_list[[case]] <- s_tau

  #-- Optimal deisgn (optimal sampling and pooling decisions)
  #-- LP optimization (Eq. (11) in the paper)
  #-- matrix of the inverses of second moments of likelihood ratios
  psi_mat <- psi_mat_lnorm(log(s_tau) + meanlog_inner, sdlog_inner)
  opt_results <- lp(
    direction = "min",
    objective.in = rep(1, n_outer), # obj. coeff, all ones
    const.mat = psi_mat, # constraint matrix
    const.dir = rep(">=", n_outer), # constraint direction
    const.rhs = rep(1, n_outer), # RHS, n_inner
    all.int = FALSE
  ) # round up after solving the LP

  n_j <- opt_results$solution # (non-integral) optimal number of inner reps
  ind_nz <- which(n_j > 1e-8) # scenario indices of non-zero inner reps
  n_nz <- length(ind_nz) # number of non-zero scenarios
  n_j_nz <- n_j[ind_nz] # optimal number of inner reps in non-zero scenarios
  n_j_opt_nested <- ceiling(n_j_nz * n_inner) # ceiling of optimal inner reps
  sim_budget_opt_nested <- sum(n_j_opt_nested) # total simulation budget
  sim_budget[case] <- sim_budget_opt_nested # save the total simulation budget

  psi_mat_nz <- psi_mat[, ind_nz] # submatrix of psi_mat for non-zero scenarios
  gamma_nz <- sweep(psi_mat_nz, MARGIN = 2, n_j_opt_nested, "*")
  gamma_nz <- gamma_nz / rowSums(gamma_nz) # optimal pooling weights

  #-- print out the optimal design
  print(sprintf("==== Optimal design"))
  print(sprintf(
    "---- %d non-zero scenarios: theta = %4.2f",
    n_nz, s_tau[ind_nz]
  ))
  print(sprintf("---- Simulation budget %d.", sim_budget_opt_nested))

  #-- Optimal LR design simulation
  # placeholder for estimated conditional means
  val_opt_nested <- matrix(0, nrow = n_outer, ncol = n_macro)
  # loop through macro replications
  for (mac in 1:n_macro) {
    print(sprintf("Optimal nested sim, macro rep %d", mac))

    # optimal nested simulation design with likelihood ratio
    val_opt_nested[, mac] <- straddle_opt_nested(
      s_tau, str_k, div, rf, vol, big_t, tau,
      s_tau[ind_nz], n_j_opt_nested, gamma_nz
    )
  }
  val_opt_nested_list[[case]] <- val_opt_nested

  #-- Standard Nested Simulation with Gamma^(1/3) inner reps and Gamma^(2/3)
  #   outer scenarios. This method is referred to as SNS in the paper
  #-- number of inner and outer replications
  n_inner_std_nested <- ceiling(sim_budget_opt_nested^(1 / 3))
  n_outer_std_nested <- ceiling(sim_budget_opt_nested^(2 / 3))
  #-- equally-spaced quantiles as outer scenarios
  qts_std_nested <- (1:n_outer_std_nested) / (n_outer_std_nested + 1)
  s_tau_std_nested <- s0 * exp(meanlog_outer +
    qnorm(qts_std_nested) * sdlog_outer)

  # placeholder for estimated conditional means
  val_std_nested <- matrix(0, nrow = n_outer_std_nested, ncol = n_macro)
  for (mac in 1:n_macro) {
    print(sprintf("Stadard nested sim with optimal budget, macro rep %d", mac))

    # standard nested simulation with similar budget and the optimal design
    val_std_nested[, mac] <- straddle_std_nested(
      s_tau_std_nested, str_k, div, rf, vol,
      big_t, tau, n_inner_std_nested
    )
  }
  val_std_nested_list[[case]] <- val_std_nested

  #-- Standard Nested Simulation with NInner inner reps per outer scenario.
  #-- This method is referred to as SNS+ in the paper
  #-- placeholder for estimated conditional means and runtime
  val_std_plus <- matrix(0, nrow = n_outer, ncol = n_macro)
  # loop through macro replications
  for (mac in 1:n_macro) {
    print(sprintf("Standard nested sim with MN budget, macro rep %d", mac))

    # standard nested simulation
    val_std_plus[, mac] <- straddle_std_nested(
      s_tau, str_k, div, rf, vol,
      big_t, tau, n_inner
    )
  }
  val_std_plus_list[[case]] <- val_std_plus

  #-- Regression-based Nested Simulation
  n_inner_reg <- 1 # number of inner replications at each design point
  # number of design points
  n_outer_reg <- ceiling(sim_budget_opt_nested / n_inner_reg)
  # equally-spaced quantiles of s_tau as design points
  qts_reg <- (1:n_outer_reg) / (n_outer_reg + 1)
  s_tau_reg <- qlnorm(qts_reg, log(s0) + meanlog_outer, sdlog_outer)

  # placeholder for estimated conditional means and runtime
  val_reg_nested <- matrix(0, nrow = n_outer, ncol = n_macro)
  # loop through macro replications
  for (mac in 1:n_macro) {
    print(sprintf("Regression nested sim, macro rep %d", mac))

    # regression-based nested simulation
    val_reg_nested[, mac] <- straddle_reg_nested(
      s_tau, str_k, div, rf, vol, big_t, tau,
      s0, s_tau_reg, n_inner_reg
    )
  }
  val_reg_nested_list[[case]] <- val_reg_nested
}

df_tbl1 <- data.frame(
  ErrorType = character(), RiskType = character(),
  Val = numeric(), NINOut = numeric(), SimBudget = numeric()
)

for (case in 1:n_cases) { # loop over different n_outer cases
  n_outer <- n_outer_set[case] # number of outer scenarios

  val_opt_nested <- val_opt_nested_list[[case]]
  risks_opt_nested <- erm_risk_measures(val_opt_nested, alpha, threshold)
  colnames(risks_opt_nested) <- paste(
    colnames(risks_opt_nested), "OptNested",
    sep = "_"
  )
  errors_opt_nested <- erm_mse_func(risks_true, risks_opt_nested)
  errors_opt_nested <- errors_opt_nested %>%
    tibble::rownames_to_column("ErrorType") %>%
    pivot_longer(-ErrorType, names_to = "RiskType", values_to = "Val")

  val_std_nested <- val_std_nested_list[[case]]
  risks_std_nested <- erm_risk_measures(val_std_nested, alpha, threshold)
  colnames(risks_std_nested) <- paste(
    colnames(risks_std_nested), "StdNested",
    sep = "_"
  )
  errors_std_nested <- erm_mse_func(risks_true, risks_std_nested)
  errors_std_nested <- errors_std_nested %>%
    tibble::rownames_to_column("ErrorType") %>%
    pivot_longer(-ErrorType, names_to = "RiskType", values_to = "Val")

  val_std_plus <- val_std_plus_list[[case]]
  risks_std_plus <- erm_risk_measures(val_std_plus, alpha, threshold)
  colnames(risks_std_plus) <- paste(
    colnames(risks_std_plus), "StdPlus",
    sep = "_"
  )
  errors_std_plus <- erm_mse_func(risks_true, risks_std_plus)
  errors_std_plus <- errors_std_plus %>%
    tibble::rownames_to_column("ErrorType") %>%
    pivot_longer(-ErrorType, names_to = "RiskType", values_to = "Val")

  val_reg_nested <- val_reg_nested_list[[case]]
  risks_reg_nested <- erm_risk_measures(val_reg_nested, alpha, threshold)
  colnames(risks_reg_nested) <- paste(
    colnames(risks_reg_nested), "RegNested",
    sep = "_"
  )
  errors_reg_nested <- erm_mse_func(risks_true, risks_reg_nested)
  errors_reg_nested <- errors_reg_nested %>%
    tibble::rownames_to_column("ErrorType") %>%
    pivot_longer(-ErrorType, names_to = "RiskType", values_to = "Val")

  new_row <- rbind(
    errors_opt_nested, errors_std_nested,
    errors_std_plus, errors_reg_nested
  )
  new_row$n_inner_outer <- n_outer
  new_row$SimBudget <- sim_budget[case]

  df_tbl1 <- rbind(df_tbl1, new_row)
}

df_tbl1 %>%
  filter(ErrorType == "MSE", grepl("VaR", RiskType)) %>%
  spread(RiskType, Val) %>%
  select(
    n_inner_outer, VaR_OptNested, VaR_StdNested,
    VaR_StdPlus, VaR_RegNested
  )
df_tbl1 %>%
  filter(ErrorType == "MSE", grepl("indicator", RiskType)) %>%
  spread(RiskType, Val) %>%
  select(
    n_inner_outer, indicator_OptNested, indicator_StdNested,
    indicator_StdPlus, indicator_RegNested
  )
df_tbl1 %>%
  filter(ErrorType == "MSE", grepl("hockeystick", RiskType)) %>%
  spread(RiskType, Val) %>%
  select(
    n_inner_outer, hockeystick_OptNested, hockeystick_StdNested,
    hockeystick_StdPlus, hockeystick_RegNested
  )
df_tbl1 %>%
  filter(ErrorType == "MSE", grepl("square", RiskType)) %>%
  spread(RiskType, Val) %>%
  select(
    n_inner_outer, square_OptNested,
    square_StdNested, square_StdPlus, square_RegNested
  )

filename <- sprintf("SingleAssetERM_2_%dCases%dMacros.RData", n_cases, n_macro)
save.image(filename)
