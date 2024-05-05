################################################################################
# Newsvendor_SimBudgetGrowth.R                                                 #
# Newsvedor example in Section 7.3 in the paper                                #
# Single-stage newsvendor problem with 10 products                             #
# Product demands X_ell, ell = 1,...,10 follow ind. Poisson distributions      #
# Total profit g(X) = sum (p_ell * min{X_ell, k_ell} - c_ell k_ell)            #
# This script produces Figure 3a in Section 7.3                                #
################################################################################

# ===== Preparation of workspace
rm(list = ls()) # clear the global environment (clear all variables)
library(lpSolve) # for solving LP (e.g., Eq. (11) in the paper)
library(ggplot2) # for plotting
source("Newsvendor_Helpers.R")
set.seed(2020) # fix random seed for reproducibility

#-- True parameter settings
n_products <- 10 # number of products, i.e., the problem dimension

ell <- 1:n_products # product index
mean_demands <- 5 + ell # avearge demand for different products
k_ell <- 9 + ell # stocking policy for each product
p_ell <- 7 + 3 * ell # unit sale price for each product
c_ell <- 2 # unit cost for each product

# sample size for each product to estimate the demand distribution parameters
sample_sizes <- 50 + 5 * ell

# calculates the posterior gamma distribution parameters for the Poisson rate
gamma_prior_params <- c(0.001, 0.001) # prior gamma parameters
gamma_params <- matrix(0, n_products, 2) # posterior gamma parameters
for (i in 1:n_products) {
  # Poisson samples from true distribution
  poisson_samples <- rpois(sample_sizes[i], mean_demands[i])
  # update the posterior gamma distribution parameters
  gamma_params[i, ] <- gamma_prior_params +
    c(sample_sizes[i], sum(poisson_samples))
}

#-- Nested simulation setting
n_outer <- 1e3 # number of outer scenarios
n_inner <- 1e3 # (desired) number of inner replications
n_macro <- 1e3 # number of macro replications

#-- simulation one fixed set of outer scenarios
lambda_outer <- matrix(0, n_outer, n_products)
for (i in 1:n_products) {
  lambda_outer[, i] <- rgamma(n_outer,
    rate = gamma_params[i, 1], shape = gamma_params[i, 2]
  )
}

# placeholder for profits calculated via different methods
profit_opt_nested <- matrix(0, n_outer, n_macro)

#-- main simulation experiment
for (mac in 1:n_macro) {
  print(sprintf("-- The %d macro rep.", mac))

  #-- Optimal LR design
  output_temp <- newsvendor_profit_opt_nested(
    lambda_outer,
    k_ell, p_ell, c_ell,
    n_inner
  )
  profit_opt_nested[, mac] <- output_temp$profits
}

# true variances at outer scenarios
var_oracle <- apply(
  lambda_outer, 1,
  function(x) newsvendor_var_oracle(x, k_ell, p_ell)
) / n_inner
var_opt_nested <- apply(profit_opt_nested, MARGIN = 1, var)

var_ratio <- var_oracle / var_opt_nested

filename <- sprintf(
  "Newsvendor_VarRatio_%dInner%dOuter%dMacro.RData",
  n_inner, n_outer, n_macro
)
save.image(filename)

df <- data.frame(VarRatio = var_ratio)

ggplot(df, aes(x = VarRatio)) +
  geom_histogram(bins = 50, color = "black", fill = "red", alpha = 0.4) +
  geom_vline(aes(xintercept = mean(VarRatio)),
    color = "blue", linetype = "dashed", linewidth = 1
  ) +
  theme(text = element_text(size = 16, face = "bold")) +
  ylab("Count") +
  xlab("Variance Ratio")

ggsave("VarRatio.pdf",
  width = 7, height = 5,
  units = "in", dpi = 600
) # save the plot as a pdf
