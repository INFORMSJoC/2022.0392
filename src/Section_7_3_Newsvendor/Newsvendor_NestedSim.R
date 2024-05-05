################################################################################
# Newsvendor_NestedSim.R                                                       #
# Newsvendor example in Section 7.3 in the paper                               #
# Single-stage newsvendor problem with 10 products                             #
# Product demands X_ell, ell = 1,...,10 follow ind. Poisson distributions      #
# Total profit g(X) = sum (p_ell * min{X_ell, k_ell} - c_ell k_ell)            #
# This script produces Table 4 in Section 7.3                                  #
################################################################################

#-- Clear the workspace and load the necessary libraries
rm(list = ls()) # clear the global environment (clear all variables)
require(lpSolve) # library for solving the sim-budget minimizing LP
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

n_true <- 1e6 # number of true samples for credible interval calculation

# simulate the Poisson mean demands from posterior gamma distribution
lambda_cred_int_true <- matrix(0, n_true, n_products)
for (i in 1:n_products) {
  lambda_cred_int_true[, i] <- rgamma(n_true,
    rate = gamma_params[i, 1], shape = gamma_params[i, 2]
  )
}

# Conditional means for n_true outer params, for credible interval calcualtion
profit_cred_int_true <- apply(
  lambda_cred_int_true, 1,
  function(x) newsvendor_profit_oracle(x, k_ell, p_ell, c_ell)
)

#-- Nested simulation setting
n_outer <- 1e3 # number of outer scenarios
n_inner <- 1e3 # (desired) number of inner replications
n_macro <- 1e3 # number of macro replications

# placeholder for profits calculated via different methods
profit_oracle <- vector("list", length = n_macro)
profit_opt_nested <- vector("list", length = n_macro)
profit_reg_nested <- vector("list", length = n_macro)
profit_std_nested <- vector("list", length = n_macro)
profit_std_plus <- vector("list", length = n_macro)

sim_budget_opt_nested <- rep(0, n_macro) # sim. budget for optimal LR design

#-- main simulation experiment
for (mac in 1:n_macro) {
  print(sprintf("-- The %d macro rep.", mac))

  #-- simulate outer parameters
  lambda_outer <- matrix(0, n_outer, n_products)
  for (i in 1:n_products) {
    lambda_outer[, i] <- rgamma(n_outer,
      rate = gamma_params[i, 1], shape = gamma_params[i, 2]
    )
  }

  #-- oracle calculation
  profit_oracle[[mac]] <- apply(
    lambda_outer, 1,
    function(x) newsvendor_profit_oracle(x, k_ell, p_ell, c_ell)
  )

  #-- Standard nested simulation with MN simulation budget
  profit_std_plus[[mac]] <- apply(
    lambda_outer, 1,
    function(x) {
      mean(newsvendor_profit_simulator(x,
        k_ell, p_ell, c_ell,
        n_reps = n_inner
      )$profits)
    }
  )

  #-- Optimal LR design
  output_temp <- newsvendor_profit_opt_nested(
    lambda_outer,
    k_ell, p_ell, c_ell,
    n_inner
  )
  profit_opt_nested[[mac]] <- output_temp$profits
  sim_budget_opt_nested[mac] <- output_temp$sim_budget

  #-- Regression estimation
  n_outer_reg <- sim_budget_opt_nested[mac]
  n_inner_reg <- 1
  lambda_outer_reg <- matrix(0, n_outer_reg, n_products)
  for (i in 1:n_products) {
    lambda_outer_reg[, i] <- rgamma(n_outer_reg,
      rate = gamma_params[i, 1], shape = gamma_params[i, 2]
    )
  }
  profit_reg_nested[[mac]] <- newsvendor_profit_reg_nested(
    lambda_outer,
    k_ell, p_ell, c_ell,
    lambda_outer_reg, n_inner_reg
  )

  #-- Standard nested simulation with same budget as optimal LR design
  n_inner_std_nested <- ceiling(sim_budget_opt_nested[mac]^(1 / 3))
  n_outer_std_nested <- ceiling(sim_budget_opt_nested[mac]^(2 / 3))
  lambda_outer_std_nested <- matrix(0, n_outer_std_nested, n_products)
  for (i in 1:n_products) {
    lambda_outer_std_nested[, i] <- rgamma(n_outer_std_nested,
      rate = gamma_params[i, 1], shape = gamma_params[i, 2]
    )
  }
  profit_std_nested[[mac]] <- apply(
    lambda_outer_std_nested, 1,
    function(x) {
      mean(newsvendor_profit_simulator(x,
        k_ell, p_ell, c_ell,
        n_reps = n_inner_std_nested
      )$profits)
    }
  )
}

#-- postprocessing and plotting
pctile <- c(0.90, 0.95, 0.99) # percentiles for credible interval calculation
ci_pctile <- c((1 - pctile) / 2, 1 - (1 - pctile) / 2) # CI percentiles
n_pctile <- length(pctile) # number of percentiles of interest

ci_oracle <- sapply(
  profit_oracle,
  function(x) quantile(x, probs = ci_pctile)
)
ci_coverage_oracle <- matrix(0, nrow = n_pctile, ncol = n_macro)
ci_width_oracle <- matrix(0, nrow = n_pctile, ncol = n_macro)

ci_std_plus <- sapply(
  profit_std_plus,
  function(x) quantile(x, probs = ci_pctile)
)
ci_coverage_std_plus <- matrix(0, nrow = n_pctile, ncol = n_macro)
ci_width_std_plus <- matrix(0, nrow = n_pctile, ncol = n_macro)

ci_opt_nested <- sapply(
  profit_opt_nested,
  function(x) quantile(x, probs = ci_pctile)
)
ci_coverage_opt_nested <- matrix(0, nrow = n_pctile, ncol = n_macro)
ci_width_opt_nested <- matrix(0, nrow = n_pctile, ncol = n_macro)

ci_reg_nested <- sapply(
  profit_reg_nested,
  function(x) quantile(x, probs = ci_pctile)
)
ci_coverage_reg_nested <- matrix(0, nrow = n_pctile, ncol = n_macro)
ci_width_reg_nested <- matrix(0, n_pctile, ncol = n_macro)

ci_std_nested <- sapply(
  profit_std_nested,
  function(x) quantile(x, probs = ci_pctile)
)
ci_coverage_std_nested <- matrix(0, nrow = n_pctile, ncol = n_macro)
ci_width_std_nested <- matrix(0, n_pctile, ncol = n_macro)

for (i in 1:n_pctile) {
  print(sprintf("-- %d%%-tile CI coverage and width.", 100 * pctile[i]))

  # oracle CI coverage and width for the current percentile
  ci_current_pctile <- cbind(ci_oracle[i, ], ci_oracle[i + n_pctile, ])
  ci_coverage_oracle[i, ] <- apply(
    ci_current_pctile, 1,
    function(x) mean(between(profit_cred_int_true, x[1], x[2]))
  )
  ci_width_oracle[i, ] <- ci_current_pctile[, 2] - ci_current_pctile[, 1]

  # standard nested plus CI coverage and width for the current percentile
  ci_current_pctile <- cbind(ci_std_plus[i, ], ci_std_plus[i + n_pctile, ])
  ci_coverage_std_plus[i, ] <- apply(
    ci_current_pctile, 1,
    function(x) mean(between(profit_cred_int_true, x[1], x[2]))
  )
  ci_width_std_plus[i, ] <- ci_current_pctile[, 2] - ci_current_pctile[, 1]

  # optimal nested CI coverage and width for the current percentile
  ci_current_pctile <- cbind(ci_opt_nested[i, ], ci_opt_nested[i + n_pctile, ])
  ci_coverage_opt_nested[i, ] <- apply(
    ci_current_pctile, 1,
    function(x) mean(between(profit_cred_int_true, x[1], x[2]))
  )
  ci_width_opt_nested[i, ] <- ci_current_pctile[, 2] - ci_current_pctile[, 1]

  # regression nested CI coverage and width for the current percentile
  ci_current_pctile <- cbind(ci_reg_nested[i, ], ci_reg_nested[i + n_pctile, ])
  ci_coverage_reg_nested[i, ] <- apply(
    ci_current_pctile, 1,
    function(x) mean(between(profit_cred_int_true, x[1], x[2]))
  )
  ci_width_reg_nested[i, ] <- ci_current_pctile[, 2] - ci_current_pctile[, 1]

  # standard nested CI coverage and width for the current percentile
  ci_current_pctile <- cbind(ci_std_nested[i, ], ci_std_nested[i + n_pctile, ])
  ci_coverage_std_nested[i, ] <- apply(
    ci_current_pctile, 1,
    function(x) mean(between(profit_cred_int_true, x[1], x[2]))
  )
  ci_width_std_nested[i, ] <- ci_current_pctile[, 2] - ci_current_pctile[, 1]
}

results_oracle <- data.frame(
  Method = "oracle",
  alpha = c(0.9, 0.95, 0.99),
  coverage_mean = rowMeans(ci_coverage_oracle),
  coverage_stderr = sqrt((rowMeans(ci_coverage_oracle^2) -
    rowMeans(ci_coverage_oracle)^2) / n_macro),
  width_mean = rowMeans(ci_width_oracle),
  width_stderr = sqrt((rowMeans(ci_width_oracle^2) -
    rowMeans(ci_width_oracle)^2)) / n_macro
)

results_std_plus <- data.frame(
  Method = "std_plus",
  alpha = c(0.9, 0.95, 0.99),
  coverage_mean = rowMeans(ci_coverage_std_plus),
  coverage_stderr = sqrt((rowMeans(ci_coverage_std_plus^2) -
    rowMeans(ci_coverage_std_plus)^2) / n_macro),
  width_mean = rowMeans(ci_width_std_plus),
  width_stderr = sqrt((rowMeans(ci_width_std_plus^2) -
    rowMeans(ci_width_std_plus)^2)) / n_macro
)

results_opt_nested <- data.frame(
  Method = "opt_nested",
  alpha = c(0.9, 0.95, 0.99),
  coverage_mean = rowMeans(ci_coverage_opt_nested),
  coverage_stderr = sqrt((rowMeans(ci_coverage_opt_nested^2) -
    rowMeans(ci_coverage_opt_nested)^2) / n_macro),
  width_mean = rowMeans(ci_width_opt_nested),
  width_stderr = sqrt((rowMeans(ci_width_opt_nested^2) -
    rowMeans(ci_width_opt_nested)^2)) / n_macro
)

results_reg_nested <- data.frame(
  Method = "reg_nested",
  alpha = c(0.9, 0.95, 0.99),
  coverage_mean = rowMeans(ci_coverage_reg_nested),
  coverage_stderr = sqrt((rowMeans(ci_coverage_reg_nested^2) -
    rowMeans(ci_coverage_reg_nested)^2) / n_macro),
  width_mean = rowMeans(ci_width_reg_nested),
  width_stderr = sqrt((rowMeans(ci_width_reg_nested^2) -
    rowMeans(ci_width_reg_nested)^2)) / n_macro
)

results_std_nested <- data.frame(
  Method = "std_nested",
  alpha = c(0.9, 0.95, 0.99),
  coverage_mean = rowMeans(ci_coverage_std_nested),
  coverage_stderr = sqrt((rowMeans(ci_coverage_std_nested^2) -
    rowMeans(ci_coverage_std_nested)^2) / n_macro),
  width_mean = rowMeans(ci_width_std_nested),
  width_stderr = sqrt((rowMeans(ci_width_std_nested^2) -
    rowMeans(ci_width_std_nested)^2)) / n_macro
)

results <- rbind(
  results_oracle, results_opt_nested,
  results_std_plus, results_std_nested, results_reg_nested
)

#-- save the results
filename <- sprintf(
  "Newsvendor_%dInner%dOuter%dMacro.RData",
  n_inner, n_outer, n_macro
)
save.image(filename)
write.csv(results, "Newsvendor_NestedSim.csv")

print("Table 4 in Section 7.3 is based on the table below.")
results

print(sprintf(
  "The average simulation budget for the optimal design is %4.2f.",
  mean(sim_budget_opt_nested)
))

print(sprintf(
  "The optimal design's budget is %4.2f times smaller than the standard
  nested design.",
  (n_outer * n_inner) / mean(sim_budget_opt_nested)
))
