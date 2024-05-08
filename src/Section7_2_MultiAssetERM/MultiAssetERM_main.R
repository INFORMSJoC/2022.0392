################################################################################
# MultiAssetERM_main.R                                                         #
# Multi-Asset ERM example in Section 7.2 in the paper                          #
# An option portfolio with 300 options written on 10 underlying assets         #
# For each underlying asset,                                                   #
# -- 10 European options: 5 puts and 5 calls                                   #
# -- 10 geometric average Asian options: 5 puts and 5 calls                    #
# -- 10 down-and-out put Barrier options                                       #
# This script produces Table 2 & Table 3 in Section 7.3 of the paper.          #
################################################################################

#-- Clear the workspace and load the necessary libraries
rm(list = ls())
library(lpSolveAPI) # for solving LP (e.g., Eq. (11) in the paper)
library(mvtnorm) # bivariate normal pdf for closed-form barrier option price
library(randcorr) # simulate random correlation matrix
source("MultiAssetERM_helpers.R") # helper functions

#-- set a working directory for saving results
setwd("/u/mbfeng/BenResearch/IJOC_Github/ERM_Multi_Assets")
#-- set seeds, run simulations, save file
args <- commandArgs() # get command line arguments
seed_number <- as.numeric(tail(args, n = 1)) # use above argument to set  seed
set.seed(seed_number) # set seed for reproducibility

#-- experiment design parameters
num_outer <- 1e3 # number of outer scenarios
num_inner <- 1e3 # number of inner reps per scenario (N0 for ESS calculation)
num_macro <- 40 # number of macro reps per seed

#-- Black-Scholes model parameters and option parameters
num_stocks <- 10 # number of stocks, problem dimension
S0 <- rep(100, num_stocks) # initial asset price
mu <- seq(
  from = 0.05, to = 0.14,
  length = num_stocks
) # annual expected return (real world drift)
rf <- 0.02 # annual risk-free rate (risk neutral drift)
vol <- seq(from = 0.3, to = 0.48, length = num_stocks) # annual volatilities
# random correlation matrix
if (num_stocks == 1) {
  Sigma <- vol^2
} else {
  corr <- randcorr(num_stocks)
  Sigma <- diag(vol) %*% corr %*% diag(vol) # covariance matrix
}
T <- 2 # time-to-maturity (at time zero)
tau <- 4 / 52 # risk horizon
h <- 1 / 52 # time step size

#-- Derived parameters
T_inner <- T - tau # time to maturity (at time tau)
num_steps_outer <- round(tau / h) # number of steps in outer loop
num_steps_inner <- round(T_inner / h) # number of steps in inner loop
num_steps_total <- num_steps_outer + num_steps_inner

meanlog_outer <- (mu - vol^2 / 2) * h # per-step normal mean in outer loops
meanlog_inner <- (rf - vol^2 / 2) * h # per-step normal mean in inner loops
sdlog_outer <- vol * sqrt(h) # per-step normal sd in outer loops
sdlog_inner <- vol * sqrt(h) # per-step normal sd in inner loops

#-- European option parameters
K_euro <- seq(from = 85, to = 130, by = 5) # European option strikes
num_euro <- length(K_euro) # number of European options
type_euro <- c(
  rep("p", num_euro / 2),
  rep("c", num_euro / 2)
) # European option types
position_euro <- rep(1, num_euro) # positions, long (+) or short (-)

#-- Asian option parameters
K_asian <- seq(from = 85, to = 130, by = 5) # Asian option strikes
num_asian <- length(K_asian) # number of Asian options
type_asian <- c(
  rep("p", num_asian / 2),
  rep("c", num_asian / 2)
) # Asian option types
position_asian <- rep(1, num_asian) # positions, long (+) or short (-)

#-- Barrier option parameters
K_barrier <- seq(from = 85, to = 130, by = 5) # Barrier option strikes
H_barrier <- seq(from = 75, to = 120, by = 5) # lower barriers
num_barrier <- length(K_barrier) # number of barrier options
U_barrier <- rep(Inf, num_barrier) # upper barriers
type_barrier <- rep("pdo", num_barrier) # pdo: down-and-out put barrier option
position_barrier <- rep(1, num_barrier) # positions, long (+) or short (-)

#-- outer simulation (real world projection from time 0 to tau)
s_tau_mac <- array(0, dim = c(num_outer, num_stocks, num_macro))
for (m in 1:num_macro) {
  print(sprintf(
    "Outer Sim: %d macro replications of %d outer scenarios",
    num_macro, num_outer
  ))

  # simulate outer scenarios
  z_tau <- rmvnorm(num_outer,
    mean = meanlog_outer * num_steps_outer,
    sigma = Sigma * tau
  ) # num_outer-by-num_stocks return matrix

  # num_stocks-by-num_outer matrix
  s_tau_mac[, , m] <- sweep(x = exp(z_tau), MARGIN = 2, STATS = S0, FUN = "*")
}

#--- closed-form portfolio value for MSE calculation
val_true <- matrix(0, nrow = num_outer, ncol = num_macro)
for (m in 1:num_macro) {
  time_start <- Sys.time()

  for (i in 1:num_outer) {
    portfolio_value <- 0 # placeholder for portfolio value

    for (s in 1:num_stocks) {
      vol_s <- vol[s] # volatility of s-th stock
      # s-th time-tau price in i-th scen (m-th macro)
      s_tau_i <- s_tau_mac[i, s, m]

      for (k in 1:num_euro) {
        # true price of k-th European option
        price_k <- EuroOptionPrice(
          s_tau_i, T_inner, vol_s, rf,
          K_euro[k], type_euro[k], position_euro[k]
        )
        portfolio_value <- portfolio_value + price_k
      }

      for (k in 1:num_asian) {
        # true price of k-th Asian option
        price_k <- AsianOptionPrice(
          s_tau_i, K_asian[k], rf, vol_s,
          T_inner, num_steps_inner, type_asian[k], position_asian[k]
        )
        portfolio_value <- portfolio_value + price_k
      }

      for (k in 1:num_barrier) {
        # true price of k-th barrier option
        price_k <- BarrierOptionPrice_DOP(s_tau_i, K_barrier[k],
          H_barrier[k], rf,
          q = 0, vol_s, h, T_inner, position_barrier[k]
        )
        portfolio_value <- portfolio_value + price_k
      }
    }

    val_true[i, m] <- portfolio_value
  }
  time_end <- Sys.time()
  time_diff <- as.numeric(difftime(time_end, time_start, units = c("secs")))
  print(sprintf(
    "True Value: This is the %d-th macro, it took %8.2f seconds.",
    m, time_diff
  ))
}

#-- Standard Nested Simulation with n_inner inner reps per outer scenario.
#-- This method is referred to as SNS+ in the paper
# placeholder for estimated conditional means and runtime
val_nested_plus <- matrix(0, nrow = num_outer, ncol = num_macro)
time_nested_plus <- rep(0, 2)
for (m in 1:num_macro) {
  time_start <- Sys.time()

  for (i in 1:num_outer) {
    portfolio_value <- 0 # placeholder for portfolio value

    for (s in 1:num_stocks) {
      vol_s <- vol[s]
      s_tau_i <- s_tau_mac[i, s, m]

      # growth factor in inner simulation
      z_inner_full <- rnorm(num_steps_inner * num_inner,
        mean = meanlog_inner[s], sd = sdlog_inner[s]
      )
      z_inner_full <- matrix(z_inner_full, nrow = num_steps_inner)
      z_inner_full <- apply(z_inner_full, 2, cumsum)

      s_inner_full <- s_tau_i * exp(z_inner_full) # full inner price paths

      for (k in 1:num_euro) {
        # simulated price of k-th European option
        price_k <- EuroOptionSim(
          S_full_paths = s_inner_full,
          T = T_inner, r = rf, K = K_euro[k], option_type = type_euro[k],
          position = position_euro[k]
        )
        portfolio_value <- portfolio_value + mean(price_k)
      }

      for (k in 1:num_asian) {
        # simulated price of k-th Asian option
        price_k <- GeoAsianOptionSim(
          S_full_paths = s_inner_full,
          T = T_inner, r = rf, K = K_asian[k], option_type = type_asian[k],
          position = position_asian[k]
        )
        portfolio_value <- portfolio_value + mean(price_k)
      }

      for (k in 1:num_barrier) {
        # simulated price of k-th barrier option
        price_k <- BarrierOptionSim(
          S_full_paths = s_inner_full,
          T = T_inner, r = rf, K = K_barrier[k], vol_s, h,
          H = H_barrier[k], U = U_barrier[k],
          option_type = type_barrier[k], position = position_barrier[k]
        )
        portfolio_value <- portfolio_value + mean(price_k)
      }
    }

    val_nested_plus[i, m] <- portfolio_value
  }
  time_end <- Sys.time()
  time_diff <- as.numeric(difftime(time_end, time_start, units = c("secs")))
  time_nested_plus[1] <- time_nested_plus[1] + time_diff
  time_nested_plus[2] <- time_nested_plus[2] + time_diff^2
  print(sprintf(
    "Standard Nested Plus: This is the %d-th macro, it took %8.2f seconds.",
    m, time_diff
  ))
}

#-- Optimal LR design simulation
# placeholders for estimated conditional means and runtime
val_opt_nested <- matrix(0, nrow = num_outer, ncol = num_macro)
sim_budgets <- matrix(0, nrow = num_stocks, ncol = num_macro)
total_num_samp <- matrix(0, nrow = num_stocks, ncol = num_macro)
time_opt_nested_opt <- rep(0, 2)
time_opt_nested_sim <- rep(0, 2)
time_opt_nested_tot <- rep(0, 2)

for (m in 1:num_macro) {
  time_diff_opt <- 0 # time for optimization
  time_diff_sim <- 0 # time for simulation

  time_start_mac <- Sys.time() # start time for the macro rep

  for (s in 1:num_stocks) {
    s_tau <- s_tau_mac[, s, m]

    time_start_opt <- Sys.time() # start time for optimization, including setup

    #-- LP optimization (Eq. (11) in the paper)
    # matrix of the reciprocals of second moments of likelihood ratios
    psi_mat <- psi_mat_lnorm(log(s_tau) + meanlog_inner[s], sdlog_inner[s])
    # round off small values to improve numerical stability
    psi_mat[psi_mat < 1e-6] <- 0

    lprec <- make.lp(num_outer, num_outer) # create LP object
    set.objfn(lprec, rep(1, num_outer)) # minimize total simulation budget

    # lpSolveAPI requires constraints to be added one by one
    for (i in 1:num_outer) {
      add.constraint(lprec, psi_mat[i, ], ">=", 1)
    }

    # set bounds for each variablem, i.e., 0<= x_i <= 1
    set.bounds(lprec,
      lower = rep(0, num_outer), upper = rep(1, num_outer),
      columns = 1:num_outer
    )

    # set control parameters for lpSolveAPI to improve speed and stability
    lp.control(lprec,
      pivoting = c("devex"),
      simplextype = c("dual"),
      scaling = c("curtisreid")
    )

    # solve the LP and catch errors if any
    status <- solve(lprec)
    if (status != 0) {
      print(sprintf("LP solve is unsuccessful: status code is %d.", status))
      print(sprintf("---- Objective value is %10.2f", get.objective(lprec)))
      quit()
    }

    n_j <- get.variables(lprec) # (non-integral) optimal number of inner reps
    ind_nz <- which(n_j > 0) # scenario indices with non-zero inner reps
    n_nz <- length(ind_nz) # number of non-zero scenarios
    s_tau_samp <- s_tau[ind_nz] # non-zero scenarios

    n_j_nz <- n_j[ind_nz] # optimal number of inner reps in non-zero scenarios
    # number of inner reps at each scenario after rounding up, minimum 10 reps
    n_j_opt_nested <- pmax(ceiling(n_j_nz * num_inner), 10)

    # calculate optimal pooling weights
    psi_mat_nz <- psi_mat[, ind_nz] # submatrix of psi_mat for non-zero scen
    gamma_nz <- sweep(psi_mat_nz, MARGIN = 2, n_j_opt_nested, "*")
    gamma_nz <- gamma_nz / rowSums(gamma_nz)

    time_end_opt <- Sys.time() # end time for optimization
    # add current optimization time to the total optimization time
    time_diff_opt <- time_diff_opt +
      as.numeric(difftime(time_end_opt, time_start_opt, units = c("secs")))

    sim_budget_opt_nested <- sum(n_j_opt_nested) # total simulation budget
    sim_budgets[s, m] <- sim_budget_opt_nested
    total_num_samp[s, m] <- n_nz

    print(sprintf("Optimal Nested: This is %d-th macro, %d-th stock.", m, s))
    print(sprintf(
      "-- There are %d non-zero scenarios, total sim budget is %d.",
      n_nz, sim_budget_opt_nested
    ))

    # Simulate and reuse
    time_start_sim <- Sys.time() # start time for simulation

    storage_s_tauplus1 <- vector("list", n_nz)
    storage_samp_dens <- vector("list", n_nz)
    storage_port_val <- vector("list", n_nz)
    storage_log_s_tauplus1 <- vector("list", n_nz)
    storage_mu_samp <- vector("list", n_nz)

    for (i in 1:n_nz) { # only simulate at non-zero scenarios
      s_tau_samp_i <- s_tau_samp[i] # i-th non-zero scenario
      num_samp_i <- n_j_opt_nested[i] # number of inner reps at i-th scenario

      # growth factor in inner simulation
      z_inner_full <- rnorm(num_steps_inner * num_samp_i,
        mean = meanlog_inner[s], sd = sdlog_inner[s]
      )
      z_inner_full <- matrix(z_inner_full, nrow = num_steps_inner)
      z_inner_full <- apply(z_inner_full, 2, cumsum)

      s_inner_full <- s_tau_samp_i * exp(z_inner_full) # full inner price paths

      portfolio_value <- rep(0, num_samp_i) # placeholder for portfolio value
      for (k in 1:num_euro) {
        # simulated price of k-th European option
        price_k <- EuroOptionSim(
          S_full_paths = s_inner_full,
          T = T_inner, r = rf, K = K_euro[k], option_type = type_euro[k],
          position = position_euro[k]
        )
        portfolio_value <- portfolio_value + price_k
      }

      for (k in 1:num_asian) {
        # simulated price of k-th Asian option
        price_k <- GeoAsianOptionSim(
          S_full_paths = s_inner_full,
          T = T_inner, r = rf, K = K_asian[k], option_type = type_asian[k],
          position = position_asian[k]
        )
        portfolio_value <- portfolio_value + price_k
      }

      for (k in 1:num_barrier) {
        # simulated price of k-th barrier option
        price_k <- BarrierOptionSim(
          S_full_paths = s_inner_full,
          T = T_inner, r = rf, K = K_barrier[k], vol[s], h,
          H = H_barrier[k], U = U_barrier[k],
          option_type = type_barrier[k], position = position_barrier[k]
        )
        portfolio_value <- portfolio_value + price_k
      }

      s_tauplus1 <- s_inner_full[1, ]
      storage_port_val[[i]] <- portfolio_value
      storage_log_s_tauplus1[[i]] <- 2 * (meanlog_inner[s] - log(s_tauplus1))
      storage_mu_samp[[i]] <- log(s_tau_samp[[i]])
    }

    # optimal pooling of self-normalized LR estimators
    for (i in 1:num_outer) {
      s_tau_i <- s_tau[i]
      mu_i <- log(s_tau_i)

      SNLR <- rep(0, n_nz) # placeholder for self-normalized LR estimators
      for (j in 1:n_nz) {
        mu_samp_j <- storage_mu_samp[[j]]

        LR_j <- exp(((mu_samp_j^2 - mu_i^2) +
          storage_log_s_tauplus1[[j]] * (mu_samp_j - mu_i)) /
          (2 * sdlog_inner[s]^2))
        portfolio_value_j <- storage_port_val[[j]]
        SNLR[j] <- mean(portfolio_value_j * LR_j) / mean(LR_j)
      }
      val_opt_nested[i, m] <- val_opt_nested[i, m] + sum(SNLR * gamma_nz[i, ])
    }
    time_end_sim <- Sys.time()
    time_diff_sim <- time_diff_sim +
      as.numeric(difftime(time_end_sim, time_start_sim, units = c("secs")))
  }

  time_opt_nested_opt[1] <- time_opt_nested_opt[1] + time_diff_opt
  time_opt_nested_opt[2] <- time_opt_nested_opt[2] + time_diff_opt^2

  time_opt_nested_sim[1] <- time_opt_nested_sim[1] + time_diff_sim
  time_opt_nested_sim[2] <- time_opt_nested_sim[2] + time_diff_sim^2

  time_end_mac <- Sys.time()
  time_diff_mac <- as.numeric(difftime(time_end_mac, time_start_mac,
    units = c("secs")
  ))
  time_opt_nested_tot[1] <- time_opt_nested_tot[1] + time_diff_mac
  time_opt_nested_tot[2] <- time_opt_nested_tot[2] + time_diff_mac^2
  print(sprintf(
    "Optimal Nested: This is the %d-th macro, it took %6.2f seconds.",
    m, time_diff_mac
  ))
  print(sprintf("---- Optimization took %6.2f secs.", time_diff_opt))
  print(sprintf("---- Simulation took %6.2f secs.", time_diff_sim))
}

#-- Regression-based Nested Simulation
#   * 1 regression model for the entire portfolio
#   * number of outer scen is the average of the optimal design’s budget

# placeholders for estimated conditional means and runtime
val_reg_nested <- matrix(0, nrow = num_outer, ncol = num_macro)
time_reg_nested_sim <- rep(0, 2)
time_reg_nested_fit <- rep(0, 2)
time_reg_nested_tot <- rep(0, 2)
s_tau_reg_nested <- vector("list", num_macro) # outer scenarios for each macro

num_inner_reg <- 1 # number of inner reps per scen for regression-based design
for (m in 1:num_macro) {
  time_start_mac <- Sys.time()
  time_diff_sim <- 0
  time_diff_fit <- 0

  # number of outer reps for regression-based design
  num_outer_reg <- ceiling(mean(sim_budgets[, m]) / num_inner_reg)

  # simulate outer scenarios in the m-th macro
  z_tau_reg_nested_m <- rmvnorm(num_outer_reg,
    mean = meanlog_outer * num_steps_outer,
    sigma = Sigma * tau
  )

  # num_stocks-by-num_outer matrix
  s_tau_design_m <- sweep(
    x = exp(z_tau_reg_nested_m),
    MARGIN = 2, STATS = S0, FUN = "*"
  )
  s_tau_reg_nested[[m]] <- s_tau_design_m # store design matrix

  time_start_sim <- Sys.time()
  y_design_m <- rep(0, num_outer_reg) # placeholder for portfolio value
  for (i in 1:num_outer_reg) {
    portfolio_value <- 0 # placeholder for portfolio value

    for (s in 1:num_stocks) {
      vol_s <- vol[s]
      s_tau_i <- s_tau_design_m[i, s]

      # growth factor in inner simulation
      z_inner_full <- rnorm(num_steps_inner * num_inner_reg,
        mean = meanlog_inner[s], sd = sdlog_inner[s]
      )
      z_inner_full <- matrix(z_inner_full, nrow = num_steps_inner)
      z_inner_full <- apply(z_inner_full, 2, cumsum)

      s_inner_full <- s_tau_i * exp(z_inner_full) # full inner price paths

      for (k in 1:num_euro) {
        # simulated price of k-th European option
        price_k <- EuroOptionSim(
          S_full_paths = s_inner_full,
          T = T_inner, r = rf, K = K_euro[k], option_type = type_euro[k],
          position = position_euro[k]
        )
        portfolio_value <- portfolio_value + mean(price_k)
      }

      for (k in 1:num_asian) {
        # simulated price of k-th Asian option
        price_k <- GeoAsianOptionSim(
          S_full_paths = s_inner_full,
          T = T_inner, r = rf, K = K_asian[k], option_type = type_asian[k],
          position = position_asian[k]
        )
        portfolio_value <- portfolio_value + mean(price_k)
      }

      for (k in 1:num_barrier) {
        # simulated price of k-th barrier option
        price_k <- BarrierOptionSim(
          S_full_paths = s_inner_full,
          T = T_inner, r = rf, K = K_barrier[k], vol_s, h,
          H = H_barrier[k], U = U_barrier[k],
          option_type = type_barrier[k], position = position_barrier[k]
        )
        portfolio_value <- portfolio_value + mean(price_k)
      }
    }

    y_design_m[i] <- portfolio_value
  }
  time_end_sim <- Sys.time()
  time_diff_sim <- as.numeric(difftime(time_end_sim, time_start_sim,
    units = c("secs")
  ))
  time_reg_nested_sim[1] <- time_reg_nested_sim[1] + time_diff_sim
  time_reg_nested_sim[2] <- time_reg_nested_sim[2] + time_diff_sim^2

  time_start_fit <- Sys.time()

  normalized_s_tau_design_m <- sweep(
    x = s_tau_design_m,
    MARGIN = 2, STATS = S0, FUN = "/"
  ) # design points (normalization for numerical stability)

  # Laguerre polynomials of degree 2 (Longstaff-Schwartz recommendation)
  x_design_m <- Laguerre_basis(normalized_s_tau_design_m)
  # regression model for the entire portfolio
  lm_reg <- lm(y_design_m ~ ., data = data.frame(x_design_m, y_design_m))

  normalized_s_tau_predict_m <- sweep(
    x = s_tau_mac[, , m],
    MARGIN = 2, STATS = S0, FUN = "/"
  ) # target scenarios/prediction points (normalization for numerical stability)
  x_predict_m <- Laguerre_basis(normalized_s_tau_predict_m) # prediction points

  val_reg_nested[, m] <- predict(lm_reg, newdata = x_predict_m)
  time_end_fit <- Sys.time()
  time_diff_fit <- as.numeric(difftime(time_end_fit, time_start_fit,
    units = c("secs")
  ))

  time_reg_nested_fit[1] <- time_reg_nested_fit[1] + time_diff_fit
  time_reg_nested_fit[2] <- time_reg_nested_fit[2] + time_diff_fit^2

  time_end_mac <- Sys.time()
  time_diff_mac <- as.numeric(difftime(time_end_mac, time_start_mac,
    units = c("secs")
  ))
  time_reg_nested_tot[1] <- time_reg_nested_tot[1] + time_diff_mac
  time_reg_nested_tot[2] <- time_reg_nested_tot[2] + time_diff_mac^2
  print(sprintf(
    "Regression: This is the %d-th macro, it took %6.2f seconds.",
    m, time_diff_mac
  ))
  print(sprintf("---- Simulation took %6.2f secs.", time_diff_sim))
  print(sprintf("---- Regression took %6.2f secs.", time_diff_fit))
}

#-- Standard Nested Simulation with Gamma^(1/3) inner reps and Gamma^(2/3)
#   outer scenarios. This method is referred to as SNS in the paper

# placeholders for estimated conditional means and runtime
val_true_nested <- vector("list", num_macro) # true value for MSE calculation
val_std_nested <- vector("list", num_macro) # standard nested estimates
time_std_nested <- rep(0, 2) # runtime for standard nested simulation

s_tau_std_nested <- vector("list", num_macro) # outer scenarios for each macro
sim_budget_std_nested <- matrix(0, nrow = 3, ncol = num_macro) # sim budget
rownames(sim_budget_std_nested) <- c("outer", "inner", "total")

for (m in 1:num_macro) {
  # calculate number of outer and inner reps based on the average budget (over
  # 10 assets) in the optimal design in the same macro
  sim_budget_m <- mean(sim_budgets[, m]) # average budget in the optimal design
  num_outer_std_nested <- ceiling(sim_budget_m^(2 / 3)) # number of outer reps
  num_inner_std_nested <- ceiling(sim_budget_m^(1 / 3)) # number of inner reps
  sim_budget_std_nested[1, m] <- num_outer_std_nested
  sim_budget_std_nested[2, m] <- num_inner_std_nested
  sim_budget_std_nested[3, m] <- num_outer_std_nested * num_inner_std_nested

  val_true_nested_m <- rep(0, num_outer_std_nested) # true value needed for MSE
  val_std_nested_m <- rep(0, num_outer_std_nested) # standard nested estimates

  # simulate outer scenarios in the m-th macro
  z_tau_std_nested_m <- rmvnorm(num_outer_std_nested,
    mean = meanlog_outer * num_steps_outer,
    sigma = Sigma * tau
  )
  # num_stocks-by-num_outer matrix
  s_tau_std_nested_m <- sweep(
    x = exp(z_tau_std_nested_m),
    MARGIN = 2, STATS = S0, FUN = "*"
  )
  s_tau_std_nested[[m]] <- s_tau_std_nested_m # num_stocks-by-num_outer matrix

  for (i in 1:num_outer_std_nested) {
    portfolio_value <- 0 # placeholder for portfolio value

    for (s in 1:num_stocks) {
      vol_s <- vol[s]
      s_tau_i <- s_tau_std_nested_m[i, s]

      for (k in 1:num_euro) {
        # true price of k-th European option
        price_k <- EuroOptionPrice(
          s_tau_i, T_inner, vol_s, rf,
          K_euro[k], type_euro[k], position_euro[k]
        )
        portfolio_value <- portfolio_value + price_k
      }

      for (k in 1:num_asian) {
        # true price of k-th Asian option
        price_k <- AsianOptionPrice(
          s_tau_i, K_asian[k], rf, vol_s,
          T_inner, num_steps_inner, type_asian[k], position_asian[k]
        )
        portfolio_value <- portfolio_value + price_k
      }

      for (k in 1:num_barrier) {
        # true price of k-th barrier option
        price_k <- BarrierOptionPrice_DOP(s_tau_i, K_barrier[k],
          H_barrier[k], rf,
          q = 0, vol_s, h, T_inner, position_barrier[k]
        )
        portfolio_value <- portfolio_value + price_k
      }
    }

    val_true_nested_m[i] <- portfolio_value
  }
  val_true_nested[[m]] <- val_true_nested_m

  time_start <- Sys.time()
  for (i in 1:num_outer_std_nested) {
    portfolio_value <- 0 # placeholder for portfolio value

    for (s in 1:num_stocks) {
      vol_s <- vol[s]
      s_tau_i <- s_tau_std_nested_m[i, s]

      # growth factor in inner simulation
      z_inner_full <- rnorm(num_steps_inner * num_inner_std_nested,
        mean = meanlog_inner[s], sd = sdlog_inner[s]
      )
      z_inner_full <- matrix(z_inner_full, nrow = num_steps_inner)
      z_inner_full <- apply(z_inner_full, 2, cumsum)

      s_inner_full <- s_tau_i * exp(z_inner_full) # full inner price paths

      for (k in 1:num_euro) {
        # simulated price of k-th European option
        price_k <- EuroOptionSim(
          S_full_paths = s_inner_full,
          T = T_inner, r = rf, K = K_euro[k], option_type = type_euro[k],
          position = position_euro[k]
        )
        portfolio_value <- portfolio_value + mean(price_k)
      }

      for (k in 1:num_asian) {
        # simulated price of k-th Asian option
        price_k <- GeoAsianOptionSim(
          S_full_paths = s_inner_full,
          T = T_inner, r = rf, K = K_asian[k], option_type = type_asian[k],
          position = position_asian[k]
        )
        portfolio_value <- portfolio_value + mean(price_k)
      }

      for (k in 1:num_barrier) {
        # simulated price of k-th barrier option
        price_k <- BarrierOptionSim(
          S_full_paths = s_inner_full,
          T = T_inner, r = rf, K = K_barrier[k], vol_s, h,
          H = H_barrier[k], U = U_barrier[k],
          option_type = type_barrier[k], position = position_barrier[k]
        )
        portfolio_value <- portfolio_value + mean(price_k)
      }
    }

    val_std_nested_m[i] <- portfolio_value
  }
  val_std_nested[[m]] <- val_std_nested_m
  time_end <- Sys.time()
  time_diff <- as.numeric(difftime(time_end, time_start, units = c("secs")))
  time_std_nested[1] <- time_std_nested[1] + time_diff
  time_std_nested[2] <- time_std_nested[2] + time_diff^2
  print(sprintf(
    "Standard Nested: This is the %d-th macro, it took %8.2f seconds.",
    m, time_diff
  ))
  print(sprintf(
    "----- There are %d and %d inner and outer reps, resp. Total budget is %d.",
    num_inner_std_nested, num_outer_std_nested,
    num_outer_std_nested * num_inner_std_nested
  ))
}

#-- save results
filename <- sprintf(
  "MultiAssetERM_%dStocks%dOuter%dMacro%dSeed.RData",
  num_stocks, num_outer, num_macro, seed_number
)
save.image(filename)
print(filename)
