# ===== Helper functions
lr2_poisson <- function(lambda_t, lambda_s) {
  # Second moment of likelihood ratios for Poission distributions whose means
  # are lambda_t and lambda_s for the target and sampling distributions, resp.
  # Input:
  # lambda_t, lambda_s: vector of (non-negative) means
  # Output:
  # second moment of the likelihood ratio between the target and the
  # sampling distribution
  # formula: exp(lambda_t^2/lambda_s+lambda_s-2*lambda_t)

  n_target <- length(lambda_t) # number of target params
  n_sample <- length(lambda_s) # number of sampling params

  lambda_t <- matrix(lambda_t, nrow = n_target, ncol = n_sample)
  lambda_s <- matrix(lambda_s, nrow = n_target, ncol = n_sample, byrow = TRUE)

  return(exp(lambda_t^2 / lambda_s + lambda_s - 2 * lambda_t))
}

lr2_poisson_sim <- function(lambda_t, lambda_s, n_reps = 1e6) {
  # compute the sample second moment of likelihood ratios for Poission distributions
  # whose means are lambda_t and lambda_s for the target and sampling
  # distributions, resp.
  # Input:
  # -- lambda_t, lambda_s: vector of (non-negative) means
  # -- n_reps: (optional) nubmer of independent replications
  # Output:
  # -- second moment of the likelihood ratio between the target and the
  # sampling distribution
  # formula: exp(lambda_t^2/lambda_s+lambda_s-2*lambda_t)
  sample <- rpois(n_reps, lambda_s)

  return(mean((dpois(sample, lambda_t) / dpois(sample, lambda_s))^2))
}

psi_mat_poisson <- function(lambda_t, lambda_s = NULL) {
  # Matrix of the reciprocals of the second moments of the likelihood ratios
  # for Poisson distributions whose means equal to lambda_t and lambda_s
  # Input:
  # lambda_t, lambda_s: matrices of (non-negative) means, n_dim values in each
  #                     row
  # Output:
  # matrix of the inverses of second moments of the likelihood ratio between
  # the target and the sampling distributions

  # if no sampling param is specify, use the target params as the default
  if (is.null(lambda_s)) {
    lambda_s <- lambda_t
  }

  n_dim <- ncol(lambda_t) # number of dimensions (input dist.)
  n_target <- nrow(lambda_t) # number of target params
  n_sample <- nrow(lambda_s) # number of sampling params

  psi_mat <- matrix(1, n_target, n_sample) # placeholder, row for target params

  for (i in 1:n_dim) {
    # loop over number of input distributions
    # outer product calculation for efficiency
    psi_mat <- psi_mat * lr2_poisson(lambda_t[, i], lambda_s[, i])
  }
  return(1 / psi_mat)
}

newsvendor_profit_oracle <- function(mean_demands, policy, price, cost) {
  # Computes the daily profit of the newsvendor given the mean demands, stocking
  # policy, sale price, and unit cost for each product
  # Input:
  # -- mean_demands: vector of mean demands of Poisson distributions
  # -- policy: vector stocking policy of each product
  # -- price: unit sale price of each product
  # -- cost: unit cost for each product
  # Output:
  # -- daily profit of the newsvendor

  n_products <- length(mean_demands)

  total <- -cost * sum(policy)

  for (i in 1:n_products) {
    total <- total +
      price[i] * ((1:policy[i]) %*% dpois((1:policy[i]), mean_demands[i]) +
        policy[i] * (1 - ppois(policy[i], mean_demands[i])))
  }

  return(total)
}

newsvendor_profit_simulator <- function(
    mean_demands, policy, price, cost,
    n_reps) {
  # This function simulates the daily profit of the newsvendor.
  # Input:
  # -- mean_demands: vector of mean demands of Poisson distributions
  # -- policy: vector stocking policy of each product
  # -- price: unit sale price of each product
  # -- cost: unit cost for each product
  # -- n_reps: number of replications
  # Output:
  # -- list of simulated profits, demands, and joint densities

  n_products <- length(mean_demands)

  # -- placeholders for inner sim. inputs, joint density, and outputs
  # inner sim. input each column corresponds to a product
  inner_input <- matrix(0, nrow = n_reps, ncol = n_products)
  joint_dens <- rep(1, n_reps) # inner sim. joint density
  inner_output <- rep(-sum(cost * policy), n_reps) # needs further calculations

  # simulate the demands for each product
  for (i in 1:n_products) {
    inner_input[, i] <- rpois(n_reps, mean_demands[i])
    joint_dens <- joint_dens * dpois(inner_input[, i], mean_demands[i])
  }

  # calculate the profit for each replication
  for (rep in 1:n_reps) {
    inner_output[rep] <- inner_output[rep] +
      sum(price %*% pmin(policy, inner_input[rep, ]))
  }

  return(list(
    profits = inner_output,
    inner_rv = inner_input,
    density = joint_dens
  ))
}

newsvendor_opt_design <- function(lambda_outer, n0) {
  # Computes the optimal design, i.e., optimal sampling and pooling decisions,
  # for the newsvendor problem. The optimal design is based on the joint
  # distribution of all products' demands, not each product individually.
  # Input:
  # -- lambda_outer: matrix of outer scenarios, each row is a scenario, each
  #                  column is a product.
  # -- n0: desired number of inner replications (as in the standard design)
  # Output:
  # -- total number of inner replications for the optimal design

  n_outer <- nrow(lambda_outer) # number of outer scenarios
  # n_outer-by-n_outer matrix of 1/E[W^2]
  psi_mat <- psi_mat_poisson(lambda_outer)
  opt_results <- lp(
    direction = "min",
    objective.in = rep(1, n_outer), # obj. coeff, all ones
    const.mat = psi_mat,
    const.dir = rep(">=", n_outer), # constraint direction
    const.rhs = rep(1, n_outer), # RHS, NInner
    all.int = FALSE # round up after solving the LP
  )
  n_j <- opt_results$solution # (non-integral) optimal number of inner reps
  ind_nz <- which(n_j > 1e-8) # indices of non-zero-sample scenarios
  n_j_nz <- n_j[ind_nz] # non-zero-sample sizes
  n_j_opt_nested <- ceiling(n_j_nz * n0) # ceiling of optimal inner reps

  # optimal pooling weights for the optimized sampling plan
  psi_mat_nz <- psi_mat[, ind_nz] # non-zero-sample columns of the Psi matrix
  gamma_nz <- sweep(psi_mat_nz, MARGIN = 2, n_j_opt_nested, "*")
  gamma_nz <- gamma_nz / rowSums(gamma_nz) # optimal pooling weights

  return(list(
    ind_nz = ind_nz,
    n_nz_inner = n_j_opt_nested,
    sim_budget = sum(n_j_opt_nested),
    gamma = gamma_nz
  ))
}

opt_pool <- function(scen, gammas, sim_output) {
  # optimal pooling of simulation outputs via self-normalized likelihood ratio
  # estimators according to the given gamma weights
  # Input:
  # -- scen: vector of target parameters
  # -- gammas: vector of pooling weights
  # -- sim_output: list of simulation outputs
  # Output:
  # -- optimally pooled self-normalized LR estimates

  n_sample <- length(sim_output)

  self_norm_lr <- rep(0, n_sample) # self-normalized LR estimator placeholder
  for (i in 1:n_sample) {
    # calculate numerator using target scenario density
    numer <- apply(
      sweep(sim_output[[i]]$inner_rv, MARGIN = 2, scen, dpois),
      1, prod
    )
    likelihood_ratio <- numer / sim_output[[i]]$density
    self_norm_lr[i] <- mean(sim_output[[i]]$profits * likelihood_ratio) /
      mean(likelihood_ratio)
  }

  return(sum(self_norm_lr * gammas))
}

newsvendor_profit_opt_nested <- function(
    lambda_outer,
    policy, price, cost,
    n0) {
  # Estimates the profits in the newsvendor problem using the optimal LR nested
  # simultion design
  # Input:
  # -- lambda_outer: matrix of outer scenarios, each row is a scenario, each
  #                  column is a product
  # -- policy: vector stocking policy of each product
  # -- price: unit sale price of each product
  # -- cost: unit cost for each product
  # -- n0: desired number of inner replications (as in the standard design)
  # Output:
  # -- list of estimated profits and the total simulation budget

  opt_design <- newsvendor_opt_design(lambda_outer, n0)

  ind_nz <- opt_design$ind_nz # indices of non-zero-sample scenarios
  n_samp <- length(ind_nz) # number of non-zero-sample scenarios
  n_inner_samp <- opt_design$n_nz_inner # non-zero-sample sizes
  gamma_samp <- opt_design$gamma # optimal pooling weights

  storage <- vector("list", length = n_samp)
  for (i in 1:n_samp) {
    storage[[i]] <- newsvendor_profit_simulator(
      lambda_outer[ind_nz[i], ],
      policy, price, cost,
      n_reps = n_inner_samp[i]
    )
  }

  # convert the lambda_outer and gamma_samp matrices to lists for mapply
  lambda_outer_lst <- as.list(data.frame(t(lambda_outer)))
  gamma_samp_lst <- as.list(data.frame(t(gamma_samp)))
  profits <- mapply(
    function(x, y) opt_pool(x, y, storage),
    lambda_outer_lst, gamma_samp_lst
  )
  names(profits) <- NULL

  return(list(profits = profits, sim_budget = opt_design$sim_budget))
}

newsvendor_profit_reg_nested <- function(
    lambda_outer_pred,
    policy, price, cost,
    lambda_outer_reg, n_inner_reg = 1) {
  # This function estimates the profits in the newsvendor problem using the
  # regression-based nested simulation design
  # Input:
  # -- lambda_outer_pred: matrix of target outer scenarios for prediction
  # -- policy: vector stocking policy of each product
  # -- price: unit sale price of each product
  # -- cost: unit cost for each product
  # -- lambda_outer_reg: matrix of outer scenarios for regression
  # -- n_inner_reg: number of inner replications for regression
  # Output:
  # -- vector of estimated profits

  # simulate profits for the regression design points
  profits_reg <- apply(
    lambda_outer_reg, 1,
    function(x) {
      mean(newsvendor_profit_simulator(x,
        policy, price, cost,
        n_reps = n_inner_reg
      )$profits)
    }
  )

  # fit a regression model
  lm_reg <- lm(profits_reg ~ 1 + I(lambda_outer_reg) + I(lambda_outer_reg^2))

  # predict the profits for the prediction design points, i.e., target scenarios
  profits <- predict(
    lm_reg,
    data.frame(lambda_outer_reg = I(lambda_outer_pred))
  )
  names(profits) <- NULL
  return(profits)
}

newsvendor_var_oracle <- function(mean_demands, policy, price) {
  # Calculates the exact conditional variance given the product mean demands.
  # Input:
  # -- mean_demands: vector of mean demands of Poisson distributions
  # -- policy: vector stocking policy of each product
  # -- price: unit sale price of each product
  # Output:
  # -- conditional variance of the profit

  n_products <- length(mean_demands)

  total <- 0
  for (i in 1:n_products) {
    var_product <- policy[i]^2
    mean_product <- policy[i]
    for (j in 0:policy[i]) {
      var_product <- var_product +
        (j^2 - policy[i]^2) * dpois(j, mean_demands[i])
      mean_product <- mean_product +
        (j - policy[i]) * dpois(j, mean_demands[i])
    }
    var_product <- var_product - mean_product^2
    total <- total + price[i]^2 * var_product
  }
  return(total)
}

between <- function(x, a, b) {
  return(x >= a & x <= b)
}
