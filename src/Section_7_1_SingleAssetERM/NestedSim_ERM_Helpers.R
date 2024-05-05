lr2_lnorm <- function(meanlog0 = 0, sdlog0 = 1, meanlog1 = 0, sdlog1 = 1) {
  # Second moment of likelihood ratio between two log normal distributions.
  # The target distribution's logarithm has mean equal to meanlog0 and
  # standard deviation equal to sdlog0. The sampling distribution's logarithm
  # has mean equal to meanlog1 and standard deviation equal to sdlog1.
  # Input:
  # meanlog0, sdlog0, meanlog1, sdlog1: mean and standard deviation of the
  #       target and the sampling distributions on the log scale with defaul
  #       values of 0 and 1 respectively.
  # Output:
  # second moment of the likelihood ratio between the target and the
  # sampling distribution
  # formula:
  # exp((mu0-mu1)^2/(2sig0^2-sig1^2))*(sig0^2/(sig1*sqrt(2sig0^2-sig1^2)))

  sig2_range <- 2 * sdlog1^2 - sdlog0^2
  if (sig2_range > 0) {
    return(exp((meanlog1 - meanlog0)^2 / sig2_range) *
      (sdlog1^2 / (sdlog0 * sqrt(sig2_range))))
  } else {
    warning("Parameter provided produces infinite second moment.")
    return(Inf)
  }
}

psi_mat_lnorm <- function(meanlogs, sdlogs) {
  # Matrix of the inverses the second moments of likelihood ratios between
  # lognormal distributions whose logarithms have means equal to meanlogs
  # and standard deviation equal to sdlogs
  # Input:
  # meanlogs: n_outer vector of means of the lognormal distributions
  # sdlogs: n_outer vector of standard deviations of the lognormal distributions
  # Output:
  # n_outer-by-n_outer matrix of the inverses of second moments of likelihood
  # ratios between the lognormal distributions. Matrix is NOT symmetric.

  if (length(sdlogs) == 1) {
    # in the ERM example, sdlogs is a scalar so psi_mat can be quickly computed
    sig2_range <- sdlogs^2
    psi_mat <- outer(
      meanlogs, meanlogs,
      function(x, y) exp(-(x - y)^2 / sig2_range)
    )
  } else {
    # if sdlogs is a vector, need to compute the matrix elements one by one
    n_outer <- length(meanlogs)
    psi_mat <- matrix(0, n_outer, n_outer)

    for (t in 1:n_outer) {
      for (s in 1:n_outer) {
        psi_mat[t, s] <- lr2_lnorm(
          meanlogs[t], sdlogs[t],
          meanlogs[s], sdlogs[s]
        )
      }
    }
    psi_mat <- 1 / psi_mat
  }
  return(psi_mat)
}

bs_price_euro_option <- function(
    s0, str_k, div, rf, vol, big_t,
    option_type = "c", position = 1) {
  # Calculates the Black-Scholes price of a European option.
  # s0: initial stock price
  # str_k: strike price
  # div: dividend rate
  # rf: risk-free rate
  # vol: volatility
  # big_t: time to maturity
  # option_type: "c" or "p" for call or put, respectively
  # position: integer, (+) and (-) for long and short, respectively

  s_disc <- s0 * exp(-div * big_t) # forward price for stock
  str_k_disc <- str_k * exp(-rf * big_t) # forward price for strike

  d1 <- (log(s_disc / str_k_disc) + (vol^2 / 2) * big_t) / (vol * sqrt(big_t))
  d2 <- d1 - vol * sqrt(big_t)

  if (option_type == "c") { # call option
    price <- s_disc * pnorm(d1) - str_k_disc * pnorm(d2)
  } else if (option_type == "p") { # put option
    price <- str_k_disc * pnorm(-d2) - s_disc * pnorm(-d1)
  } else {
    stop("Invalid option type. Please use 'c' or 'p'.")
  }

  return(price * position)
}

bs_price_straddle <- function(s0, str_k, div, rf, vol, big_t) {
  # Black-Scholes price for straddle option portfolio
  # Inputs:
  # s0: initial stock price
  # str_k: strike price
  # div: dividend rate
  # rf: risk-free rate
  # vol: volatility
  # big_t: time to maturity
  # Output:
  # price of the straddle option portfolio

  # call and put prices
  call_price <- bs_price_euro_option(s0, str_k, div, rf, vol, big_t, "c")
  put_price <- bs_price_euro_option(s0, str_k, div, rf, vol, big_t, "p")

  return(call_price + put_price) # return straddle price
}

bs_sim_straddle <- function(s0, str_k, div, rf, vol, big_t, n_rep) {
  # Simulated price for straddle option strategy under the Black-Scholes model
  # Inputs:
  # s0: initial stock price
  # str_k: strike price
  # div: dividend rate
  # rf: risk-free rate
  # vol: volatility
  # big_t: time to maturity
  # n_rep: number of indepenDent replications
  # Output:
  # simulated price of the straddle option portfolio

  meanlog <- (rf - div - vol^2 / 2) * (big_t) # normal mean in lognormal price
  sdlog <- vol * sqrt(big_t) # normal sd in lognormal price
  s_big_t <- rlnorm(n_rep, log(s0) + meanlog, sdlog) # projected terminal price
  val <- exp(-rf * big_t) * abs(s_big_t - str_k) # discounted terminal payoff
  den <- dlnorm(s_big_t, log(s0) + meanlog, sdlog) # sampling Density

  return(list(val = val, rv = s_big_t, den = den))
}

straddle_std_nested <- function(
    s_tau, str_k, div, rf,
    vol, big_t, tau, n_inner) {
  # Standard nested simulation for valuation of straddle option portfolio
  # Inputs:
  # s_tau: vector of projected stock prices at time tau
  # str_k: strike price
  # div: dividend rate
  # rf: risk-free rate
  # vol: volatility
  # big_t: time to maturity
  # tau: risk horizon
  # n_inner: number of inner replications per outer scenario
  # Output:
  # vector of estimated conditional means for s_tau

  n_outer <- length(s_tau) # number of outer scenarios
  val <- rep(0, n_outer) # estimated conditional means for s_tau

  for (i in 1:n_outer) { # loop over outer scenarios
    val[i] <- mean(bs_sim_straddle(
      s_tau[i], str_k, div, rf,
      vol, (big_t - tau), n_inner
    )$val) # estimated conditional mean
  }

  return(val) # vector of estimated conditional means for s_tau
}

straddle_opt_nested <- function(
    s_tau, str_k, div, rf,
    vol, big_t, tau,
    s_tau_samp, n_inner_samp, gamma_samp) {
  # Optimal nested simulation for valuation of straddle option portfolio
  # Inputs:
  # s_tau: vector of projected stock prices at time tau, i.e., target scenarios
  # str_k: strike price
  # div: dividend rate
  # rf: risk-free rate
  # vol: volatility
  # big_t: time to maturity
  # tau: risk horizon
  # s_tau_samp: vector of sampling scens, i.e., where inner sims are performed
  # n_inner_samp: vector of number of inner reps in each sampling scen
  # gamma_samp: vector of pooling weights for each sampling scen
  # Output:
  # vector of estimated conditional means for s_tau

  n_outer <- length(s_tau) # number of outer scenarios
  # inner sim. normal mean
  meanlog_inner <- (rf - div - vol^2 / 2) * (big_t - tau)
  sdlog_inner <- vol * sqrt(big_t - tau) # inner sim. normal sd
  val <- rep(0, n_outer) # placeholder of estimated conditional means

  n_samp <- length(s_tau_samp) # number of sampling scenarios
  storage_rv <- vector("list", n_samp) # placeholder for rv in LR calculation
  storage_val <- vector("list", n_samp) # placeholder for portfolio value
  storage_den <- vector("list", n_samp) # placeholder for sampling density
  for (i in 1:n_samp) { # loop over sampling scenarios
    # simulate straddle option portfolio
    storage_i <- bs_sim_straddle(
      s_tau_samp[i], str_k, div, rf,
      vol, (big_t - tau), n_inner_samp[i]
    )

    storage_rv[[i]] <- storage_i$rv # store the rv in LR calculation
    storage_val[[i]] <- storage_i$val # store the portfolio value
    storage_den[[i]] <- storage_i$den # store the sampling density
  }

  # optimal pooling of self-normalized LR estimators
  for (str_k in 1:n_outer) {
    s_tau_k <- s_tau[str_k] # str_k-th (target) outer scenario
    sn_lr <- rep(0, n_samp) # placeholder for self-normalized LR estimators
    for (j in 1:n_samp) {
      # numerator of LR, i.e., evaluate sampling rv's using target pdf
      numer_j <- dlnorm(
        storage_rv[[j]],
        log(s_tau_k) + meanlog_inner, sdlog_inner
      )
      lr_j <- numer_j / storage_den[[j]] # likelihood ratio

      # self-normalized LR estimator for j-th sampling scenario
      sn_lr[j] <- mean(storage_val[[j]] * lr_j) / mean(lr_j)
    }
    # pooling of self-normalized LR estimators using optimal weights
    val[str_k] <- sum(sn_lr * gamma_samp[str_k, ])
  }

  return(val)
}

straddle_reg_nested <- function(s_tau, str_k, div, rf, vol, big_t, tau,
                                s0, s_tau_design, n_inner_design = 1) {
  # Regression-based nested simulation for valuation of straddle option
  # portfolio
  # Inputs:
  # s_tau: vector of projected stock prices at time tau, i.e., target scenarios
  # str_k: strike price
  # div: dividend rate
  # rf: risk-free rate
  # vol: volatility
  # big_t: time to maturity
  # tau: risk horizon
  # s0: initial stock price
  # s_tau_design: vector of design points, can be different than s_tau
  # n_inner_design: vector of number of inner reps in each design point
  # Output:
  # vector of estimated conditional means for s_tau

  n_outer <- length(s_tau) # number of outer scenarios (prediction points)
  val <- rep(0, n_outer) # placeholder for estimated conditional means

  n_design <- length(s_tau_design) # number of design points
  # design points (normalization for numerical stability)
  x_design <- s_tau_design / s0
  y_design <- rep(0, n_design) # estimated conditional means at design points

  for (i in 1:n_design) {
    y_design[i] <- mean(bs_sim_straddle(
      s_tau_design[i], str_k, div, rf, vol, (big_t - tau),
      n_inner_design
    )$val)
  }

  # Laguerre polynomials of degree 2 (Longstaff-Schwartz recommendation)
  x_0 <- exp(-x_design / 2)
  x_1 <- exp(-x_design / 2) * (1 - x_design)
  x_2 <- exp(-x_design / 2) * (1 - 2 * x_design + x_design^2 / 2)

  # regression coefficients
  coef_reg <- summary(lm(y_design ~ x_0 + x_1 + x_2))$coefficients

  # prediction points (normalization for numerical stability)
  x_predict <- s_tau / s0
  x_0 <- exp(-x_predict / 2)
  x_1 <- exp(-x_predict / 2) * (1 - x_predict)
  x_2 <- exp(-x_predict / 2) * (1 - 2 * x_predict + x_predict^2 / 2)
  val <- coef_reg[1, 1] + coef_reg[2, 1] * x_0 +
    coef_reg[3, 1] * x_1 + coef_reg[4, 1] * x_2

  return(val) # vector of estimated conditional means
}

erm_mse_func <- function(val_true, val_est) {
  # Calculate bias, variance, and MSE for risk measures, such as VaR,
  # indicator function, hockey-stick function, and squared losses,
  # estimated by different nested simulation designs such as OptNested,
  # StdNested, StdPlus, and RegNested.
  # Inputs:
  #   val_true: a vector of true values for different risk measures
  #   val_est: a matrix of estimated values for different risk measures
  # Outputs:
  #   a data frame of bias, variance, and MSE for different risk measures

  return(data.frame(rbind(
    Bias = colMeans(sweep(val_est, 2, val_true, "-")),
    Var = colMeans((sweep(val_est, 2, colMeans(val_est), "-"))^2),
    MSE = colMeans((sweep(val_est, 2, val_true, "-"))^2)
  )))
}

erm_risk_measures <- function(val_mat, alpha, threshold) {
  # Calculate different risk measures, such as VaR, indicator function,
  # hockey-stick function, and squared losses, estimated by different
  # nested simulation designs such as OptNested, StdNested, StdPlus,
  # and RegNested.
  # Inputs:
  #   val_mat: an NOut-by-n_macro matrix of estimated conditional means
  #   alpha: a scalar for the quantile level for VaR calculation
  #   threshold: a scalar for the threshold level for tail risk measures
  # Outputs:
  #   A n_outer-by-4 matrix of different risk measures

  # VaR
  VaR <- apply(val_mat, 2, function(x) quantile(x, alpha))
  # indicator function
  indicator <- colMeans(sweep(val_mat, 2, threshold, ">"))
  # hockey-stick function
  hockeystick <- colMeans(pmax(sweep(val_mat, 2, threshold, "-"), 0))
  # squared loss
  square <- colMeans((pmax(sweep(val_mat, 2, threshold, "-"), 0))^2)

  # return a matrix
  return(cbind(VaR, indicator, hockeystick, square))
}
