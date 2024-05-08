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

EuroOptionPrice <- function(S, T, sigma, r, K, option_type, position) {
  # Calculates the Black-Scholes price of a European option.
  # Input:
  # S: current asset price
  # T: time to maturity
  # sigma: volatility
  # r: risk-free rate
  # K: strike price
  # option_type: "C" or "P"
  # position: integer, positive for long, negative for short
  # Output:
  # price of the option

  d1 <- (1 / (sigma * sqrt(T))) * (log(S / K) + (r + sigma^2 / 2) * T)
  d2 <- d1 - (sigma * sqrt(T))

  if (option_type == "c") {
    price <- pnorm(d1) * S - pnorm(d2) * K * exp(-r * T)
  } else {
    price <- pnorm(-d2) * K * exp(-r * T) - pnorm(-d1) * S
  }

  return(price * position)
}

AsianOptionPrice <- function(
    S, K, r, sigma, T, n,
    option_type = "c", position = 1) {
  # Calculates the Black-Scholes price of a fixed-strike geometric Asian option.
  # Geometric average without the initial price.
  # Inputs:
  # S: initial stock price
  # K: strike price
  # r: risk-free interest rate
  # sigma: annualized volatility
  # T: time to maturity
  # n: number of subintervals
  # option_type: "c" and "p" for fixed strike call and put, respectively.
  # position: integer, positive for long, negative for short
  # Output:
  # price of the option

  sigma_Z <- (sigma / n) * sqrt((n + 1) * (2 * n + 1) / 6)
  rho <- (r - 0.5 * sigma^2) * (n + 1) / (2 * n) + 0.5 * sigma_Z^2

  d1 <- (log(S / K) + (rho + 0.5 * sigma_Z^2) * T) / (sigma_Z * sqrt(T))
  d2 <- d1 - sigma_Z * sqrt(T)

  if (option_type == "c") {
    price <- S * exp((rho - r) * T) * pnorm(d1) - K * exp(-r * T) * pnorm(d2)
  } else if (option_type == "p") {
    price <- K * exp(-r * T) * pnorm(-d2) - S * exp((rho - r) * T) * pnorm(-d1)
  } else {
    stop("type must be either 'c' or 'p'")
  }

  return(price * position)
}

compute_abcd <- function(upper, lower, r, q, t, sigma) {
  res_minus <- (log(upper / lower) + (r - q) * t) / (sigma * sqrt(t)) -
    (1 / 2) * sigma * sqrt(t)
  res_plus <- res_minus + sigma * sqrt(t)
  return(list(res_minus, res_plus))
}

BarrierOptionPrice_DOP <- function(S, K, H, r, q, sigma, t, T, position) {
  # Formula are shown above
  at_minus <- compute_abcd(S, K, r, q, t, sigma)[[1]]
  at_plus <- compute_abcd(S, K, r, q, t, sigma)[[2]]
  aT_minus <- compute_abcd(S, K, r, q, T, sigma)[[1]]
  aT_plus <- compute_abcd(S, K, r, q, T, sigma)[[2]]
  bt_minus <- compute_abcd(S, H, r, q, t, sigma)[[1]]
  bt_plus <- compute_abcd(S, H, r, q, t, sigma)[[2]]
  bT_minus <- compute_abcd(S, H, r, q, T, sigma)[[1]]
  bT_plus <- compute_abcd(S, H, r, q, T, sigma)[[2]]
  ct_minus <- compute_abcd(H, S, r, q, t, sigma)[[1]]
  ct_plus <- compute_abcd(H, S, r, q, t, sigma)[[2]]
  cT_minus <- compute_abcd(H, S, r, q, T, sigma)[[1]]
  cT_plus <- compute_abcd(H, S, r, q, T, sigma)[[2]]
  dt_minus <- compute_abcd(H^2, S * K, r, q, t, sigma)[[1]]
  dt_plus <- compute_abcd(H^2, S * K, r, q, t, sigma)[[2]]
  dT_minus <- compute_abcd(H^2, S * K, r, q, T, sigma)[[1]]
  dT_plus <- compute_abcd(H^2, S * K, r, q, T, sigma)[[2]]

  rho <- sqrt(t / T)
  k <- 2 * (r - q) / (sigma^2)

  price <- K * exp(-r * T) * (pmvnorm(
    lower = c(-Inf, -Inf), upper = c(bt_minus, -aT_minus), mean = c(0, 0),
    matrix(c(1, -rho, -rho, 1), nrow = 2, byrow = TRUE)
  )
  - pmvnorm(
      lower = c(-Inf, -Inf), upper = c(bt_minus, -bT_minus), mean = c(0, 0),
      matrix(c(1, -rho, -rho, 1), nrow = 2, byrow = TRUE)
    )) -
    S * exp(-q * T) * (pmvnorm(
      lower = c(-Inf, -Inf), upper = c(bt_plus, -aT_plus), mean = c(0, 0),
      matrix(c(1, -rho, -rho, 1), nrow = 2, byrow = TRUE)
    )
    - pmvnorm(
        lower = c(-Inf, -Inf), upper = c(bt_plus, -bT_plus), mean = c(0, 0),
        matrix(c(1, -rho, -rho, 1), nrow = 2, byrow = TRUE)
      )) -
    (H / S)^(k - 1) * K * exp(-r * T) * (pmvnorm(
      lower = c(-Inf, -Inf), upper = c(-ct_minus, -dT_minus), mean = c(0, 0),
      matrix(c(1, rho, rho, 1), nrow = 2, byrow = TRUE)
    )
    - pmvnorm(
        lower = c(-Inf, -Inf), upper = c(-ct_minus, -cT_minus), mean = c(0, 0),
        matrix(c(1, rho, rho, 1), nrow = 2, byrow = TRUE)
      )) +
    (H / S)^(k + 1) * S * exp(-q * T) * (pmvnorm(
      lower = c(-Inf, -Inf), upper = c(-ct_plus, -dT_plus), mean = c(0, 0),
      matrix(c(1, rho, rho, 1), nrow = 2, byrow = TRUE)
    )
    - pmvnorm(
        lower = c(-Inf, -Inf), upper = c(-ct_plus, -cT_plus), mean = c(0, 0),
        matrix(c(1, rho, rho, 1), nrow = 2, byrow = TRUE)
      ))

  return(price * position)
}

EuroOptionSim <- function(S_full_paths, T, r, K, option_type, position) {
  # Simulate the payoff of a European option for given stock price paths
  # Inputs:
  # S_full_paths: matrix of stock paths, each column is a path, each row is a
  #               time point
  # T: time to maturity
  # r: risk-free rate
  # K: strike price
  # option_type: option types "c" for call, "p" for put
  # position: position of the option, positive for long, negative for short
  # Output:
  # vector of discounted payoffs of the specified option

  num_steps <- nrow(S_full_paths)

  # S_T: stock price at maturity
  S_T <- S_full_paths[num_steps, ]

  if (option_type == "c") {
    # call option
    payoff <- pmax(S_T - K, 0) * position
  } else if (option_type == "p") {
    # put option
    payoff <- pmax(K - S_T, 0) * position
  } else {
    # give a warning message that the type should be either "c" or "p"
    warning("The type of the option should be either 'c' or 'p'.")
  }

  return(exp(-r * T) * payoff)
}

GeoAsianOptionSim <- function(S_full_paths, T, r, K, option_type, position) {
  # Simulate the payoff of a fixed-strike geometric Asian option for given stock
  # price paths
  # Inputs:
  # S_full_paths: matrix of stock paths, each column is a path, each row is a
  # time point
  # T: time to maturity
  # r: risk-free rate
  # K: strike price
  # option_type: option types "c" for call, "p" for put
  # position: position of the option, positive for long, negative for short
  # Output:
  # vector of discounted payoffs of the specified option

  num_steps <- nrow(S_full_paths)

  # S_bar: geometric average stock price
  S_bar <- apply(S_full_paths, 2, prod)^(1 / num_steps)

  if (option_type == "c") {
    # call option
    payoff <- pmax(S_bar - K, 0) * position
  } else if (option_type == "p") {
    # put option
    payoff <- pmax(K - S_bar, 0) * position
  } else {
    # give a warning message that the type should be either "c" or "p"
    warning("The type of the option should be either 'c' or 'p'.")
  }

  return(exp(-r * T) * payoff)
}

BarrierOptionSim <- function(
    S_full_paths, T, r, K, vol, h,
    H = -Inf, U = Inf, option_type, position) {
  # Simulate the payoff of a barrier option for given stock price paths
  # Inputs:
  # S_full_paths: matrix of stock paths, each column is a path, each row is a
  # time point
  # T: time to maturity
  # r: risk-free rate
  # K: strike price
  # option_type:
  # pdo: down-and-out put
  # cdo: down-and-out call
  # pdi: down-and-in put
  # cdi: down-and-in call
  # puo: up-and-out put
  # cuo: up-and-out call
  # pui: up-and-in put
  # cui: up-and-in call
  # position: position of the option, positive for long, negative for short
  # Output:
  # vector of discounted payoffs of the specified option

  num_steps <- nrow(S_full_paths)
  num_paths <- ncol(S_full_paths)

  # S_T: stock price at maturity
  S_T <- S_full_paths[num_steps, ]

  # check if letter d is in each of option_type
  if (grepl("d", option_type)) {
    # minimum price of the underlying asset if lower barrier is present
    s_min <- matrix(0, nrow = num_steps - 1, ncol = num_paths)
    for (s in 2:num_steps) {
      b1 <- log(S_full_paths[s, ])
      b2 <- log(S_full_paths[s - 1, ])
      U <- runif(num_paths)
      s_min[s - 1, ] <- exp(((b1 + b2) -
        sqrt((b1 - b2)^2 - 2 * vol^2 * h * log(U))) / 2)
    }
    s_min <- apply(s_min, 2, min)
  } else if (grepl("u", option_type)) {
    # maximum price of the underlying asset if upper barrier is present
    s_max <- matrix(0, nrow = num_steps - 1, ncol = num_paths)
    for (s in 2:num_steps) {
      b1 <- log(S_full_paths[s, ])
      b2 <- log(S_full_paths[s - 1, ])
      U <- runif(num_paths)
      s_max[s - 1, ] <- exp(((b1 + b2) +
        sqrt((b1 - b2)^2 - 2 * vol^2 * h * log(U))) / 2)
    }
    s_max <- apply(s_max, 2, max)
  }

  if (option_type == "pdo") {
    # down-and-out put
    payoff <- pmax(K - S_T, 0) * (s_min > H) * position
  } else if (option_type == "cdo") {
    # down-and-out call
    payoff <- pmax(S_T - K, 0) * (s_min > H) * position
  } else if (option_type == "pdi") {
    # down-and-in put
    payoff <- pmax(K - S_T, 0) * (s_min < H) * position
  } else if (option_type == "cdi") {
    # down-and-in call
    payoff <- pmax(S_T - K, 0) * (s_min < H) * position
  } else if (option_type == "puo") {
    # up-and-out put
    payoff <- pmax(K - S_T, 0) * (s_max < U) * position
  } else if (option_type == "cuo") {
    # up-and-out call
    payoff <- pmax(S_T - K, 0) * (s_max < U) * position
  } else if (option_type == "pui") {
    # up-and-in put
    payoff <- pmax(K - S_T, 0) * (s_max > U) * position
  } else if (option_type == "cui") {
    # up-and-in call
    payoff <- pmax(S_T - K, 0) * (s_max > U) * position
  } else {
    # give a warning message that the type should be either "pdo", "cdo", "pdi",
    # "cdi", "puo", "cuo", "pui", or "cui"
    warning("The type of the option should be either 'pdo', 'cdo', 'pdi', 'cdi',
     'puo', 'cuo', 'pui', or 'cui'.")
  }

  return(exp(-r * T) * payoff)
}

Laguerre_basis <- function(X) {
  # Calculate the first three Laguerre polynomials
  # Input:
  # X: input values
  # Output:
  # data frame with columns L0, L1, L2

  L0 <- exp(-X / 2)
  L1 <- exp(-X / 2) * (1 - X)
  L2 <- exp(-X / 2) * (1 / 2) * (X^2 - 4 * X + 2)
  # L3 <- (1 / 6) * (-X^3 + 9 * X^2 - 18 * X + 6)
  # L4 <- (1/24)*(X^4 - 16*X^3 + 72*X^2 - 96*X + 24)

  return(as.data.frame(cbind(L0, L1, L2)))
}

Qntile <- function(df, probs) {
  return(t(apply(df, 1, function(x) quantile(x, probs))))
}

MSE <- function(df, dfTrue) {
  return(rowMeans((df - dfTrue)^2))
}

Variance <- function(df, dfTrue) {
  return(rowMeans((df - rowMeans(df))^2))
}

MAE <- function(df, dfTrue) {
  return(rowMeans(abs(df - dfTrue)))
}
