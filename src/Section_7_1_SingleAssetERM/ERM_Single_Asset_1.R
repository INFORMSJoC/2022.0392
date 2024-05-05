################################################################################
# ERM_Single_Asset_1.R                                                         #
# Single-Asset ERM example in Section 7.1 in the paper                         #
# Straddle option portfolio that consists of a long call and a long put        #
# Optimized nested simulation for enterprise risk management (ERM) example     #
# -- straddle option (long call + long put,the same strike for both options)   #
# This script produces Figure 1 & Figure 2 in Section 7.2 of the paper.        #
################################################################################

#-- Clear the workspace and load the necessary libraries
rm(list = ls())
library(lpSolve) # for solving LP (e.g., Eq. (11) in the paper)
library(latex2exp) # for TeX expressions
library(ggplot2) # for plotting
library(dplyr) # for data manipulations
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

#-- Plot the pdf of the outer scenarios and the true conditional mean
#-- Quantiles, s_tau, outer scen. density, and true conditional mean
qts_plt <- seq(0.001, 0.999, 0.001) # quantiles for plotting
s_tau_plt <- qlnorm(qts_plt, log(s0) + meanlog_outer, sdlog_outer)
s_den_plt <- dlnorm(s_tau_plt, log(s0) + meanlog_outer, sdlog_outer)
mu_s_tau_plt <- sapply(
  s_tau_plt,
  function(x) bs_price_straddle(x, str_k, div, rf, vol, (big_t - tau))
) # true conditional mean

#-- data frame for plotting
df_plt <- data.frame(s_tau = s_tau_plt, SDen = s_den_plt, Val = mu_s_tau_plt)
#-- pdf of the outer scenarios (Figure 1a in the paper)
ggplot(data = df_plt, aes(x = s_tau_plt, y = s_den_plt)) +
  geom_area(fill = "deeppink3", alpha = 0.5) +
  labs(
    x = TeX("Projected stock prices $S_{\\tau}=\\theta$"),
    y = "Outer scenario density function",
    title = ""
  ) +
  theme(
    plot.title = element_text(size = 20, hjust = 0.5),
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 16),
    text = element_text(size = 16)
  )
ggsave("OuterDist.pdf",
  width = 6, height = 4,
  units = "in", dpi = 600
) # save the plot as a pdf

#-- mu(S_tau) vs. S_tau (Figure 1b in the paper)
ggplot(data = df_plt, aes(x = s_tau_plt, y = mu_s_tau_plt)) +
  geom_line(linewidth = 2, color = "steelblue3") +
  labs(
    x = TeX("Projected stock prices $S_{\\tau}=\\theta$"),
    y = TeX("Conditional mean $\\mu(\\theta)$"),
    title = ""
  ) +
  theme(
    plot.title = element_text(size = 20, hjust = 0.5),
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 16),
    text = element_text(size = 16)
  )
ggsave("ConditionalMean.pdf",
  width = 6, height = 4,
  units = "in", dpi = 600
) # save the plot as a pdf

#-- Experiment parameters
n_outer <- 1e3 # number of outer scenarios
n_inner <- 1e3 # number of inner replications (N0 for ESS calculation)
n_macro <- 1e4 # nubmer of macro replications

#-- Outer simulation (real world projection from time 0 to tau)
#-- Deterministic outer scenarios
qts <- (1:n_outer) / (n_outer + 1) # equally-spaced quantiles
s_tau <- qlnorm(qts, log(s0) + meanlog_outer, sdlog_outer)
# mean(s_tau) - s0*exp(mu*tau) # unit test, result should be close to 0
#-- sStochastic outer scenarios (not used in the paper)
# s_tau <- s0 * exp(meanlog_outer + rnorm(n_outer) * sdlog_outer)
# s_tau <- rlnorm(n_outer, log(s0) + meanlog_outer, sdlog_outer)

#-- closed-form calculations (the true conditional means)
val_true <- sapply(
  s_tau,
  function(x) bs_price_straddle(x, str_k, div, rf, vol, (big_t - tau))
)

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

psi_mat_nz <- psi_mat[, ind_nz] # submatrix of psi_mat for non-zero scenarios
gamma_nz <- sweep(psi_mat_nz, MARGIN = 2, n_j_opt_nested, "*")
gamma_nz <- gamma_nz / rowSums(gamma_nz) # optimal pooling weights

#-- print out the optimal design
print(sprintf("==== Optimal design"))
print(sprintf("---- %d non-zero scenarios: theta = %4.2f", n_nz, s_tau[ind_nz]))
print(sprintf("---- Simulation budget %d.", sim_budget_opt_nested))

#-- Optimal LR design simulation
# placeholder for estimated conditional means and runtime
val_opt_nested <- matrix(0, nrow = n_outer, ncol = n_macro)
time_opt_nested <- rep(0, 2) # first and second moment for runtime
# loop through macro replications
for (mac in 1:n_macro) {
  print(sprintf("Optimal nested sim, macro rep %d", mac))

  time_start <- Sys.time()
  # optimal nested simulation design with likelihood ratio
  val_opt_nested[, mac] <- straddle_opt_nested(
    s_tau, str_k, div, rf, vol, big_t, tau,
    s_tau[ind_nz], n_j_opt_nested, gamma_nz
  )
  time_end <- Sys.time()
  time_temp <- as.numeric(difftime(time_end, time_start, units = c("secs")))
  time_opt_nested[1] <- time_opt_nested[1] + time_temp
  time_opt_nested[2] <- time_opt_nested[2] + time_temp^2
}
time_opt_nested <- time_opt_nested / n_macro

#-- Standard Nested Simulation with Gamma^(1/3) inner reps and Gamma^(2/3)
#   outer scenarios. This method is referred to as SNS in the paper
#-- number of inner and outer replications
n_inner_std_nested <- ceiling(sim_budget_opt_nested^(1 / 3))
n_outer_std_nested <- ceiling(sim_budget_opt_nested^(2 / 3))
#-- equally-spaced quantiles as outer scenarios
qts_std_nested <- (1:n_outer_std_nested) / (n_outer_std_nested + 1)
s_tau_std_nested <- s0 * exp(meanlog_outer +
  qnorm(qts_std_nested) * sdlog_outer)

# placeholder for estimated conditional means and runtime
val_std_nested <- matrix(0, nrow = n_outer_std_nested, ncol = n_macro)
time_std_nested <- rep(0, 2) # first and second moment for runtime
for (mac in 1:n_macro) {
  print(sprintf("Stadard nested sim with optimal budget, macro rep %d", mac))

  time_start <- Sys.time()
  # standard nested simulation with similar budget and the optimal design
  val_std_nested[, mac] <- straddle_std_nested(
    s_tau_std_nested, str_k, div, rf, vol,
    big_t, tau, n_inner_std_nested
  )

  time_end <- Sys.time()
  time_temp <- as.numeric(difftime(time_end, time_start, units = c("secs")))
  time_std_nested[1] <- time_std_nested[1] + time_temp
  time_std_nested[2] <- time_std_nested[2] + time_temp^2
}
time_std_nested <- time_std_nested / n_macro

#-- Standard Nested Simulation with NInner inner reps per outer scenario.
#-- This method is referred to as SNS+ in the paper
#-- placeholder for estimated conditional means and runtime
val_std_plus <- matrix(0, nrow = n_outer, ncol = n_macro)
time_std_plus <- rep(0, 2) # first and second moment for runtime
# loop through macro replications
for (mac in 1:n_macro) {
  print(sprintf("Standard nestd sim, macro rep %d", mac))

  time_start <- Sys.time()
  # standard nested simulation
  val_std_plus[, mac] <- straddle_std_nested(
    s_tau, str_k, div, rf, vol,
    big_t, tau, n_inner
  )
  time_end <- Sys.time()
  time_temp <- as.numeric(difftime(time_end, time_start, units = c("secs")))
  time_std_plus[1] <- time_std_plus[1] + time_temp
  time_std_plus[2] <- time_std_plus[2] + time_temp^2
}
time_std_plus <- time_std_plus / n_macro

#-- Regression-based Nested Simulation
n_inner_reg <- 1 # number of inner replications at each design point
# number of design points
n_outer_reg <- ceiling(sim_budget_opt_nested / n_inner_reg)
# equally-spaced quantiles of s_tau as design points
qts_reg <- (1:n_outer_reg) / (n_outer_reg + 1)
s_tau_reg <- qlnorm(qts_reg, log(s0) + meanlog_outer, sdlog_outer)

# placeholder for estimated conditional means and runtime
val_reg_nested <- matrix(0, nrow = n_outer, ncol = n_macro)
time_reg_nested <- rep(0, 2) # first and second moment for runtime
# loop through macro replications
for (mac in 1:n_macro) {
  print(sprintf("Regression nested sim, macro rep %d", mac))

  time_start <- Sys.time()
  # regression-based nested simulation
  val_reg_nested[, mac] <- straddle_reg_nested(
    s_tau, str_k, div, rf, vol, big_t, tau,
    s0, s_tau_reg, n_inner_reg
  )
  time_end <- Sys.time()
  time_temp <- as.numeric(difftime(time_end, time_start, units = c("secs")))
  time_reg_nested[1] <- time_reg_nested[1] + time_temp
  time_reg_nested[2] <- time_reg_nested[2] + time_temp^2
}
time_reg_nested <- time_reg_nested / n_macro

#-- Table of runtime
time_table <- cbind(
  time_opt_nested, time_reg_nested,
  time_std_nested, time_std_plus
)
time_table[2, ] <- sqrt(time_table[2, ] - time_table[1, ]^2) / sqrt(n_macro)
rownames(time_table) <- c("mean", "std.err")
colnames(time_table) <- c(
  "Opt.Nested", "Reg.Nested",
  "Std.Nested", "Std.Nested+"
)
time_table

#-- Figure 2 in the paper
#-- Data preparation for plotting
quantile_fun <- function(df, probs) {
  return(t(apply(df, 1, function(x) quantile(x, probs))))
}

#-- Data preparation for plotting
df_std_plus <- data.frame(
  x = s_tau, true = val_true, type = "Std.Nested",
  mean = rowMeans(val_std_plus),
  quantile_fun(val_std_plus, c(0.025, 0.975))
)
df_opt_nested <- data.frame(
  x = s_tau, true = val_true, type = "Opt.Nested",
  mean = rowMeans(val_opt_nested),
  quantile_fun(val_opt_nested, c(0.025, 0.975))
)
df_reg_nested <- data.frame(
  x = s_tau, true = val_true, type = "Regression",
  mean = rowMeans(val_reg_nested),
  quantile_fun(val_reg_nested, c(0.025, 0.975))
)

df_plt2 <- rbind(df_opt_nested, df_std_plus, df_reg_nested)
colnames(df_plt2) <- c("Stau", "True", "Type", "Val", "LB", "UB")
df_plt2 <- mutate(df_plt2,
  Type = factor(Type,
    levels = c("Opt.Nested", "Std.Nested", "Regression")
  )
)

fig3 <- ggplot(df_plt2) +
  geom_ribbon(
    aes(
      ymin = LB, ymax = UB, x = Stau,
      colour = Type, linetype = Type
    ),
    linewidth = 1.5, alpha = 0
  ) +
  geom_line(linewidth = 1, aes(x = Stau, y = True), colour = "black") +
  scale_linetype_manual(
    values = c(1, 2, 3),
    labels = c(
      expression("Opt. Design"),
      expression("SNS+       "),
      expression("Regression ")
    )
  ) +
  scale_colour_manual(
    values = c("#d7191c", "#7fc97f", "black"),
    labels = c(
      expression("Opt. Design"),
      expression("SNS+       "),
      expression("Regression ")
    )
  ) +
  guides(
    colour = guide_legend(title = ""),
    linetype = guide_legend(title = "")
  ) +
  theme_bw() +
  theme(
    legend.position = c(0.2, 0.8),
    legend.direction = "vertical",
    legend.background = element_rect(fill = "transparent"),
    legend.key.size = unit(1, "cm"),
    legend.key.width = unit(2, "cm"),
    plot.title = element_text(hjust = 0.5),
    text = element_text(size = 18)
  ) +
  labs(
    title = "",
    x = TeX("outer scenarios $\\theta$"),
    y = TeX("conditional mean, $\\mu(\\theta)$")
  )
fig3
ggsave(fig3,
  file = "ERM_All_RedGreen.pdf",
  width = 7, height = 5, units = "in", dpi = 600
)

fig4 <- fig3 + xlim(80, 120) + ylim(30, 43)
fig4
ggsave(fig4,
  file = "ERM_Zoomed_RedGreen.pdf",
  width = 7, height = 5, units = "in", dpi = 600
)

filename <- sprintf(
  "SingleAssetERM_1_%dInner%dOuter%dMacro.RData",
  n_inner, n_outer, n_macro
)
save.image(filename)
