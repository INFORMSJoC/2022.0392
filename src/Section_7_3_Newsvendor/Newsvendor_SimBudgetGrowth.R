################################################################################
# Newsvendor_SimBudgetGrowth.R                                                 #
# Newsvedor example in Section 7.3 in the paper                                #
# Single-stage newsvendor problem with 10 products                             #
# Product demands X_ell, ell = 1,...,10 follow ind. Poisson distributions      #
# Total profit g(X) = sum (p_ell * min{X_ell, k_ell} - c_ell k_ell)            #
# This script produces Figure 3b in Section 7.3                                #
################################################################################

#-- Clear the workspace and load the necessary libraries
rm(list = ls())
source("Newsvendor_Helpers.R")
library(lpSolve) # for solving LP (e.g., Eq. (11) in the paper)
library(ggplot2) # for plotting

#-- Newsvendor problem setting
n_products <- 10 # number of products
# gamma_params: shape and rate parameters of the posterior gamma distribution
# QUESTION: Where do these parameters come from?
gamma_params <- cbind(
  5 * (11:20) + 0.001,
  0.001 + c(343, 396, 547, 619, 736, 939, 983, 1152, 1361, 1544)
)

#-- Simulation design settings
n_outer <- (1:10) * 1e3 # number of outer scenarios
n_settings <- length(n_outer) # number of different settings/outer scenarios
n_macro <- 100 # number of macro replications

#-- placeholder for optimal simulation in different settings
sim_budget <- matrix(0, n_macro, n_settings)

#-- main simulation experiment
for (b in 1:n_settings) {
  # set the same number of outer and inner replications
  n_outer_b <- n_outer[b]
  n_inner_b <- n_outer[b]

  for (mac in 1:n_macro) {
    print(sprintf(
      "Running the %d macro replication with %d outer scenarios.",
      mac, n_outer_b
    ))

    #-- simulate outer scenarios
    lambda_outer <- matrix(0, n_outer_b, n_products)
    for (i in 1:n_products) {
      lambda_outer[, i] <- rgamma(n_outer_b,
        rate = gamma_params[i, 1], shape = gamma_params[i, 2]
      )
    }

    #-- compute the optimal simulation budget given the outer scenarios
    sim_budget[mac, b] <- newsvendor_opt_design(
      lambda_outer,
      n_inner_b
    )$sim_budget
  }
}

# average simulation budget over macro replications
avg_sim_budget <- apply(sim_budget, 2, mean)

filename <- sprintf("Newsvendor_SimBudgetGrowth.RData")
save.image(filename)

#-- plot the results
df_sim_budget <- data.frame(opt_sim_budget = avg_sim_budget, n_outer = n_outer)

ggplot(df_sim_budget, aes(x = n_outer, y = opt_sim_budget)) +
  geom_point(size = 5) +
  ylab("") +
  xlab("Number of Outer Scenarios, M") +
  ggtitle("Minimized Simulation Budget with M=N") +
  theme(
    plot.title = element_text(size = 14, face = "bold.italic", hjust = 0.5),
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 14)
  )
ggsave(
  file = "SimBudgetGrowth.pdf",
  width = 7, height = 5, units = "in", dpi = 600
)
