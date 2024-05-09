################################################################################
# ERM_Single_Asset_1.R                                                         #
# Single-Asset ERM example in Section 7.1 in the paper                         #
# Straddle option portfolio that consists of a long call and a long put        #
# Optimized nested simulation for enterprise risk management (ERM) example     #
# -- straddle option (long call + long put,the same strike for both options)   #
# This script loads simulation results directly then produces Figure 1 &       #
# Figure 2 in Section 7.2 of the paper.                                        #
################################################################################

rm(list = ls())
library(ggplot2) # for plotting
load("SingleAssetERM_Figure1_Figure2.RData")

#-- pdf of the outer scenarios (Figure 1a in the paper)
ggplot(data = df_plt, aes(x = s_tau, y = SDen)) +
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
ggplot(data = df_plt, aes(x = s_tau, y = Val)) +
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

#-- 95% confidence bands
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

#-- 95% confidence bands (zoomed)
fig4 <- fig3 + xlim(80, 120) + ylim(30, 43)
fig4
ggsave(fig4,
       file = "ERM_Zoomed_RedGreen.pdf",
       width = 7, height = 5, units = "in", dpi = 600
)
