library(ggplot2)
library(ggsci)
library(dplyr)
library(tidyr)

rm(list = ls())
setwd("D:/project/NMA/CEA/point_plot")

data <- read.csv("point_plot.csv",fileEncoding = "UTF-8-BOM")

ggplot(data, aes(x = incremental_utility_cam, y = incremental_cost_cam)) +
  geom_point(aes(color = "cam"),alpha = 0.7) +
  geom_point(aes(x = incremental_utility_niv, y = incremental_cost_niv, color = "niv"),alpha = 0.7) +
  geom_point(aes(x = incremental_utility_sin, y = incremental_cost_sin, color = "sin"),alpha = 0.7) +
  geom_point(aes(x = incremental_utility_tor, y = incremental_cost_tor, color = "tor"),alpha = 0.7) +
  geom_point(aes(x = incremental_utility_pem, y = incremental_cost_pem, color = "pem"),alpha = 0.7) +
  geom_point(aes(x = incremental_utility_tis, y = incremental_cost_tis, color = "tis"),alpha = 0.7) +
  scale_color_jama() +
  labs(color = "group") +
  theme(legend.position = "bottom", panel.background = element_rect(fill = "white", color = "black"), axis.line = element_line(color = "black")) +
  labs(x = "incredential_utility", y = "incredential_cost")

ggplot(data, aes(x = incremental_utility_cam, y = incremental_cost_cam)) +
  geom_point(aes(color = "cam"),alpha = 0.7) +
  geom_point(aes(x = incremental_utility_sin, y = incremental_cost_sin, color = "sin"),alpha = 0.7) +
  geom_point(aes(x = incremental_utility_tor, y = incremental_cost_tor, color = "tor"),alpha = 0.7) +
  geom_point(aes(x = incremental_utility_tis, y = incremental_cost_tis, color = "tis"),alpha = 0.7) +
  scale_color_jama() +
  labs(color = "group") +
  theme(legend.position = "bottom", panel.background = element_rect(fill = "white", color = "black"), axis.line = element_line(color = "black")) +
  labs(x = "incredential_utility", y = "incredential_cost")
