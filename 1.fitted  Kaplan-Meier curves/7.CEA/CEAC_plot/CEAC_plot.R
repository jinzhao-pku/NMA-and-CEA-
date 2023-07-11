library(ggplot2)
library(ggsci)
library(dplyr)
library(tidyr)
library(reshape2)

rm(list = ls())
setwd("D:/project/NMA/CEA/CEAC_plot")

df <- read.csv("CEAC.csv",fileEncoding = "UTF-8-BOM")
df_long <- melt(df, id.vars = "WTP")

ggplot(df_long, aes(x = WTP, y = value, color = variable)) +
  geom_smooth(method="loess", se=FALSE) +
  scale_y_continuous(limits = c(0, 1), name = "possibility acceptance") +
  scale_x_continuous(name = "WTP threshold($/QALY)")+
  scale_color_jama() +
  labs(color = "group") +
  theme(legend.position = "bottom", panel.background = element_rect(fill = "white", color = "black"), axis.line = element_line(color = "black"))+
  geom_vline(xintercept = 38184, linetype = "dashed", color = "grey", size = 1) +
  annotate("text", x = 38184, y = 0.8, label = "WTP threshold = $38184/QALY")
