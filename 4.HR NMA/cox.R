install.packages(c("remotes", "knitr"))
remotes::install_github("audrey-b/BUGSnet@v1.1.0", upgrade = TRUE,
                        build_vignettes = TRUE, dependencies = TRUE)
install.packages("ggpattern")

library(BUGSnet)
library(ggplot2)
library(ggsci)
library(readr)
library(dplyr)
library(tidyr)
library(export)

#OS----
rm(list = ls())
setwd("D:/project/NMA/cox/os")
data <- read_csv("cox_os.csv")
data <- data.prep(arm.data = data,
                  varname.t = "Treatment", 
                  varname.s = "Study")

#构建contrast_based模型
fixed_effects_model <- nma.model.contrast(data_contrast = data,
                                           differences = "logHR",
                                           se.diffs = "std.err",
                                           reference = "che",
                                           type = "consistency",
                                           effects = "fixed",
                                           scale = "Log-odds Ratio")
set.seed(20232023)
fixed_effects_results <- nma.run(fixed_effects_model, 
                                 n.adapt=1000,
                                 n.burnin=1000,
                                 n.iter=10000)
#NMA分析的效应模型评估，即评估模型的拟合和识别潜在异常值，并输出杠杆图
nma.fit(fixed_effects_results)
export::graph2tif(file="杠杆图",dpi=300,width=8,height=6)

#trace and diag plot
nma.diag(
  fixed_effects_results,
  trace = TRUE,
  gelman.rubin = TRUE,
  geweke = TRUE,
  params = "all",
  thin = 1,
  ncol = 1,
  nrow = 3,
  plot_prompt = TRUE,
  geweke_frac1 = 0.1,
  geweke_frac2 = 0.5
)
export::graph2tif(file="trace and density 2",dpi=300,width=8,height=6)

#治疗排名结果图
sucra_out<- nma.rank(fixed_effects_results, 
                     largerbetter=FALSE)
r<- sucra_out$rankogram+theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1))
r
export::graph2tif(file="治疗效果叠图",dpi=300,width=8,height=6)

sucra_out$sucraplot +
  scale_color_jama() +
  ggtitle("OS")
export::graph2tif(file="治疗效果折线图",dpi=300,width=8,height=6)

#排名表热图
league.out<- nma.league(fixed_effects_results,
                        central.tdcy="median",
                        order =sucra_out$order,
                        log.scale =FALSE,
                        low.colour ="springgreen4",
                        mid.colour = "white", 
                        high.colour ="red", 
                        digits = 2)
m = league.out$heatplot

m[["plot_env"]][["midpoint"]] = 1
m[["data"]][["ct.stat"]] = exp(m[["data"]][["ct.stat"]])
m[["data"]][["lci"]] = exp(m[["data"]][["lci"]])
m[["data"]][["uci"]] = exp(m[["data"]][["uci"]])
m[["layers"]][[2]][["mapping"]][["label"]][[2]][[3]][[3]][[2]]=""
m[["layers"]][[2]][["mapping"]][["label"]][[2]][[3]][[3]][[4]]=""
m +
  scale_fill_gradient2(low="red", mid="white", high="springgreen4", midpoint=1)
export::graph2tif(file="联赛图",dpi=300,width=8,height=6)


#森林图
f<- nma.forest(fixed_effects_results,
               central.tdcy="median",
               log.scale =FALSE,
               comparator ='che')
f
#手动画图
# 创建数据框
df <- data.frame(
  scheme = f[["data"]][["trt"]],
  hr = exp(f[["data"]][["mean"]]),
  lower = exp(f[["data"]][["lci"]]),
  upper = exp(f[["data"]][["uci"]])
)
# 创建一个新的数据框
rect_df <- data.frame(scheme = df$scheme,
                      ymin = seq(0.5, 5.5, by = 1),
                      ymax = seq(1.5, 6.5, by = 1),
                      fill = c("grey", "white"))
# 绘制森林图
df$hr <- round(df$hr, 2)
df$lower <- round(df$lower, 2)
df$upper <- round(df$upper, 2)
df$text <- paste0("HR = ", df$hr, "\n95% CI = (", df$lower, ", ", df$upper, ")")
ggplot(df) +
  geom_pointrange(aes(x = hr,
                      y = scheme,
                      xmin = lower,
                      xmax = upper),
                  color = "black") +
  geom_vline(xintercept = 1,
             linetype = "dashed",
             color = "black") +
  xlab("HR (95% CI) compared with chemotherapy") +
  ylab("Scheme") +
  ggtitle("OS")+
  geom_text(aes(x = 0.85,
                y = scheme,
                label = text),
            hjust = 0)+
  scale_color_jama()+
  theme_bw()
export::graph2tif(file="森林图",dpi=300,width=8,height=6)

#PFS----
rm(list = ls())
setwd("D:/project/NMA/cox/pfs")
data <- read_csv("cox_pfs.csv")
data <- data.prep(arm.data = data,
                  varname.t = "Treatment", 
                  varname.s = "Study")

#构建contrast_based模型
fixed_effects_model <- nma.model.contrast(data_contrast = data,
                                          differences = "logHR",
                                          se.diffs = "std.err",
                                          reference = "che",
                                          type = "consistency",
                                          effects = "fixed",
                                          scale = "Log-odds Ratio")
set.seed(20232023)
fixed_effects_results <- nma.run(fixed_effects_model, 
                                 n.adapt=1000,
                                 n.burnin=1000,
                                 n.iter=10000)
#NMA分析的效应模型评估，即评估模型的拟合和识别潜在异常值，并输出杠杆图
nma.fit(fixed_effects_results)
export::graph2tif(file="杠杆图",dpi=300,width=8,height=6)

#trace and diag plot
nma.diag(
  fixed_effects_results,
  trace = TRUE,
  gelman.rubin = TRUE,
  geweke = TRUE,
  params = "all",
  thin = 1,
  ncol = 1,
  nrow = 3,
  plot_prompt = TRUE,
  geweke_frac1 = 0.1,
  geweke_frac2 = 0.5
)
export::graph2tif(file="trace and density 2",dpi=300,width=8,height=6)


#治疗排名结果图
sucra_out<- nma.rank(fixed_effects_results, 
                     largerbetter=FALSE)
r<- sucra_out$rankogram+theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1))
r
export::graph2tif(file="治疗效果叠图",dpi=300,width=8,height=6)

#折线图
sucra_out$sucraplot+
  scale_color_jama()+
  ggtitle("PFS")
export::graph2tif(file="治疗效果折线图",dpi=300,width=8,height=6)

#排名表热图
league.out<- nma.league(fixed_effects_results,
                        central.tdcy="median",
                        order =sucra_out$order,
                        log.scale =FALSE,
                        low.colour ="springgreen4",
                        mid.colour = "white", 
                        high.colour ="red", 
                        digits = 2)
m = league.out$heatplot

m[["plot_env"]][["midpoint"]] = 1
m[["data"]][["ct.stat"]] = exp(m[["data"]][["ct.stat"]])
m[["data"]][["lci"]] = exp(m[["data"]][["lci"]])
m[["data"]][["uci"]] = exp(m[["data"]][["uci"]])
m[["layers"]][[2]][["mapping"]][["label"]][[2]][[3]][[3]][[2]]=""
m[["layers"]][[2]][["mapping"]][["label"]][[2]][[3]][[3]][[4]]=""
m +
  scale_fill_gradient2(low="red", mid="white", high="springgreen4", midpoint=1)
export::graph2tif(file="联赛图",dpi=300,width=8,height=6)

#森林图
f<- nma.forest(fixed_effects_results,
               central.tdcy="median",
               log.scale =FALSE,
               comparator ='che')
f
#手动画图
# 创建数据框
df <- data.frame(
  scheme = f[["data"]][["trt"]],
  hr = exp(f[["data"]][["mean"]]),
  lower = exp(f[["data"]][["lci"]]),
  upper = exp(f[["data"]][["uci"]])
)
# 创建一个新的数据框
rect_df <- data.frame(scheme = df$scheme,
                      ymin = seq(0.5, 5.5, by = 1),
                      ymax = seq(1.5, 6.5, by = 1),
                      fill = c("grey", "white"))
# 绘制森林图
df$hr <- round(df$hr, 2)
df$lower <- round(df$lower, 2)
df$upper <- round(df$upper, 2)
df$text <- paste0("HR = ", df$hr, "\n95% CI = (", df$lower, ", ", df$upper, ")")
ggplot(df) +
  geom_pointrange(aes(x = hr,
                      y = scheme,
                      xmin = lower,
                      xmax = upper),
                  color = "black") +
  geom_vline(xintercept = 1,
             linetype = "dashed",
             color = "black") +
  xlab("HR (95% CI) compared with chemotherapy") +
  ylab("Scheme") +
  ggtitle("PFS")+
  geom_text(aes(x = 0.85,
                y = scheme,
                label = text),
            hjust = 0)+
  scale_color_jama()+
  theme_bw()
export::graph2tif(file="森林图",dpi=300,width=8,height=6)

#OS_sub----
rm(list = ls())
setwd("D:/project/NMA/cox/os_sub")
data <- read_csv("cox_os_sub.csv")
data <- data.prep(arm.data = data,
                  varname.t = "Treatment", 
                  varname.s = "Study")

#构建contrast_based模型
fixed_effects_model <- nma.model.contrast(data_contrast = data,
                                          differences = "logHR",
                                          se.diffs = "std.err",
                                          reference = "che",
                                          type = "consistency",
                                          effects = "fixed",
                                          scale = "Log-odds Ratio")
set.seed(20232023)
fixed_effects_results <- nma.run(fixed_effects_model, 
                                 n.adapt=1000,
                                 n.burnin=1000,
                                 n.iter=10000)
#NMA分析的效应模型评估，即评估模型的拟合和识别潜在异常值，并输出杠杆图
nma.fit(fixed_effects_results)
export::graph2tif(file="杠杆图",dpi=300,width=8,height=6)

#trace and diag plot
nma.diag(
  fixed_effects_results,
  trace = TRUE,
  gelman.rubin = TRUE,
  geweke = TRUE,
  params = "all",
  thin = 1,
  ncol = 1,
  nrow = 3,
  plot_prompt = TRUE,
  geweke_frac1 = 0.1,
  geweke_frac2 = 0.5
)
export::graph2tif(file="trace and density 2",dpi=300,width=8,height=6)

#治疗排名结果图
sucra_out<- nma.rank(fixed_effects_results, 
                     largerbetter=FALSE)
r<- sucra_out$rankogram+theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1))
r
export::graph2tif(file="治疗效果叠图",dpi=300,width=8,height=6)
#折线图
sucra_out$sucraplot+
  scale_color_jama()+
  ggtitle("OS(CPS≥10)")
export::graph2tif(file="治疗效果折线图",dpi=300,width=8,height=6)

#排名表热图
league.out<- nma.league(fixed_effects_results,
                        central.tdcy="median",
                        order =sucra_out$order,
                        log.scale =FALSE,
                        low.colour ="springgreen4",
                        mid.colour = "white", 
                        high.colour ="red", 
                        digits = 2)
m = league.out$heatplot

m[["plot_env"]][["midpoint"]] = 1
m[["data"]][["ct.stat"]] = exp(m[["data"]][["ct.stat"]])
m[["data"]][["lci"]] = exp(m[["data"]][["lci"]])
m[["data"]][["uci"]] = exp(m[["data"]][["uci"]])
m[["layers"]][[2]][["mapping"]][["label"]][[2]][[3]][[3]][[2]]=""
m[["layers"]][[2]][["mapping"]][["label"]][[2]][[3]][[3]][[4]]=""
m +
  scale_fill_gradient2(low="red", mid="white", high="springgreen4", midpoint=1)
export::graph2tif(file="联赛图",dpi=300,width=8,height=6)

#森林图
f<- nma.forest(fixed_effects_results,
               central.tdcy="median",
               log.scale =FALSE,
               comparator ='che')
f
#手动画图
# 创建数据框
df <- data.frame(
  scheme = f[["data"]][["trt"]],
  hr = exp(f[["data"]][["mean"]]),
  lower = exp(f[["data"]][["lci"]]),
  upper = exp(f[["data"]][["uci"]])
)
# 创建一个新的数据框
rect_df <- data.frame(scheme = df$scheme,
                      ymin = seq(0.5, 5.5, by = 1),
                      ymax = seq(1.5, 6.5, by = 1),
                      fill = c("grey", "white"))
# 绘制森林图
df$hr <- round(df$hr, 2)
df$lower <- round(df$lower, 2)
df$upper <- round(df$upper, 2)
df$text <- paste0("HR = ", df$hr, "\n95% CI = (", df$lower, ", ", df$upper, ")")
ggplot(df) +
  geom_pointrange(aes(x = hr,
                      y = scheme,
                      xmin = lower,
                      xmax = upper),
                  color = "black") +
  geom_vline(xintercept = 1,
             linetype = "dashed",
             color = "black") +
  xlab("HR (95% CI) compared with chemotherapy") +
  ylab("Scheme") +
  ggtitle("OS (CPS≥10)")+
  geom_text(aes(x = 0.85,
                y = scheme,
                label = text),
            hjust = 0)+
  scale_color_jama()+
  theme_bw()
export::graph2tif(file="森林图",dpi=300,width=8,height=6)

#PFS_sub----
rm(list = ls())
setwd("D:/project/NMA/cox/pfs_sub")
data <- read_csv("cox_pfs_sub.csv")
data <- data.prep(arm.data = data,
                  varname.t = "Treatment", 
                  varname.s = "Study")

#构建contrast_based模型
fixed_effects_model <- nma.model.contrast(data_contrast = data,
                                          differences = "logHR",
                                          se.diffs = "std.err",
                                          reference = "che",
                                          type = "consistency",
                                          effects = "fixed",
                                          scale = "Log-odds Ratio")
set.seed(20232023)
fixed_effects_results <- nma.run(fixed_effects_model, 
                                 n.adapt=1000,
                                 n.burnin=1000,
                                 n.iter=10000)
#NMA分析的效应模型评估，即评估模型的拟合和识别潜在异常值，并输出杠杆图
nma.fit(fixed_effects_results)
export::graph2tif(file="杠杆图",dpi=300,width=8,height=6)

#trace and diag plot
nma.diag(
  fixed_effects_results,
  trace = TRUE,
  gelman.rubin = TRUE,
  geweke = TRUE,
  params = "all",
  thin = 1,
  ncol = 1,
  nrow = 3,
  plot_prompt = TRUE,
  geweke_frac1 = 0.1,
  geweke_frac2 = 0.5
)
export::graph2tif(file="trace and density 2",dpi=300,width=8,height=6)

#治疗排名结果图
sucra_out<- nma.rank(fixed_effects_results, 
                     largerbetter=FALSE)
r<- sucra_out$rankogram+theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1))
r
export::graph2tif(file="治疗效果叠图",dpi=300,width=8,height=6)
#折线图
sucra_out$sucraplot+
  scale_color_jama()+
  ggtitle("PFS(CPS≥10)")
export::graph2tif(file="治疗效果折线图",dpi=300,width=8,height=6)

#排名表热图
league.out<- nma.league(fixed_effects_results,
                        central.tdcy="median",
                        order =sucra_out$order,
                        log.scale =FALSE,
                        low.colour ="springgreen4",
                        mid.colour = "white", 
                        high.colour ="red", 
                        digits = 2)
m = league.out$heatplot

m[["plot_env"]][["midpoint"]] = 1
m[["data"]][["ct.stat"]] = exp(m[["data"]][["ct.stat"]])
m[["data"]][["lci"]] = exp(m[["data"]][["lci"]])
m[["data"]][["uci"]] = exp(m[["data"]][["uci"]])
m[["layers"]][[2]][["mapping"]][["label"]][[2]][[3]][[3]][[2]]=""
m[["layers"]][[2]][["mapping"]][["label"]][[2]][[3]][[3]][[4]]=""
m +
  scale_fill_gradient2(low="red", mid="white", high="springgreen4", midpoint=1)
export::graph2tif(file="联赛图",dpi=300,width=8,height=6)

#森林图
f<- nma.forest(fixed_effects_results,
               central.tdcy="median",
               log.scale =FALSE,
               comparator ='che')
f
#手动画图
# 创建数据框
df <- data.frame(
  scheme = f[["data"]][["trt"]],
  hr = exp(f[["data"]][["mean"]]),
  lower = exp(f[["data"]][["lci"]]),
  upper = exp(f[["data"]][["uci"]])
)
# 创建一个新的数据框
rect_df <- data.frame(scheme = df$scheme,
                      ymin = seq(0.5, 3.5, by = 1),
                      ymax = seq(1.5, 4.5, by = 1),
                      fill = c("grey", "white"))
# 绘制森林图
df$hr <- round(df$hr, 2)
df$lower <- round(df$lower, 2)
df$upper <- round(df$upper, 2)
df$text <- paste0("HR = ", df$hr, "\n95% CI = (", df$lower, ", ", df$upper, ")")
ggplot(df) +
  geom_pointrange(aes(x = hr,
                      y = scheme,
                      xmin = lower,
                      xmax = upper),
                  color = "black") +
  geom_vline(xintercept = 1,
             linetype = "dashed",
             color = "black") +
  xlab("HR (95% CI) compared with chemotherapy") +
  ylab("Scheme") +
  ggtitle("PFS (CPS≥10)")+
  geom_text(aes(x = 0.85,
                y = scheme,
                label = text),
            hjust = 0)+
  scale_color_jama()+
  theme_bw()
export::graph2tif(file="森林图",dpi=300,width=8,height=6)
