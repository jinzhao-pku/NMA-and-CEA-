# Test for non-PH

# Load libraries
library(survival)
library(broom)
library(metafor)
library(survminer)
library(patchwork)


#pfs----
# Start with an empty environment
rm(list=ls())
setwd("D:/project/NMA/extra_check")
data <- read.csv("RP_PFS.csv")
# Set the working directory
path <- "D:/project/NMA/extra_check/plot/pfs"
setwd(path)
# Import the data and split into data frames for each trial
name<-c("ESCORT-1st.tiff","cm648.tiff","ju06.tiff","kn590.tiff","ra306.tiff","or15.tiff")

#S test
for (i in 1:6) {
  df<-data[data$studycode==i,]
  cox <- coxph(formula = Surv(time, event) ~ arm, data=df)
  test<-cox.zph(cox)
  tiff(name[i],width = 600,height=400)
  fig<-ggcoxzph(test)
  print(fig)
  dev.off()
}

#拼图
library(ggpubr)
library(stringr)

fig=list()
for (i in 1:6) {
  df<-data[data$studycode==i,]
  cox <- coxph(formula = Surv(time, event) ~ arm, data=df)
  test<-cox.zph(cox)
  fig[i]<-ggcoxzph(test)
  print(fig[i])
  dev.off()
}

p=list()
for(i in 1:6){
  p[[i]]=fig[[i]]$plot_env$gplot
  a=str_extract_all(p[[i]]$labels$title,"(?<= )p: 0.[0-9]+")
  p[[i]] = p[[i]]+annotate("text", x = layer_scales(p[[i]])$x$range$range[2]/2, y = 13, label = a)
  p[[i]]$labels$title = NULL
  p[[i]] = p[[i]] + ylim(c(-25,25))
}

mix = ggarrange(p[[1]],p[[2]],p[[3]],p[[4]],p[[5]],p[[6]],ncol = 3, nrow = 2,labels = c("A","B","C","D","E","F"))
mix$labels$title = "Schoenfeld Individual Test"
mix

#log cumulative hazard plot
fig1 = list()
for (i in 1:6) {
  df<-data[data$studycode==i,]
  cox = survfit(Surv(time, event) ~ arm, data=df)
  coxplot = ggsurvplot(cox, data=df,
                       fun = "cumhaz")
  m= coxplot$plot 
  pg <- ggplot_build(m)
  x = log(pg$data[[1]]$x[c(-1,-2)],exp(1))
  y = -log(pg$data[[1]]$y[c(-1,-2)],exp(1))
  k = pg$data[[1]]$group[c(-1,-2)]
  g = if_else(k=="2","PD-1+che","che")
  df1 = tibble(x,y,g)
  df1$g = as.factor(df1$g)
  
  fig1[[i]] = ggplot(data=df1, aes(x=x, y=y,group=g,color = g)) +
    geom_line()+
    scale_x_continuous(name="ln(Analysis time)")+
    scale_y_continuous(name="-ln(-ln(Survival probability))")+
    scale_color_manual("treatment", values= c("#006000","#FF4500"))+
    theme_bw()
}
fig1[[1]]


mix1 = ggarrange(fig1[[1]],fig1[[2]],fig1[[3]],fig1[[4]],fig1[[5]],fig1[[6]],ncol = 3, nrow = 2,labels = c("A","B","C","D","E","F"))
mix1$labels$title = "Log Cumulative Hazard Plot"
mix1

#os----
# Start with an empty environment
rm(list=ls())
setwd("D:/project/NMA/extra_check")
data <- read.csv("RP_OS.csv")
# Set the working directory
path <- "D:/project/NMA/extra_check/plot/os"
setwd(path)
# Import the data and split into data frames for each trial
name<-c("ESCORT-1st.tiff","cm648.tiff","ju06.tiff","kn590.tiff","ra306.tiff","or15.tiff")

#S test
for (i in 1:6) {
  df<-data[data$studycode==i,]
  cox <- coxph(formula = Surv(time, event) ~ arm, data=df)
  test<-cox.zph(cox)
  tiff(name[i],width = 600,height=400)
  fig<-ggcoxzph(test)
  print(fig)
  dev.off()
}

#拼图
library(ggpubr)
library(stringr)

fig=list()
for (i in 1:6) {
  df<-data[data$studycode==i,]
  cox <- coxph(formula = Surv(time, event) ~ arm, data=df)
  test<-cox.zph(cox)
  fig[i]<-ggcoxzph(test)
  print(fig[i])
  dev.off()
}

p=list()
for(i in 1:6){
  p[[i]]=fig[[i]]$plot_env$gplot
  a=str_extract_all(p[[i]]$labels$title,"(?<= )p: 0.[0-9]+")
  p[[i]] = p[[i]]+annotate("text", x = layer_scales(p[[i]])$x$range$range[2]/2, y = 13, label = a)
  p[[i]]$labels$title = NULL
  p[[i]] = p[[i]] + ylim(c(-25,25))
}

mix = ggarrange(p[[1]],p[[2]],p[[3]],p[[4]],p[[5]],p[[6]],ncol = 3, nrow = 2,labels = c("A","B","C","D","E","F"))
mix$labels$title = "Schoenfeld Individual Test"
mix
#log cumulative hazard plot
fig1 = list()
for (i in 1:6) {
  df<-data[data$studycode==i,]
  cox = survfit(Surv(time, event) ~ arm, data=df)
  coxplot = ggsurvplot(cox, data=df,
                       fun = "cumhaz")
  m= coxplot$plot 
  pg <- ggplot_build(m)
  x = log(pg$data[[1]]$x[c(-1,-2)],exp(1))
  y = -log(pg$data[[1]]$y[c(-1,-2)],exp(1))
  k = pg$data[[1]]$group[c(-1,-2)]
  g = if_else(k=="2","PD-1+che","che")
  df1 = tibble(x,y,g)
  df1$g = as.factor(df1$g)
  
  fig1[[i]] = ggplot(data=df1, aes(x=x, y=y,group=g,color = g)) +
    geom_line()+
    scale_x_continuous(name="ln(Analysis time)")+
    scale_y_continuous(name="-ln(-ln(Survival probability))")+
    scale_color_manual("treatment", values= c("#006000","#FF4500"))+
    theme_bw()
  }
fig1[[1]]

mix1 = ggarrange(fig1[[1]],fig1[[2]],fig1[[3]],fig1[[4]],fig1[[5]],fig1[[6]],ncol = 3, nrow = 2,labels = c("A","B","C","D","E","F"))
mix1$labels$title = "Log Cumulative Hazard Plot"
mix1
