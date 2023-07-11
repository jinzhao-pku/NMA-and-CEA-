# Create a network plot

# Load library
library(netmeta)
library(tidyverse)
library(patchwork)
# Start with empty environment
rm(list=ls())

# Set the workign directory
setwd("C:/Users/28629/Documents/OS")

# OS AND PFS
data <- read.csv("data_os and pfs.csv")
data1 <- data%>%filter(class=="os")
data1$lhr<-log(data1$hr)
data1$se<-(log(data1$upper)-log(data1$lower))/(2*1.96)
data2 <- data%>%filter(class=="pfs")
data2$lhr<-log(data2$hr)
data2$se<-(log(data2$upper)-log(data2$lower))/(2*1.96)
a1 <- netmeta(data1$lhr, data1$se, treat1=data1$t1, treat2=data1$t2, studlab=data1$study, reference=1)
a2 <- netmeta(data2$lhr, data2$se, treat1=data2$t1, treat2=data2$t2, studlab=data2$study, reference=1)
# Treatment labels
lab <- c( "che","cam+che", "nivo+che", "sinti+che", "tori+che", "pemb+che","tisle+che")
 netgraph(a1, labels=lab, offset=0.02, plastic=F, col="#00A087FF", points=T,
         bg.points="#DC0000FF", col.point="#DC0000FF",number.of.studies = T, cex=2,rotate=0, 
         cex.number.of.studies =1.5,col.number.of.studies = "black",bg.number.of.studies = "white",
         cex.points=2*c(6,1,1,1,1,1,1))
title("OS",adj=0.5,cex.main=2.5)
export::graph2tif(file= "figure3-1",,width=10,height=10)
netgraph(a2, labels=lab, offset=0.02, plastic=F, col="#00A087FF", points=T,
               bg.points="#DC0000FF", col.point="#DC0000FF",number.of.studies = T, cex=2,rotate=-90, 
               cex.number.of.studies =1.5,col.number.of.studies = "black",bg.number.of.studies = "white",
               cex.points=2*c(6,1,1,1,1,1,1))
title("PFS",adj=0.5,cex.main=2.5)
export::graph2tif(file = "figure3-2",width=10,height=10)

meta1<-metafor::rma(yi=lhr,sei=se,data=data1, method="FE")
b1 <- metaviz::viz_sunset(meta1,contours = TRUE,method="FE")
meta2<-metafor::rma(yi=lhr,sei=se,data=data2, method="FE")
b2 <- metaviz::viz_sunset(meta2,contours = TRUE,method="FE")
b1+b2+plot_annotation(title= "sunset funnel plot for OS and PFS",tag_levels = "A")
export::graph2tif(file="figure3-3",width=20,height=10)