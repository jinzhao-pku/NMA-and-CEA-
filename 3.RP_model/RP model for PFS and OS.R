# Royston-Parmar model in R - fixed effect

# Empty environment
rm(list = ls()) 

# Load libraries
library(survival)
library(foreach)
library(R2WinBUGS)
library(here)
library(tidyverse)
library(ggh4x)


####--------------------------------------------------------------------
rm(list=ls())
data <- read.csv("RP_PFS.csv")
path <- getwd()
setwd(path)
bugs.directory <- "C:/Users/28629/Documents/WinBUGS14"
# Set the location for WinBUGS
num.sims <- 5000
burn.in <- 5000
# Number of studies
num_studies <- length(unique(data$studycode))
# Sort data into trialid order
data <- data[order(data$studycode),]
# Number of patients in each study
# Note: if data is re-ordered will now to re-order this vector
num_patients <- c(596,645,514,749,649,659)####
# Number of treatments
nt <- length(unique(data$txcode))
# Create treatment indicator variables - caution needed to ensure consistency equations hold
data$trt2 <- 0
data$trt2[data$txcode==2] <- 1
data$trt3 <- 0
data$trt3[data$txcode==3] <- 1
data$trt4 <- 0
data$trt4[data$txcode==4] <- 1
data$trt5 <- 0
data$trt5[data$txcode==5] <- 1
data$trt6 <- 0
data$trt6[data$txcode==6] <- 1
data$trt7 <- 0
data$trt7[data$txcode==7] <- 1
# Add  trialid variable - easier than changing trialid to studyCode throughout this file!
data$trialid <- data$studycode
# Put eventtime onto the ln scale
data$eventtime <- data$time
data$lnt <- log(data$eventtime)
# Create trt*lnt variables
data$trt2lnt <- data$trt2*data$lnt
data$trt3lnt <- data$trt3*data$lnt
data$trt4lnt <- data$trt4*data$lnt
data$trt5lnt <- data$trt5*data$lnt
data$trt6lnt <- data$trt6*data$lnt
data$trt7lnt <- data$trt7*data$lnt
# Total number of patients
pts <- nrow(data)
# Set location of knots for each trial - 33rd and 67th percentiles of uncensored survival times
knots <- data.frame(trialid=NA, k1=NA, k2=NA, k3=NA, k4=NA)
foreach(i=1:num_studies) %do% {
  k <- quantile(data$lnt[data$event==1 & data$trialid==i], c(0, 0.33, 0.67, 1))
  knots <- rbind(knots, data.frame(trialid=i, k1=k[1], k2=k[2], k3=k[3], k4=k[4]))
}
# Remove the empty top row
knots <- knots[-1,]
# Add knot values to data
foreach(i=1:num_studies) %do% {
  data$k1[data$trialid==i] <- knots$k1[i]
  data$k2[data$trialid==i] <- knots$k2[i]
  data$k3[data$trialid==i] <- knots$k3[i]
  data$k4[data$trialid==i] <- knots$k4[i]
}  
# Orthogonalisation
data$t1 <- data$lnt-data$k2
data$t2 <- data$lnt-data$k1
data$t3 <- data$lnt-data$k4
data$t4 <- data$lnt-data$k3
foreach(i=1:pts) %do% {
  data$a[i] <- max(c(data$t1[i], 0))
  data$b[i] <- max(c(data$t2[i], 0))
  data$c[i] <- max(c(data$t3[i], 0))
  data$e[i] <- max(c(data$t4[i], 0))
  
  data$phi1[i] <- (data$k4[i]-data$k2[i])/(data$k4[i]-data$k1[i])
  data$phi2[i] <- (data$k4[i]-data$k3[i])/(data$k4[i]-data$k1[i])
  
  data$v0[i] <- data$lnt[i]
  data$v1[i] <- data$a[i]^3 - (data$phi1[i]*(data$b[i]^3)) - ((1-data$phi1[i])*(data$c[i]^3))
  data$v2[i] <- data$e[i]^3 - (data$phi2[i]*(data$b[i]^3)) - ((1-data$phi2[i])*(data$c[i]^3))
}
##--------------##
df1 <- data.frame(trialid=NA, mean=NA, sd=NA)
foreach(i=1:num_studies) %do% {
  mv0 <- mean(data$v0[data$trialid==i])
  sdv0 <- sd(data$v0[data$trialid==i])
  df1 <- rbind(df1, data.frame(trialid=i, mean=mv0, sd=sdv0))
}
df1 <- df1[-1,]
# Add values back to data
foreach(i=1:num_studies) %do% {
  data$mv0[data$trialid==i] <- df1$mean[i]
  data$sdv0[data$trialid==i] <- df1$sd[i]
} 
foreach(i=1:pts) %do% {
  data$u0[i] <- (data$v0[i]-data$mv0[i])/data$sdv0[i]
  data$f[i] <- data$v1[i]*data$u0[i]
  data$g[i] <- data$u0[i]*data$u0[i]
}
df2 <- data.frame(trialid=NA, v1u0=NA, u0u0=NA)
foreach(i=1:num_studies) %do% {
  v1u0 <- sum(data$f[data$trialid==i])
  u0u0 <- sum(data$g[data$trialid==i])
  df2 <- rbind(df2, data.frame(trialid=i, v1u0=v1u0, u0u0=u0u0))
}
df2 <- df2[-1,]
# Add values back to data
foreach(i=1:num_studies) %do% {
  data$v1u0[data$trialid==i] <- df2$v1u0[i]
  data$u0u0[data$trialid==i] <- df2$u0u0[i]
}
foreach(i=1:pts) %do% {
  data$u1[i] <- data$v1[i] - ((data$v1u0[i]/data$u0u0[i])*data$u0[i])
}
##--------------##
df3 <- data.frame(trialid=NA, mean=NA, sd=NA)
foreach(i=1:num_studies) %do% {
  mu1 <- mean(data$u1[data$trialid==i])
  sdu1 <- sd(data$u1[data$trialid==i])
  df3 <- rbind(df3, data.frame(trialid=i, mean=mu1, sd=sdu1))
}
df3 <- df3[-1,]
# Add values back to data
foreach(i=1:num_studies) %do% {
  data$mu1[data$trialid==i] <- df3$mean[i]
  data$sdu1[data$trialid==i] <- df3$sd[i]
} 
foreach(i=1:pts) %do% {
  data$u1n[i] <- (data$u1[i]-data$mu1[i])/data$sdu1[i]
  data$h[i] <- data$v2[i]*data$u0[i]
  data$aa[i] <- data$v2[i]*data$u1n[i]
  data$j[i] <- data$u1n[i]*data$u1n[i]
}
df4 <- data.frame(trialid=NA, v2u0=NA, v2u1n=NA, u1nu1n=NA)
foreach(i=1:num_studies) %do% {
  v2u0 <- sum(data$h[data$trialid==i])
  v2u1n <- sum(data$aa[data$trialid==i])
  u1nu1n <- sum(data$j[data$trialid==i])
  df4 <- rbind(df4, data.frame(trialid=i, v2u0=v2u0, v2u1n=v2u1n, u1nu1n=u1nu1n))
}
df4 <- df4[-1,]

# Add values back to data
foreach(i=1:num_studies) %do% {
  data$v2u0[data$trialid==i] <- df4$v2u0[i]
  data$v2u1n[data$trialid==i] <- df4$v2u1n[i]
  data$u1nu1n[data$trialid==i] <- df4$u1nu1n[i]
}
foreach(i=1:pts) %do% {
  data$u2[i] <- data$v2[i] - ((data$v2u0[i]/data$u0u0[i])*data$u0[i]) - ((data$v2u1n[i]/data$u1nu1n[i])*data$u1n[i])
}
##--------------##
df5 <- data.frame(trialid=NA, mean=NA, sd=NA)
foreach(i=1:num_studies) %do% {
  mu2 <- mean(data$u2[data$trialid==i])
  sdu2 <- sd(data$u2[data$trialid==i])
  df5 <- rbind(df5, data.frame(trialid=i, mean=mu2, sd=sdu2))
}
df5 <- df5[-1,]
# Add values back to data
foreach(i=1:num_studies) %do% {
  data$mu2[data$trialid==i] <- df5$mean[i]
  data$sdu2[data$trialid==i] <- df5$sd[i]
} 
foreach(i=1:pts) %do% {
  data$u2n[i] <- (data$u2[i]-data$mu2[i])/data$sdu2[i]
  data$du0[i] <- 1/data$sdv0[i]
  data$k[i] <- (3*(data$a[i]^2)) -((3*data$phi1[i])*(data$b[i]^2))- (3*(1-data$phi1[i])*(data$c[i]^2))
  data$l[i] <- data$sdu1[i]
  data$du1n[i] <- (1/data$l[i])*data$k[i] - (((data$v1u0[i]/data$u0u0[i])/data$l[i])*data$du0[i])
  data$m[i] <- (3*(data$e[i]^2)) - ((3*data$phi2[i])*(data$b[i]^2))- (3*(1-data$phi2[i])*(data$c[i]^2))
  data$n[i] <- data$sdu2[i]
  data$du2n[i] <- ((1/data$n[i])*data$m[i]) - (((data$v2u0[i]/data$u0u0[i])/data$n[i])*data$du0[i]) - (((data$v2u1n[i]/data$u1nu1n[[i]])/data$n[i])*data$du1n[i])
}
offset <- 0
foreach(i=1:num_studies) %do% {
  offset[i+1] <- offset[i] + num_patients[i]
}
# Prepare to fit model in WinBUGS
bugs_data <- list(d=data$event, trt2=data$trt2, trt3=data$trt3, trt4=data$trt4, trt5=data$trt5, trt6=data$trt6,trt7=data$trt7,
                  trt2lnt=data$trt2lnt, trt3lnt=data$trt3lnt, trt4lnt=data$trt4lnt, trt5lnt=data$trt5lnt,  trt6lnt=data$trt6lnt,trt7lnt=data$trt7lnt,
                  u0=data$u0, u1n=data$u1n, u2n=data$u2n, du0=data$du0,
                  du1n=data$du1n, du2n=data$du2n, Ntrials=num_studies, offset=offset,
                  k0=knots$k1[1], k1=knots$k2[1], k2=knots$k3[1], k3=knots$k4[1],
                  meanv0=df1$mean[1], sdv0=df1$sd[1], v1w0=df2$v1u0[1], w0w0=df2$u0u0[1],
                  meanw1=df3$mean[1], sdw1=df3$sd[1], v2w0=df4$v2u0[1], v2w1n=df4$v2u1n[1],
                  w1nw1n=df4$u1nu1n[1], meanw2=df5$mean[1], sdw2=df5$sd[1],
                  nt=nt)
# Create initial values for model
beta1 <- c(rep(0.1, 12))
beta2 <- c(rep(0.2, 12))
beta3 <- c(rep(0.3, 12))
gamma1 <- array(rep(c(0.1), 4*num_studies), dim=c(4,num_studies))
gamma2 <- array(rep(c(0.2), 4*num_studies), dim=c(4,num_studies))
gamma3 <- array(rep(c(0.3), 4*num_studies), dim=c(4,num_studies))
inits <- list(list(beta=beta1, gamma=gamma1), 
              list(beta=beta2, gamma=gamma2),
              list(beta=beta3, gamma=gamma3))
# Fit model in WinBUGS
bugs.object<-bugs(data=bugs_data, inits=inits, 
                  parameters.to.save=c("beta", "surv", "rk", "gamma"), 
                  model.file="model_pfs.txt", clearWD=F, 
                  summary.only=FALSE, n.iter=(num.sims+burn.in), n.sims=num.sims,
                  n.burnin=burn.in, n.chains=3, 
                  bugs.seed=1, bugs.directory=bugs.directory, 
                  debug=T, working.directory=path)
results <- bugs.object$summary
# Save results as a csv file
write.csv(results, "results_pfs.csv")


OS
OS
####----------------------------------------------------------------------
rm(list=ls())
data <- read.csv("RP_OS.csv")
path <- getwd()
setwd(path)
bugs.directory <- "C:/Users/28629/Documents/WinBUGS14"
# Set the location for WinBUGS
num.sims <- 5000
burn.in <- 5000
# Number of studies
num_studies <- length(unique(data$studycode))
# Sort data into trialid order
data <- data[order(data$studycode),]
# Number of patients in each study
# Note: if data is re-ordered will now to re-order this vector
num_patients <- c(596,645,514,747,649,659)####
# Number of treatments
nt <- length(unique(data$txcode))
# Create treatment indicator variables - caution needed to ensure consistency equations hold
data$trt2 <- 0
data$trt2[data$txcode==2] <- 1
data$trt3 <- 0
data$trt3[data$txcode==3] <- 1
data$trt4 <- 0
data$trt4[data$txcode==4] <- 1
data$trt5 <- 0
data$trt5[data$txcode==5] <- 1
data$trt6 <- 0
data$trt6[data$txcode==6] <- 1
data$trt7 <- 0
data$trt7[data$txcode==7] <- 1
# Add  trialid variable - easier than changing trialid to studyCode throughout this file!
data$trialid <- data$studycode
# Put eventtime onto the ln scale
data$eventtime <- data$time
data$lnt <- log(data$eventtime)
# Create trt*lnt variables
data$trt2lnt <- data$trt2*data$lnt
data$trt3lnt <- data$trt3*data$lnt
data$trt4lnt <- data$trt4*data$lnt
data$trt5lnt <- data$trt5*data$lnt
data$trt6lnt <- data$trt6*data$lnt
data$trt7lnt <- data$trt7*data$lnt
# Total number of patients
pts <- nrow(data)
# Set location of knots for each trial - 33rd and 67th percentiles of uncensored survival times
knots <- data.frame(trialid=NA, k1=NA, k2=NA, k3=NA, k4=NA)
foreach(i=1:num_studies) %do% {
  k <- quantile(data$lnt[data$event==1 & data$trialid==i], c(0, 0.33, 0.67, 1))
  knots <- rbind(knots, data.frame(trialid=i, k1=k[1], k2=k[2], k3=k[3], k4=k[4]))
}
# Remove the empty top row
knots <- knots[-1,]
# Add knot values to data
foreach(i=1:num_studies) %do% {
  data$k1[data$trialid==i] <- knots$k1[i]
  data$k2[data$trialid==i] <- knots$k2[i]
  data$k3[data$trialid==i] <- knots$k3[i]
  data$k4[data$trialid==i] <- knots$k4[i]
}  
# Orthogonalisation
data$t1 <- data$lnt-data$k2
data$t2 <- data$lnt-data$k1
data$t3 <- data$lnt-data$k4
data$t4 <- data$lnt-data$k3
foreach(i=1:pts) %do% {
  data$a[i] <- max(c(data$t1[i], 0))
  data$b[i] <- max(c(data$t2[i], 0))
  data$c[i] <- max(c(data$t3[i], 0))
  data$e[i] <- max(c(data$t4[i], 0))
  
  data$phi1[i] <- (data$k4[i]-data$k2[i])/(data$k4[i]-data$k1[i])
  data$phi2[i] <- (data$k4[i]-data$k3[i])/(data$k4[i]-data$k1[i])
  
  data$v0[i] <- data$lnt[i]
  data$v1[i] <- data$a[i]^3 - (data$phi1[i]*(data$b[i]^3)) - ((1-data$phi1[i])*(data$c[i]^3))
  data$v2[i] <- data$e[i]^3 - (data$phi2[i]*(data$b[i]^3)) - ((1-data$phi2[i])*(data$c[i]^3))
}
##--------------##
df1 <- data.frame(trialid=NA, mean=NA, sd=NA)
foreach(i=1:num_studies) %do% {
  mv0 <- mean(data$v0[data$trialid==i])
  sdv0 <- sd(data$v0[data$trialid==i])
  df1 <- rbind(df1, data.frame(trialid=i, mean=mv0, sd=sdv0))
}
df1 <- df1[-1,]
# Add values back to data
foreach(i=1:num_studies) %do% {
  data$mv0[data$trialid==i] <- df1$mean[i]
  data$sdv0[data$trialid==i] <- df1$sd[i]
} 
foreach(i=1:pts) %do% {
  data$u0[i] <- (data$v0[i]-data$mv0[i])/data$sdv0[i]
  data$f[i] <- data$v1[i]*data$u0[i]
  data$g[i] <- data$u0[i]*data$u0[i]
}
df2 <- data.frame(trialid=NA, v1u0=NA, u0u0=NA)
foreach(i=1:num_studies) %do% {
  v1u0 <- sum(data$f[data$trialid==i])
  u0u0 <- sum(data$g[data$trialid==i])
  df2 <- rbind(df2, data.frame(trialid=i, v1u0=v1u0, u0u0=u0u0))
}
df2 <- df2[-1,]
# Add values back to data
foreach(i=1:num_studies) %do% {
  data$v1u0[data$trialid==i] <- df2$v1u0[i]
  data$u0u0[data$trialid==i] <- df2$u0u0[i]
}
foreach(i=1:pts) %do% {
  data$u1[i] <- data$v1[i] - ((data$v1u0[i]/data$u0u0[i])*data$u0[i])
}
##--------------##
df3 <- data.frame(trialid=NA, mean=NA, sd=NA)
foreach(i=1:num_studies) %do% {
  mu1 <- mean(data$u1[data$trialid==i])
  sdu1 <- sd(data$u1[data$trialid==i])
  df3 <- rbind(df3, data.frame(trialid=i, mean=mu1, sd=sdu1))
}
df3 <- df3[-1,]
# Add values back to data
foreach(i=1:num_studies) %do% {
  data$mu1[data$trialid==i] <- df3$mean[i]
  data$sdu1[data$trialid==i] <- df3$sd[i]
} 
foreach(i=1:pts) %do% {
  data$u1n[i] <- (data$u1[i]-data$mu1[i])/data$sdu1[i]
  data$h[i] <- data$v2[i]*data$u0[i]
  data$aa[i] <- data$v2[i]*data$u1n[i]
  data$j[i] <- data$u1n[i]*data$u1n[i]
}
df4 <- data.frame(trialid=NA, v2u0=NA, v2u1n=NA, u1nu1n=NA)
foreach(i=1:num_studies) %do% {
  v2u0 <- sum(data$h[data$trialid==i])
  v2u1n <- sum(data$aa[data$trialid==i])
  u1nu1n <- sum(data$j[data$trialid==i])
  df4 <- rbind(df4, data.frame(trialid=i, v2u0=v2u0, v2u1n=v2u1n, u1nu1n=u1nu1n))
}
df4 <- df4[-1,]

# Add values back to data
foreach(i=1:num_studies) %do% {
  data$v2u0[data$trialid==i] <- df4$v2u0[i]
  data$v2u1n[data$trialid==i] <- df4$v2u1n[i]
  data$u1nu1n[data$trialid==i] <- df4$u1nu1n[i]
}
foreach(i=1:pts) %do% {
  data$u2[i] <- data$v2[i] - ((data$v2u0[i]/data$u0u0[i])*data$u0[i]) - ((data$v2u1n[i]/data$u1nu1n[i])*data$u1n[i])
}
##--------------##
df5 <- data.frame(trialid=NA, mean=NA, sd=NA)
foreach(i=1:num_studies) %do% {
  mu2 <- mean(data$u2[data$trialid==i])
  sdu2 <- sd(data$u2[data$trialid==i])
  df5 <- rbind(df5, data.frame(trialid=i, mean=mu2, sd=sdu2))
}
df5 <- df5[-1,]
# Add values back to data
foreach(i=1:num_studies) %do% {
  data$mu2[data$trialid==i] <- df5$mean[i]
  data$sdu2[data$trialid==i] <- df5$sd[i]
} 
foreach(i=1:pts) %do% {
  data$u2n[i] <- (data$u2[i]-data$mu2[i])/data$sdu2[i]
  data$du0[i] <- 1/data$sdv0[i]
  data$k[i] <- (3*(data$a[i]^2)) -((3*data$phi1[i])*(data$b[i]^2))- (3*(1-data$phi1[i])*(data$c[i]^2))
  data$l[i] <- data$sdu1[i]
  data$du1n[i] <- (1/data$l[i])*data$k[i] - (((data$v1u0[i]/data$u0u0[i])/data$l[i])*data$du0[i])
  data$m[i] <- (3*(data$e[i]^2)) - ((3*data$phi2[i])*(data$b[i]^2))- (3*(1-data$phi2[i])*(data$c[i]^2))
  data$n[i] <- data$sdu2[i]
  data$du2n[i] <- ((1/data$n[i])*data$m[i]) - (((data$v2u0[i]/data$u0u0[i])/data$n[i])*data$du0[i]) - (((data$v2u1n[i]/data$u1nu1n[[i]])/data$n[i])*data$du1n[i])
}
offset <- 0
foreach(i=1:num_studies) %do% {
  offset[i+1] <- offset[i] + num_patients[i]
}
# Prepare to fit model in WinBUGS
bugs_data <- list(d=data$event, trt2=data$trt2, trt3=data$trt3, trt4=data$trt4, trt5=data$trt5, trt6=data$trt6,trt7=data$trt7,
                  trt2lnt=data$trt2lnt, trt3lnt=data$trt3lnt, trt4lnt=data$trt4lnt, trt5lnt=data$trt5lnt,  trt6lnt=data$trt6lnt,trt7lnt=data$trt7lnt,
                  u0=data$u0, u1n=data$u1n, u2n=data$u2n, du0=data$du0,
                  du1n=data$du1n, du2n=data$du2n, Ntrials=num_studies, offset=offset,
                  k0=knots$k1[1], k1=knots$k2[1], k2=knots$k3[1], k3=knots$k4[1],
                  meanv0=df1$mean[1], sdv0=df1$sd[1], v1w0=df2$v1u0[1], w0w0=df2$u0u0[1],
                  meanw1=df3$mean[1], sdw1=df3$sd[1], v2w0=df4$v2u0[1], v2w1n=df4$v2u1n[1],
                  w1nw1n=df4$u1nu1n[1], meanw2=df5$mean[1], sdw2=df5$sd[1],
                  nt=nt)
# Create initial values for model
beta1 <- c(rep(0.1, 12))
beta2 <- c(rep(0.2, 12))
beta3 <- c(rep(0.3, 12))
gamma1 <- array(rep(c(0.1), 4*num_studies), dim=c(4,num_studies))
gamma2 <- array(rep(c(0.2), 4*num_studies), dim=c(4,num_studies))
gamma3 <- array(rep(c(0.3), 4*num_studies), dim=c(4,num_studies))
inits <- list(list(beta=beta1, gamma=gamma1), 
              list(beta=beta2, gamma=gamma2),
              list(beta=beta3, gamma=gamma3))
# Fit model in WinBUGS
bugs.object<-bugs(data=bugs_data, inits=inits, 
                  parameters.to.save=c("beta", "surv", "rk", "gamma"), 
                  model.file="model_os.txt", clearWD=F, 
                  summary.only=FALSE, n.iter=(num.sims+burn.in), n.sims=num.sims,
                  n.burnin=burn.in, n.chains=3, 
                  bugs.seed=1, bugs.directory=bugs.directory, 
                  debug=T, working.directory=path)
results <- bugs.object$summary
# Save results as a csv file
write.csv(results, "results_os.csv")

#################################################
data <- read.csv("results_pfs.csv")

data1 <- data[853:2616,1:2]%>%
  mutate(months=rep(1:36,each =49),rank=rep(rep(1:7,times=7),times=36),treatment=rep(rep(c("che","cam+che","niv+che","tori+che","pemb+che","tisle+che","sinti+che"),each= 7),times=36))%>%select(-X)%>%filter(rank==7)

data1 <- data1%>%filter(months<=36)%>%ggplot(aes(x=months, y=mean,group=treatment, color=treatment)) +
  geom_point()+geom_line()+
  scale_x_continuous(name="Time(months)", limits=c(0, 36), breaks=seq(0,36,3))+
  scale_y_continuous(name="Progression-free Survival",limits = c(0,1.0))+
  geom_vline(xintercept = 12,color="gray",linetype='dashed')+
  ggsci::scale_color_jama()+theme_bw()+
  labs(title= "Probability to be best for PFS")

data2 <- data[853:2616,1:2]%>%
  mutate(months=rep(1:36,each =49),rank=rep(rep(1:7,times=7),times=36),treatment=rep(rep(c("che","cam+che","niv+che","tori+che","pemb+che","tisle+che","sinti+che"),each = 7),times=36))%>%select(-X)%>%filter(rank==6)

data2 <- data2%>%filter(months<=36)%>%ggplot(aes(x=months, y=mean,group=treatment, color=treatment)) +
  geom_point()+geom_line()+
  scale_x_continuous(name="Time(months)", limits=c(0, 36), breaks=seq(0,36,3))+
  scale_y_continuous(name="Progression-free Survival",limits = c(0,1.0))+
  geom_vline(xintercept = 12,color="gray",linetype='dashed')+
  ggsci::scale_color_jama()+theme_bw()+
  labs(title= "Probability to be second-best for PFS")

#################################################
data <- read.csv("results_os.csv")


data3 <- data[853:2616,1:2]%>%
  mutate(months=rep(1:36,each =49),rank=rep(rep(1:7,times=7),times=36),treatment=rep(rep(c("che","cam+che","niv+che","tori+che","pemb+che","tisle+che","sinti+che"),each=7),times=36))%>%select(-X)%>%filter(rank==7)

data3 <- data3%>%filter(months<=36)%>%ggplot(aes(x=months, y=mean,group=treatment, color=treatment)) +
  geom_point()+geom_line()+
  scale_x_continuous(name="Time(months)", limits=c(0, 36), breaks=seq(0,36,3))+
  scale_y_continuous(name="Overall Survival",limits = c(0,1.0))+
  geom_vline(xintercept = 12,color="gray",linetype='dashed')+
  ggsci::scale_color_jama()+theme_bw()+
  labs(title= "Probability to be best for OS")

data4 <- data[853:2616,1:2]%>%
  mutate(months=rep(1:36,each =49),rank=rep(rep(1:7,times=7),times=36),treatment=rep(rep(c("che","cam+che","niv+che","tori+che","pemb+che","tisle+che","sinti+che"),each=7),times=36))%>%select(-X)%>%filter(rank==6)

data4 <- data4%>%filter(months<=36)%>%ggplot(aes(x=months, y=mean,group=treatment, color=treatment)) +
  geom_point()+geom_line()+
  scale_x_continuous(name="Time(months)", limits=c(0, 36), breaks=seq(0,36,3))+
  scale_y_continuous(name="Overall Survival",limits = c(0,1.0))+
  geom_vline(xintercept = 12,color="gray",linetype='dashed')+
  ggsci::scale_color_jama()+theme_bw()+
  labs(title= "Probability to be second-best for OS")
##############################################################
library(patchwork)

a <- data1+data2+data3+data4+plot_annotation(tag_levels = "A")
a
export::graph2tif(file="figure 2",width=12,height=12)