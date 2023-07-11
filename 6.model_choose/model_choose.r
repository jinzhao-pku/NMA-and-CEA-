rm(list=ls())
library(survival)
library(tidyverse)
library(flexsurv)
library(survminer)
data = list ()
fig1 = list ()
fig2 = list ()
fig3 = list ()
fig4 = list ()
AIC  = list ()
BIC  = list ()
OS   = list ()
PFS  = list ()
dt <- read.csv("RP_OS.csv")
for (i in 1:7) {
  k = i
  data[[i]] <- dt%>%filter(txcode==k)
  data1 = data[[i]]
  out_km <- survfit(Surv(time,event)~1,data=data1)
  fit_wb <- flexsurvreg(Surv(time,event)~1,data= data1,dist="weibull")
  t <- seq(0, 120, by = 7/30)
  s <- predict(fit_wb, type = "survival", times = t, newdata = data1[1, ])
  out_wb <- data.frame(time=t,surv=s[[1]][[1]]$.pred_survival,class="Weibull",AIC=-2*fit_wb$loglik+2*fit_wb$npars,BIC=-2*fit_wb$loglik+fit_wb$npars*log(fit_wb$N))
  
  fit_exp <- flexsurvreg(Surv(time,event)~1,data= data1,dist="Exponential")
  t <- seq(0, 36, by = 7/30)
  s <- predict(fit_exp, type = "survival", times = t, newdata = data1[1, ])
  out_exp <- data.frame(time=t,surv=s[[1]][[1]]$.pred_survival,class="Exponential",AIC=-2*fit_exp$loglik+2*fit_exp$npars,BIC=-2*fit_exp$loglik+fit_exp$npars*log(fit_exp$N))
  
  fit_lnorm <- flexsurvreg(Surv(time,event)~1,data= data1,dist="lnorm")
  t <- seq(0, 120, by = 7/30)
  s <- predict(fit_lnorm, type = "survival", times = t, newdata = data1[1, ])
  out_lnorm <- data.frame(time=t,surv=s[[1]][[1]]$.pred_survival,class="LogNormal",AIC=-2*fit_lnorm$loglik+2*fit_lnorm$npars,BIC=-2*fit_lnorm$loglik+fit_lnorm$npars*log(fit_lnorm$N))
  
  fit_gengamma <- flexsurvreg(Surv(time,event)~1,data= data1,dist="gengamma")
  t <- seq(0, 120, by = 7/30)
  s <- predict(fit_gengamma, type = "survival", times = t, newdata = data1[1, ])
  out_gengamma <- data.frame(time=t,surv=s[[1]][[1]]$.pred_survival,class="Generalized gamma",AIC=-2*fit_gengamma$loglik+2*fit_gengamma$npars,BIC=-2*fit_gengamma$loglik+fit_gengamma$npars*log(fit_gengamma$N))
  
  fit_llogis <- flexsurvreg(Surv(time,event)~1,data= data1,dist="llogis")
  t <- seq(0, 120, by = 7/30)
  s <- predict(fit_llogis, type = "survival", times = t, newdata = data1[1, ])
  out_llogis <- data.frame(time=t,surv=s[[1]][[1]]$.pred_survival,class="Loglogistic",AIC=-2*fit_llogis$loglik+2*fit_llogis$npars,BIC=-2*fit_llogis$loglik+fit_llogis$npars*log(fit_llogis$N))
  
  fit_gompertz <- flexsurvreg(Surv(time,event)~1,data= data1,dist="gompertz")
  t <- seq(0, 120, by = 7/30)
  s <- predict(fit_gompertz, type = "survival", times = t, newdata = data1[1, ])
  out_gompertz <- data.frame(time=t,surv=s[[1]][[1]]$.pred_survival,class="Gompertz",AIC=-2*fit_gompertz$loglik+2*fit_gompertz$npars,BIC=-2*fit_gompertz$loglik+fit_gompertz$npars*log(fit_gompertz$N))
  
  fit_RP2odds<- flexsurvspline(Surv(time,event)~1, data=data1, k=2, scale="odds")
  t <- seq(0, 120, by = 7/30)
  s <- predict(fit_RP2odds, type = "survival", times = t, newdata = data1[1, ])
  out_RP2odds <- data.frame(time=t,surv=s[[1]][[1]]$.pred_survival,class="Royston and Parmar-odds 2",AIC=-2*fit_RP2odds$loglik+2*fit_RP2odds$npars,BIC=-2*fit_RP2odds$loglik+fit_RP2odds$npars*log(fit_RP2odds$N))
  
  fit_RP3odds<- flexsurvspline(Surv(time,event)~1, data=data1, k=3, scale="odds")
  t <- seq(0, 120, by = 7/30)
  s <- predict(fit_RP3odds, type = "survival", times = t, newdata = data1[1, ])
  out_RP3odds <- data.frame(time=t,surv=s[[1]][[1]]$.pred_survival,class="Royston and Parmar-odds 3",AIC=-2*fit_RP3odds$loglik+2*fit_RP3odds$npars,BIC=-2*fit_RP3odds$loglik+fit_RP3odds$npars*log(fit_RP3odds$N))
  
  fit_RP2hazard<- flexsurvspline(Surv(time,event)~1, data=data1, k=2, scale="hazard")
  t <- seq(0, 120, by = 7/30)
  s <- predict(fit_RP2hazard, type = "survival", times = t, newdata = data1[1, ])
  out_RP2hazard <- data.frame(time=t,surv=s[[1]][[1]]$.pred_survival,class="Royston and Parmar-hazard 2",AIC=-2*fit_RP2hazard$loglik+2*fit_RP2hazard$npars,BIC=-2*fit_RP2hazard$loglik+fit_RP2hazard$npars*log(fit_RP2hazard$N))
  
  fit_RP3hazard<- flexsurvspline(Surv(time,event)~1, data=data1, k=3, scale="hazard")
  t <- seq(0, 120, by = 7/30)
  s <- predict(fit_RP3hazard, type = "survival", times = t, newdata = data1[1, ])
  out_RP3hazard <- data.frame(time=t,surv=s[[1]][[1]]$.pred_survival,class="Royston and Parmar-hazard 3",AIC=-2*fit_RP3hazard$loglik+2*fit_RP3hazard$npars,BIC=-2*fit_RP3hazard$loglik+fit_RP3hazard$npars*log(fit_RP3hazard$N))
  
  fit_RP2normal<- flexsurvspline(Surv(time,event)~1, data=data1, k=2, scale="normal")
  t <- seq(0, 120, by = 7/30)
  s <- predict(fit_RP2normal, type = "survival", times = t, newdata = data1[1, ])
  out_RP2normal <- data.frame(time=t,surv=s[[1]][[1]]$.pred_survival,class="Royston and Parmar-normal 2",AIC=-2*fit_RP2normal$loglik+2*fit_RP2normal$npars,BIC=-2*fit_RP2normal$loglik+fit_RP2normal$npars*log(fit_RP2normal$N))
  
  fit_RP3normal<- flexsurvspline(Surv(time,event)~1, data=data1, k=3, scale="normal")
  t <- seq(0, 120, by = 7/30)
  s <- predict(fit_RP3normal, type = "survival", times = t, newdata = data1[1, ])
  out_RP3normal <- data.frame(time=t,surv=s[[1]][[1]]$.pred_survival,class="Royston and Parmar-normal 3",AIC=-2*fit_RP3normal$loglik+2*fit_RP3normal$npars,BIC=-2*fit_RP3normal$loglik+fit_RP3normal$npars*log(fit_RP3normal$N))
  
  out1 <- rbind(out_exp,out_lnorm,out_llogis,out_gengamma,out_gompertz,out_wb)
  out2 <- rbind(out_RP2odds,out_RP3odds,out_RP2hazard,out_RP3hazard,out_RP2normal,out_RP3normal)
  out <-  rbind(out1,out2)
  
  
  fig1[[i]]= ggsurvplot(out_km,xlim=c(0,120),ci=T,conf.int=F,legend.title="",legend.labs="KaplanMeier curves",break.time.by = 6)$plot+
    geom_line(data = out1, aes(x = time,y=surv,color=class),linewidth=0.02) + theme(legend.position = "right") +
    scale_color_manual(data1$trt[1], values= c("#ffe119", "#4363d8", "#800000", "#dcbeff", "#a9a9a9", "#000075", "#f58231"))
  
  fig2[[i]]= ggsurvplot(out_km,xlim=c(0,120),ci=T,conf.int=F,legend.title="",legend.labs="KaplanMeier curves",break.time.by = 6)$plot+
    geom_line(data = out2, aes(x = time,y=surv,color=class),linewidth=0.02) + theme(legend.position = "right") +
    scale_color_manual(data1$trt[1], values= c("#ffe119", "#4363d8", "#800000", "#dcbeff", "#a9a9a9", "#000075", "#f58231"))
  if(i==1)
    {OS[[i]] <- out_RP2hazard[1:2]}
  if(i!=1){OS[[i]] <- out_llogis[1:2]}
  AIC[[i]] <-out%>%distinct(class,AIC) 
  BIC[[i]] <-out%>%distinct(class,BIC) 
}
mix1 = ggarrange(fig1[[1]],fig2[[1]],
                 fig1[[2]],fig2[[2]],
                 fig1[[3]],fig2[[3]],
                 fig1[[4]],fig2[[4]],
                 fig1[[5]],fig2[[5]],
                 fig1[[6]],fig2[[6]],
                 fig1[[7]],fig2[[7]],
                 ncol = 2, nrow = 7)
mix1
export::graph2pdf(file="figure 2",width=18,height=24)
write.csv(AIC,"AIC-OS.csv")
write.csv(BIC,"BIC-OS.csv")


dt <- read.csv("RP_PFS.csv")
for (i in 1:7) {
  k = i
  data[[i]] <- dt%>%filter(txcode==k)
  data1 = data[[i]]
  out_km <- survfit(Surv(time,event)~1,data=data1)
  fit_wb <- flexsurvreg(Surv(time,event)~1,data= data1,dist="weibull")
  t <- seq(0, 120, by = 7/30)
  s <- predict(fit_wb, type = "survival", times = t, newdata = data1[1, ])
  out_wb <- data.frame(time=t,surv=s[[1]][[1]]$.pred_survival,class="Weibull",AIC=-2*fit_wb$loglik+2*fit_wb$npars,BIC=-2*fit_wb$loglik+fit_wb$npars*log(fit_wb$N))
  
  fit_exp <- flexsurvreg(Surv(time,event)~1,data= data1,dist="Exponential")
  t <- seq(0, 120, by = 7/30)
  s <- predict(fit_exp, type = "survival", times = t, newdata = data1[1, ])
  out_exp <- data.frame(time=t,surv=s[[1]][[1]]$.pred_survival,class="Exponential",AIC=-2*fit_exp$loglik+2*fit_exp$npars,BIC=-2*fit_exp$loglik+fit_exp$npars*log(fit_exp$N))
  
  fit_lnorm <- flexsurvreg(Surv(time,event)~1,data= data1,dist="lnorm")
  t <- seq(0, 120, by = 7/30)
  s <- predict(fit_lnorm, type = "survival", times = t, newdata = data1[1, ])
  out_lnorm <- data.frame(time=t,surv=s[[1]][[1]]$.pred_survival,class="LogNormal",AIC=-2*fit_lnorm$loglik+2*fit_lnorm$npars,BIC=-2*fit_lnorm$loglik+fit_lnorm$npars*log(fit_lnorm$N))
  
  fit_gengamma <- flexsurvreg(Surv(time,event)~1,data= data1,dist="gengamma")
  t <- seq(0, 120, by = 7/30)
  s <- predict(fit_gengamma, type = "survival", times = t, newdata = data1[1, ])
  out_gengamma <- data.frame(time=t,surv=s[[1]][[1]]$.pred_survival,class="Generalized gamma",AIC=-2*fit_gengamma$loglik+2*fit_gengamma$npars,BIC=-2*fit_gengamma$loglik+fit_gengamma$npars*log(fit_gengamma$N))
  
  fit_llogis <- flexsurvreg(Surv(time,event)~1,data= data1,dist="llogis")
  t <- seq(0, 120, by = 7/30)
  s <- predict(fit_llogis, type = "survival", times = t, newdata = data1[1, ])
  out_llogis <- data.frame(time=t,surv=s[[1]][[1]]$.pred_survival,class="Loglogistic",AIC=-2*fit_llogis$loglik+2*fit_llogis$npars,BIC=-2*fit_llogis$loglik+fit_llogis$npars*log(fit_llogis$N))
  
  fit_gompertz <- flexsurvreg(Surv(time,event)~1,data= data1,dist="gompertz")
  t <- seq(0, 120, by = 7/30)
  s <- predict(fit_gompertz, type = "survival", times = t, newdata = data1[1, ])
  out_gompertz <- data.frame(time=t,surv=s[[1]][[1]]$.pred_survival,class="Gompertz",AIC=-2*fit_gompertz$loglik+2*fit_gompertz$npars,BIC=-2*fit_gompertz$loglik+fit_gompertz$npars*log(fit_gompertz$N))
  
  fit_RP2odds<- flexsurvspline(Surv(time,event)~1, data=data1, k=2, scale="odds")
  t <- seq(0, 120, by = 0.035)
  s <- predict(fit_RP2odds, type = "survival", times = t, newdata = data1[1, ])
  out_RP2odds <- data.frame(time=t,surv=s[[1]][[1]]$.pred_survival,class="Royston and Parmar-odds 2",AIC=-2*fit_RP2odds$loglik+2*fit_RP2odds$npars,BIC=-2*fit_RP2odds$loglik+fit_RP2odds$npars*log(fit_RP2odds$N))
  
  fit_RP3odds<- flexsurvspline(Surv(time,event)~1, data=data1, k=3, scale="odds")
  t <- seq(0, 120, by = 7/30)
  s <- predict(fit_RP3odds, type = "survival", times = t, newdata = data1[1, ])
  out_RP3odds <- data.frame(time=t,surv=s[[1]][[1]]$.pred_survival,class="Royston and Parmar-odds 3",AIC=-2*fit_RP3odds$loglik+2*fit_RP3odds$npars,BIC=-2*fit_RP3odds$loglik+fit_RP3odds$npars*log(fit_RP3odds$N))
  
  fit_RP2hazard<- flexsurvspline(Surv(time,event)~1, data=data1, k=2, scale="hazard")
  t <- seq(0, 120, by = 7/30)
  s <- predict(fit_RP2hazard, type = "survival", times = t, newdata = data1[1, ])
  out_RP2hazard <- data.frame(time=t,surv=s[[1]][[1]]$.pred_survival,class="Royston and Parmar-hazard 2",AIC=-2*fit_RP2hazard$loglik+2*fit_RP2hazard$npars,BIC=-2*fit_RP2hazard$loglik+fit_RP2hazard$npars*log(fit_RP2hazard$N))
  
  fit_RP3hazard<- flexsurvspline(Surv(time,event)~1, data=data1, k=3, scale="hazard")
  t <- seq(0, 120, by = 7/30)
  s <- predict(fit_RP3hazard, type = "survival", times = t, newdata = data1[1, ])
  out_RP3hazard <- data.frame(time=t,surv=s[[1]][[1]]$.pred_survival,class="Royston and Parmar-hazard 3",AIC=-2*fit_RP3hazard$loglik+2*fit_RP3hazard$npars,BIC=-2*fit_RP3hazard$loglik+fit_RP3hazard$npars*log(fit_RP3hazard$N))
  
  fit_RP2normal<- flexsurvspline(Surv(time,event)~1, data=data1, k=2, scale="normal")
  t <- seq(0, 120, by = 7/30)
  s <- predict(fit_RP2normal, type = "survival", times = t, newdata = data1[1, ])
  out_RP2normal <- data.frame(time=t,surv=s[[1]][[1]]$.pred_survival,class="Royston and Parmar-normal 2",AIC=-2*fit_RP2normal$loglik+2*fit_RP2normal$npars,BIC=-2*fit_RP2normal$loglik+fit_RP2normal$npars*log(fit_RP2normal$N))
  
  fit_RP3normal<- flexsurvspline(Surv(time,event)~1, data=data1, k=3, scale="normal")
  t <- seq(0, 120, by = 7/30)
  s <- predict(fit_RP3normal, type = "survival", times = t, newdata = data1[1, ])
  out_RP3normal <- data.frame(time=t,surv=s[[1]][[1]]$.pred_survival,class="Royston and Parmar-normal 3",AIC=-2*fit_RP3normal$loglik+2*fit_RP3normal$npars,BIC=-2*fit_RP3normal$loglik+fit_RP3normal$npars*log(fit_RP3normal$N))
  
  out1 <- rbind(out_exp,out_lnorm,out_llogis,out_gengamma,out_gompertz,out_wb)
  out2 <- rbind(out_RP2odds,out_RP3odds,out_RP2hazard,out_RP3hazard,out_RP2normal,out_RP3normal)
  out <-  rbind(out1,out2)
  
  
  fig1[[i]]= ggsurvplot(out_km,xlim=c(0,120),ci=T,conf.int=F,legend.title="",legend.labs="KaplanMeier curves",break.time.by = 6)$plot+
    geom_line(data = out1, aes(x = time,y=surv,color=class),linewidth=0.02) + theme(legend.position = "right") +
    scale_color_manual(data1$trt[1], values= c("#ffe119", "#4363d8", "#800000", "#dcbeff", "#a9a9a9", "#000075", "#f58231"))
  
  fig2[[i]]= ggsurvplot(out_km,xlim=c(0,120),ci=T,conf.int=F,legend.title="",legend.labs="KaplanMeier curves",break.time.by = 6)$plot+
    geom_line(data = out2, aes(x = time,y=surv,color=class),linewidth=0.02) + theme(legend.position = "right") +
    scale_color_manual(data1$trt[1], values= c("#ffe119", "#4363d8", "#800000", "#dcbeff", "#a9a9a9", "#000075", "#f58231"))
  
  AIC[[i]] <-out%>%distinct(class,AIC) 
  BIC[[i]] <-out%>%distinct(class,BIC) 
  if(i==1|i==4|i==5|i==6)
    {PFS[[i]] <- out_RP3normal[1:2]}
  if(i==2)
    {PFS[[i]] <- out_RP2hazard[1:2]}
  if(i==3)
    {PFS[[i]] <- out_llogis[1:2]}
  if(i==7)
    {PFS[[i]] <- out_RP2normal[1:2]}
}
mix1 = ggarrange(fig1[[1]],fig2[[1]],
                 fig1[[2]],fig2[[2]],
                 fig1[[3]],fig2[[3]],
                 fig1[[4]],fig2[[4]],
                 fig1[[5]],fig2[[5]],
                 fig1[[6]],fig2[[6]],
                 fig1[[7]],fig2[[7]],
                 ncol = 2, nrow = 7)
mix1
export::graph2pdf(file="figure 1",width=18,height=24)
write.csv(AIC,"AIC-PFS.csv")
write.csv(BIC,"BIC-PFS.csv")
write.csv(OS,"OS.csv")
write.csv(PFS,"PFS.csv")






