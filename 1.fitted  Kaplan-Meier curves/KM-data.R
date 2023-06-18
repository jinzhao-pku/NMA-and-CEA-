rm(list=ls())
library(survHE)
library(survival)
library(survminer)
library(readr)
setwd("C:/Users/28629/Documents/OS")


###1
##OS for ESCORT-1st
surv_inp <- "C_os_surv.txt"
nrisk_inp <- "C_os_table.txt"
digitise(surv_inp,nrisk_inp,km_output = "KM_camel_os0.txt",ipd_output = "IPD_camel_os0.txt")
surv_inp <- "I_os_surv.txt"
nrisk_inp <- "I_os_table.txt"
digitise(surv_inp,nrisk_inp,km_output = "KM_camel_os1.txt",ipd_output = "IPD_camel_os1.txt")
ipd_files <- list("IPD_camel_os0.txt","IPD_camel_os1.txt")
data <- make.ipd(ipd_files, ctr=1 ,var.labs = c("time","event","arm"))
#fig
max_xlim<-as.integer(max(data$time))+1

###Cox regression and KM curves
fit <- survfit(Surv(time,event==1)~arm,data=data)
res.cox <- coxph(Surv(time, event==1) ~ arm, data =data)
summary(res.cox)
ggsurvplot(fit,data=data,pval = TRUE,risk.table = TRUE,
           palette = c("#D95F02","steelblue"),
           legend.tittle="Intervention",
           legend.labs=c("Control","Intervention"),
           risk.table.height=0.2,
           xlab=c("Time(Months))"),
           break.y.by=0.2,
           break.x.by=2,
           xlim=c(0,max_xlim),
           ylab=c("OS"),
           risk.table.col="strata",
           tables.theme=theme_cleantable(),
           ggtheme=theme_classic(),
           surv.median.line = "hv")

###PFS for ESCORT-1st
surv_inp <- "C_pfs_surv.txt"
nrisk_inp <- "C_pfs_table.txt"
digitise(surv_inp,nrisk_inp,km_output = "KM_camel_pfs0.txt",ipd_output = "IPD_camel_pfs0.txt")
surv_inp <- "I_pfs_surv.txt"
nrisk_inp <- "I_pfs_table.txt"
digitise(surv_inp,nrisk_inp,km_output = "KM_camel_pfs1.txt",ipd_output = "IPD_camel_pfs1.txt")
ipd_files <- list("IPD_camel_pfs0.txt","IPD_camel_pfs1.txt")
data <- make.ipd(ipd_files, ctr=1 ,var.labs = c("time","event","arm"))
#fig
max_xlim<-as.integer(max(data$time))+1

#Cox regression and KM curves
fit <- survfit(Surv(time,event==1)~arm,data=data)
res.cox <- coxph(Surv(time, event==1) ~ arm, data =data)
summary(res.cox)
ggsurvplot(fit,data=data,pval = TRUE,risk.table = TRUE,
           palette = c("#D95F02","steelblue"),
           legend.tittle="Intervention",
           legend.labs=c("Control","Intervention"),
           risk.table.height=0.2,
           xlab=c("Time(Months))"),
           break.y.by=0.2,
           break.x.by=2,
           xlim=c(0,max_xlim),
           ylab=c("PFS"),
           risk.table.col="strata",
           tables.theme=theme_cleantable(),
           ggtheme=theme_classic(),
           surv.median.line = "hv")
###2
##OS for CheckMate 648
surv_inp <- "C_os_surv.txt"
nrisk_inp <- "C_os_table.txt"
digitise(surv_inp,nrisk_inp,km_output = "KM_nivol_os0.txt",ipd_output = "IPD_nivol_os0.txt")
surv_inp <- "I_os_surv.txt"
nrisk_inp <- "I_os_table.txt"
digitise(surv_inp,nrisk_inp,km_output = "KM_nivol_os1.txt",ipd_output = "IPD_nivol_os1.txt")
ipd_files <- list("IPD_nivol_os0.txt","IPD_nivol_os1.txt")
data <- make.ipd(ipd_files, ctr=1 ,var.labs = c("time","event","arm"))
#fig
max_xlim<-as.integer(max(data$time))+1

###Cox regression and KM curves
fit <- survfit(Surv(time,event==1)~arm,data=data)
res.cox <- coxph(Surv(time, event==1) ~ arm, data =data)
summary(res.cox)
ggsurvplot(fit,data=data,pval = TRUE,risk.table = TRUE,
           palette = c("#D95F02","steelblue"),
           legend.tittle="Intervention",
           legend.labs=c("Control","Intervention"),
           risk.table.height=0.2,
           xlab=c("Time(Months))"),
           break.y.by=0.2,
           break.x.by=2,
           xlim=c(0,max_xlim),
           ylab=c("OS"),
           risk.table.col="strata",
           tables.theme=theme_cleantable(),
           ggtheme=theme_classic(),
           surv.median.line = "hv")

###PFS for CheckMate 648
surv_inp <- "C_pfs_surv.txt"
nrisk_inp <- "C_pfs_table.txt"
digitise(surv_inp,nrisk_inp,km_output = "KM_nivol_pfs0.txt",ipd_output = "IPD_nivol_pfs0.txt")
surv_inp <- "I_pfs_surv.txt"
nrisk_inp <- "I_pfs_table.txt"
digitise(surv_inp,nrisk_inp,km_output = "KM_nivol_pfs1.txt",ipd_output = "IPD_nivol_pfs1.txt")
ipd_files <- list("IPD_nivol_pfs0.txt","IPD_nivol_pfs1.txt")
data <- make.ipd(ipd_files, ctr=1 ,var.labs = c("time","event","arm"))
#fig
max_xlim<-as.integer(max(data$time))+1

#Cox regression and KM curves
fit <- survfit(Surv(time,event==1)~arm,data=data)
res.cox <- coxph(Surv(time, event==1) ~ arm, data =data)
summary(res.cox)
ggsurvplot(fit,data=data,pval = TRUE,risk.table = TRUE,
           palette = c("#D95F02","steelblue"),
           legend.tittle="Intervention",
           legend.labs=c("Control","Intervention"),
           risk.table.height=0.2,
           xlab=c("Time(Months))"),
           break.y.by=0.2,
           break.x.by=2,
           xlim=c(0,max_xlim),
           ylab=c("PFS"),
           risk.table.col="strata",
           tables.theme=theme_cleantable(),
           ggtheme=theme_classic(),
           surv.median.line = "hv")

###3
##OS for ORIENT-15
surv_inp <- "C_os_surv.txt"
nrisk_inp <- "C_os_table.txt"
digitise(surv_inp,nrisk_inp,km_output = "KM_sinti_os0.txt",ipd_output = "IPD_sinti_os0.txt")
surv_inp <- "I_os_surv.txt"
nrisk_inp <- "I_os_table.txt"
digitise(surv_inp,nrisk_inp,km_output = "KM_sinti_os1.txt",ipd_output = "IPD_sinti_os1.txt")
ipd_files <- list("IPD_sinti_os0.txt","IPD_sinti_os1.txt")
data <- make.ipd(ipd_files, ctr=1 ,var.labs = c("time","event","arm"))
#fig
max_xlim<-as.integer(max(data$time))+1

###Cox regression and KM curves
fit <- survfit(Surv(time,event==1)~arm,data=data)
res.cox <- coxph(Surv(time, event==1) ~ arm, data =data)
summary(res.cox)
ggsurvplot(fit,data=data,pval = TRUE,risk.table = TRUE,
           palette = c("#D95F02","steelblue"),
           legend.tittle="Intervention",
           legend.labs=c("Control","Intervention"),
           risk.table.height=0.2,
           xlab=c("Time(Months))"),
           break.y.by=0.2,
           break.x.by=2,
           xlim=c(0,max_xlim),
           ylab=c("OS"),
           risk.table.col="strata",
           tables.theme=theme_cleantable(),
           ggtheme=theme_classic(),
           surv.median.line = "hv")

###PFS for ORIENT-15
surv_inp <- "C_pfs_surv.txt"
nrisk_inp <- "C_pfs_table.txt"
digitise(surv_inp,nrisk_inp,km_output = "KM_sinti_pfs0.txt",ipd_output = "IPD_sinti_pfs0.txt")
surv_inp <- "I_pfs_surv.txt"
nrisk_inp <- "I_pfs_table.txt"
digitise(surv_inp,nrisk_inp,km_output = "KM_sinti_pfs1.txt",ipd_output = "IPD_sinti_pfs1.txt")
ipd_files <- list("IPD_sinti_pfs0.txt","IPD_sinti_pfs1.txt")
data <- make.ipd(ipd_files, ctr=1 ,var.labs = c("time","event","arm"))
#fig
max_xlim<-as.integer(max(data$time))+1

#Cox regression and KM curves
fit <- survfit(Surv(time,event==1)~arm,data=data)
res.cox <- coxph(Surv(time, event==1) ~ arm, data =data)
summary(res.cox)
ggsurvplot(fit,data=data,pval = TRUE,risk.table = TRUE,
           palette = c("#D95F02","steelblue"),
           legend.tittle="Intervention",
           legend.labs=c("Control","Intervention"),
           risk.table.height=0.2,
           xlab=c("Time(Months))"),
           break.y.by=0.2,
           break.x.by=2,
           xlim=c(0,max_xlim),
           ylab=c("PFS"),
           risk.table.col="strata",
           tables.theme=theme_cleantable(),
           ggtheme=theme_classic(),
           surv.median.line = "hv")

###4
##OS for KEYNOTE-590
surv_inp <- "C_os_surv.txt"
nrisk_inp <- "C_os_table.txt"
digitise(surv_inp,nrisk_inp,km_output = "KM_pembr_os0.txt",ipd_output = "IPD_pembr_os0.txt")
surv_inp <- "I_os_surv.txt"
nrisk_inp <- "I_OS_table.txt"
digitise(surv_inp,nrisk_inp,km_output = "KM_pembr_os1.txt",ipd_output = "IPD_pembr_os1.txt")
ipd_files <- list("IPD_pembr_os0.txt","IPD_pembr_os1.txt")
data <- make.ipd(ipd_files, ctr=1 ,var.labs = c("time","event","arm"))
#fig
max_xlim<-as.integer(max(data$time))+1

###Cox regression and KM curves
fit <- survfit(Surv(time,event==1)~arm,data=data)
res.cox <- coxph(Surv(time, event==1) ~ arm, data =data)
summary(res.cox)
ggsurvplot(fit,data=data,pval = TRUE,risk.table = TRUE,
           palette = c("#D95F02","steelblue"),
           legend.tittle="Intervention",
           legend.labs=c("Control","Intervention"),
           risk.table.height=0.2,
           xlab=c("Time(Months))"),
           break.y.by=0.2,
           break.x.by=2,
           xlim=c(0,max_xlim),
           ylab=c("OS"),
           risk.table.col="strata",
           tables.theme=theme_cleantable(),
           ggtheme=theme_classic(),
           surv.median.line = "hv")

##PFS for KEYNOTE-590
surv_inp <- "C_pfs_surv.txt"
nrisk_inp <- "C_pfs_table.txt"
digitise(surv_inp,nrisk_inp,km_output = "KM_pembr_pfs0.txt",ipd_output = "IPD_pembr_pfs0.txt")
surv_inp <- "I_pfs_surv.txt"
nrisk_inp <- "I_pfs_table.txt"
digitise(surv_inp,nrisk_inp,km_output = "KM_pembr_pfs1.txt",ipd_output = "IPD_pembr_pfs1.txt")
ipd_files <- list("IPD_pembr_pfs0.txt","IPD_pembr_pfs1.txt")
data <- make.ipd(ipd_files, ctr=1 ,var.labs = c("time","event","arm"))
#fig
max_xlim<-as.integer(max(data$time))+1

#Cox regression and KM curves
fit <- survfit(Surv(time,event==1)~arm,data=data)
res.cox <- coxph(Surv(time, event==1) ~ arm, data =data)
summary(res.cox)
ggsurvplot(fit,data=data,pval = TRUE,risk.table = TRUE,
           palette = c("#D95F02","steelblue"),
           legend.tittle="Intervention",
           legend.labs=c("Control","Intervention"),
           risk.table.height=0.2,
           xlab=c("Time(Months))"),
           break.y.by=0.2,
           break.x.by=2,
           xlim=c(0,max_xlim),
           ylab=c("PFS"),
           risk.table.col="strata",
           tables.theme=theme_cleantable(),
           ggtheme=theme_classic(),
           surv.median.line = "hv")

###5
##OS for JUPITER-06
surv_inp <- "C_os_surv.txt"
nrisk_inp <- "C_os_table.txt"
digitise(surv_inp,nrisk_inp,km_output = "KM_torip_os0.txt",ipd_output = "IPD_torip_os0.txt")
surv_inp <- "I_os_surv.txt"
nrisk_inp <- "I_OS_table.txt"
digitise(surv_inp,nrisk_inp,km_output = "KM_torip_os1.txt",ipd_output = "IPD_torip_os1.txt")
ipd_files <- list("IPD_torip_os0.txt","IPD_torip_os1.txt")
data <- make.ipd(ipd_files, ctr=1 ,var.labs = c("time","event","arm"))
#fig
max_xlim<-as.integer(max(data$time))+1

###Cox regression and KM curves
fit <- survfit(Surv(time,event==1)~arm,data=data)
res.cox <- coxph(Surv(time, event==1) ~ arm, data =data)
summary(res.cox)
ggsurvplot(fit,data=data,pval = TRUE,risk.table = TRUE,
           palette = c("#D95F02","steelblue"),
           legend.tittle="Intervention",
           legend.labs=c("Control","Intervention"),
           risk.table.height=0.2,
           xlab=c("Time(Months))"),
           break.y.by=0.2,
           break.x.by=2,
           xlim=c(0,max_xlim),
           ylab=c("OS"),
           risk.table.col="strata",
           tables.theme=theme_cleantable(),
           ggtheme=theme_classic(),
           surv.median.line = "hv")

###PFS for JUPITER-06
surv_inp <- "C_pfs_surv.txt"
nrisk_inp <- "C_pfs_table.txt"
digitise(surv_inp,nrisk_inp,km_output = "KM_torip_pfs0.txt",ipd_output = "IPD_torip_pfs0.txt")
surv_inp <- "I_pfs_surv.txt"
nrisk_inp <- "I_pfs_table.txt"
digitise(surv_inp,nrisk_inp,km_output = "KM_torip_pfs1.txt",ipd_output = "IPD_torip_pfs1.txt")
ipd_files <- list("IPD_torip_pfs0.txt","IPD_torip_pfs1.txt")
data <- make.ipd(ipd_files, ctr=1 ,var.labs = c("time","event","arm"))
#fig
max_xlim<-as.integer(max(data$time))+1

#Cox regression and KM curves
fit <- survfit(Surv(time,event==1)~arm,data=data)
res.cox <- coxph(Surv(time, event==1) ~ arm, data =data)
summary(res.cox)
ggsurvplot(fit,data=data,pval = TRUE,risk.table = TRUE,
           palette = c("#D95F02","steelblue"),
           legend.tittle="Intervention",
           legend.labs=c("Control","Intervention"),
           risk.table.height=0.2,
           xlab=c("Time(Months))"),
           break.y.by=0.2,
           break.x.by=2,
           xlim=c(0,max_xlim),
           ylab=c("PFS"),
           risk.table.col="strata",
           tables.theme=theme_cleantable(),
           ggtheme=theme_classic(),
           surv.median.line = "hv")

