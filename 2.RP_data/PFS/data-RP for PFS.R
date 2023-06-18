rm(list=ls())
library("stringr")
library(survHE)
library(survival)
library(survminer)
library("dplyr")
setwd("C:/Users/28629/Documents/OS")
filename<- read.delim("Filename-list.txt",header=F)
filename$trt<-c("che","cam+che","che","niv+che","che","tori+che","che","pemb+che","che","tisle+che","che","sinti+che")
filename$study <- c("ESCORT-1st","ESCORT-1st","CheckMate-648","CheckMate-648","Jupiter-06","Jupiter-06","KEYNOTE-590","KEYNOTE-590","RATIONALE-306","RATIONALE-306","ORIENT-15","ORIENT-15")
long<-as.integer(length(filename[[1]]))
IPDfilename<-data.frame(matrix(nrow=long,ncol=1))
for (i in 1:long) {
  temp<-filename[i,1]
  temp<-str_sub(temp,-nchar(temp),-5)
  IPDfilename[i,1]<-temp
  data_temp<-read.table(filename[i,1],header=T,quote = "\"")
  data_temp$patid<-seq(1,length(data_temp$event))
  data_temp$trt<-filename[i,2]
  data_temp$study <- case_when(i==1|i==2~"ESCORT-1st",
                               i==3|i==4~"CheckMate-648",
                               i==5|i==6~"Jupiter-06",
                               i==7|i==8~"KEYNOTE-590",
                               i==9|i==10~"RATIONALE-306",
                               i==11|i==12~"ORIENT-15"
                               )
  data_temp$studycode <- case_when(i==1|i==2~"1",
                               i==3|i==4~"2",
                               i==5|i==6~"3",
                               i==7|i==8~"4",
                               i==9|i==10~"5",
                               i==11|i==12~"6")
  data_temp$txcode <- case_when(data_temp$trt=="che"~"1",
                                data_temp$trt=="cam+che"~"2",
                                data_temp$trt=="niv+che"~"3",
                                data_temp$trt=="tori+che"~"4",
                                data_temp$trt=="pemb+che"~"5",
                                data_temp$trt=="tisle+che"~"6",
                                data_temp$trt=="sinti+che"~"7",
                                )
  assign(IPDfilename[i,1],data_temp)
  write.table(data_temp,filename[i,1],row.names = FALSE)
}
combined_IPD<-rbind(IPD_camel_pfs0,IPD_camel_pfs1,IPD_nivol_pfs0,IPD_nivol_pfs1,IPD_torip_pfs0,IPD_torip_pfs1,IPD_pembr_pfs0,IPD_pembr_pfs1,IPD_tisle_pfs0,IPD_tisle_pfs1,IPD_sinti_pfs0,IPD_sinti_pfs1)
combined_IPD$i<-seq(1,length(combined_IPD$event))
combined_IPD <- combined_IPD%>%select('i','patid',everything())
write.csv(combined_IPD,"RP_PFS.csv",row.names = FALSE)
