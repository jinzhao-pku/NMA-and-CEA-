library(tidyverse)
library(ggplot2)
rm(list = ls())
setwd("D:/project/NMA/CEA/one_way_analysis")
data = read.csv("one_way_data.csv")
PD1 = list()
name = list("cam","niv","sin","tor","pem","tis")
for (i in 1:6) {
  PD1[[i]] = data[c(((i-1)*11+1):((i-1)*11+11)),]
  
}

for (i in 1:6) {
  a = PD1[[i]][order(PD1[[i]]$diff)[4:11],]
  INMB_base <- a[["base.value"]][1] #INMB基础值
  
  # 对INMB的差值绝对值进行排序
  order.var <- a %>% 
    mutate(INMB_diff=abs(upper.limit-lower.limit))%>% 
    arrange(INMB_diff) %>%
    mutate(var=factor(x=var, levels=var)) %>%
    select(var) %>% 
    unlist() %>% 
    levels()
  
  
  width <- 0.45 # 条形的宽度设定 值越大 条形越宽
  
  ## 数据规整
  a2 <- a %>% 
    mutate(INMB_diff=abs(upper.limit-lower.limit))%>% 
    gather(key='type', value='output.value', lower.limit:upper.limit) %>%
    select(var, type, output.value,INMB_diff) %>%
    mutate(var=factor(var, levels=order.var),
           ymin=pmin(output.value, INMB_base),
           ymax=pmax(output.value, INMB_base),
           xmin=as.numeric(var)-width/2,
           xmax=as.numeric(var)+width/2)
  
  m = paste("INMB",name[[i]])
  
  p = ggplot() + 
    geom_rect(data = a2, 
              aes(ymax=ymax, 
                  ymin=ymin, 
                  xmax=xmax, 
                  xmin=xmin, 
                  fill=type)) +
    theme_bw() +
    theme(axis.title.y=element_blank(), 
          legend.position = c(0.85,0.15),
          legend.title = element_blank()) + 
    geom_hline(yintercept = INMB_base) +
    scale_x_continuous(breaks = c(1:length(order.var)), 
                       labels = order.var) +
    coord_flip()+
    scale_fill_manual(values = c("#5C8286",  "#BFBFBF"))+
    labs(x="", y=m)
  ggsave(paste0("plot_", name[[i]], ".tiff"), p,width = 8, height = 6,dpi = 150)
}
