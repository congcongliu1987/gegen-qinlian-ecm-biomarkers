rm(list = ls())
setwd('/home/data/project/30.JNPJ-6001-1/')
file.dir = '14_KM_CRC'
if (! dir.exists(file.dir)){
  dir.create(file.dir)
}
setwd(file.dir)
gene <- read.csv('../06_113model/bestModel.gene.csv')
library(tidyverse)
dat<-read.csv('../00_rawdata/dat.fpkm.tcga.csv',row.names = 1,check.names = F)
km.dat <- dat[gene$symbol,]
km.dat <- log2(km.dat+1)

survival <- read.csv('../11_CRCsubtype/survival.csv')%>%column_to_rownames(var = 'sample')
dat.km <- t(km.dat[,colnames(km.dat)%in%rownames(survival)])%>%as.data.frame()
survival <- survival[rownames(dat.km),]
dat.km <- cbind(dat.km,survival)
dat.km$OS <- as.numeric(dat.km$OS)
dat.km$OS.time <- as.numeric(dat.km$OS.time)
library(survival)
library(survminer)
res <- data.frame(symbol = gene$symbol,pvalue = NA)
for (i in c(1:6)) {
  library(survival)
  dat.km$group <- factor(ifelse(dat.km[,i]>median(dat.km[,i]),'High','Low'),levels = c('High','Low'))
  kmfit <- survfit(Surv(OS.time,OS) ~ group, data = dat.km)
  res$pvalue[i] <- surv_pvalue(kmfit)$pval
  KM <- ggsurvplot(kmfit,
                   pval = TRUE, 
                   conf.int = F,
                   legend.labs=c("High","Low" ),
                   legend.title="group",
                   title=colnames(dat.km[i]),
                   font.main = c(15,"bold"),
                   risk.table = TRUE, 
                   risk.table.col = "strata", 
                   linetype = "strata", 
                   surv.median.line = "hv", 
                   ggtheme = theme_bw(), 
                   palette = c("#A73030FF", "#0073C2FF"))
  print(KM)
  fn1 = paste0(i,'.',colnames(dat.km)[i],'.pdf')
  fn2 = paste0(i,'.',colnames(dat.km)[i],'.png')
  pdf(file = fn1,w = 6, h=5,onefile = F)
  print(KM)
  dev.off()
  png(file = fn2,w = 6, h=5,units = 'in',res = 600)
  print(KM)
  dev.off()
  i <- i+1
}
