rm(list = ls())
setwd('/home/data/project/30.JNPJ-6001-1/')
file.dir = '13_ROC_CRC'
if (! dir.exists(file.dir)){
  dir.create(file.dir)
}
setwd(file.dir)

library(tidyverse)
maingene <- read.csv('../06_113model/bestModel.gene.csv')
dat <- read.csv('../00_rawdata/dat.fpkm.tcga.csv',row.names = 1)
dat <- log2(dat + 1)
group<-read.csv('../00_rawdata/group.tcga.csv')
table(group$group)
group$group = factor(group$group, levels = c("Tumor", "Normal"))
hub_exp<-dat[maingene$symbol,]
library(ROCR)
library(ggplot2)

hub_exp2<-t(hub_exp)
library(pROC)
res.auc <- data.frame(symbol = maingene$symbol,AUC = NA)
for (i in c(1:6)) {
  roc<-roc(group$group,hub_exp2[,i],levels=c("Tumor", "Normal"))
  auc <- gsub('Area under the curve: ','',roc$auc,fixed = T)
  res.auc$AUC[i] <- auc
  png(paste0(i, ".", colnames(hub_exp2)[i],".png"),width = 300,height = 300)
  plot(roc,
       print.auc=T,
       print.auc.x=0.4,print.auc.y=0.5,
       print.auc.pattern='AUC=%.2f',
       grid=c(0.5,0.2),
       grid.col=c("black","black"),
       main=paste0(colnames(hub_exp2)[i]),
       col="#E05B5B",
       legacy.axes=T)
  dev.off()
  pdf(paste0(i, ".", colnames(hub_exp2)[i],".pdf"),width = 4,height = 4)
  plot(roc,
       print.auc=T,
       print.auc.x=0.4,print.auc.y=0.5,
       print.auc.pattern='AUC=%.2f',
       grid=c(0.5,0.2),
       grid.col=c("black","black"),
       main=paste0(colnames(hub_exp2)[i]),
       col="#E05B5B",
       legacy.axes=T)
  dev.off()
  i<-i+1
}
res.auc <- data.frame(symbol = maingene$symbol, AUC = NA)
roc_list <- list() 

for (i in 1:6) {
  roc <- roc(group$group, hub_exp2[,i], levels = c("Tumor", "Normal"))
  
  auc <- as.numeric(roc$auc)  
  res.auc$AUC[i] <- auc  # 存AUC到结果表
  
  # ggroc绘图
  p <- ggroc(roc, legacy.axes = TRUE, colour = "#E05B5B", size = 1) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey") +
    annotate("text", x = 0.4, y = 0.5, label = sprintf("AUC=%.2f", auc), size = 4, color = "#E05B5B") +
    ggtitle(colnames(hub_exp2)[i]) +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
  
  roc_list[[i]] <- p  # 存入列表
}

# 拼成2行3列大图
combined_plot <- plot_grid(plotlist = roc_list, nrow = 2, ncol = 3, labels = NULL)

# 保存成pdf或png
ggsave("TCGA_combined_ROC_plot.pdf", combined_plot, width = 9, height = 6)
ggsave("TCGA_combined_ROC_plot.png", combined_plot, width = 9, height = 6, dpi = 300)
