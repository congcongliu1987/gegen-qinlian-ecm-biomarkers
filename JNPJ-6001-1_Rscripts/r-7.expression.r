rm(list = ls()) 
setwd('/home/data/project/30.JNPJ-6001-1/')
file.dir = '07_expression'
if (! dir.exists(file.dir)){
  dir.create(file.dir)
}
setwd(file.dir)
library(tidyverse)
## 训练集-----
dat <- read.csv('../00_rawdata/dat.GSE87466.csv',row.names = 1,check.names = F)
group <- read.csv('../00_rawdata/group.GSE87466.csv')
table(group$group)
hubgene <- read.csv('../06_113model/bestModel.gene.csv')
control.sample <- group$sample[which(group$group=='Normal')]
hub_exp<-dat[hubgene$symbol,]
hub_exp2 <- as.data.frame(hub_exp) %>% na.omit()
hub_exp2$Symbol<-rownames(hub_exp2)
hub_exp2<-gather(hub_exp2,
                 key = sample,
                 value = expr,
                 -c("Symbol"))

hub_exp2$Group<-ifelse(hub_exp2$sample%in%control.sample,'Normal','UC')

hub_exp2$Group <- factor(hub_exp2$Group,levels = c('Normal','UC'))
library(ggplot2)
library(tidyverse)
library(ggpubr)
library(ggsci)
library(rstatix)
stat.test<-hub_exp2%>%
  group_by(Symbol)%>%
  wilcox_test(expr ~ Group)%>%
  adjust_pvalue(method = 'fdr')
stat.test$p<-ifelse(stat.test$p<0.0001,"****",
                    ifelse(stat.test$p<0.001,"***",
                           ifelse(stat.test$p<0.01,"**",
                                  ifelse(stat.test$p<0.05,"*",'ns'))))
colnames(hub_exp2)
exp_plot <- ggplot(hub_exp2,aes(x = Group, y = expr, color = Group)) +
  geom_violin(trim=F,color="black",aes(fill=Group)) + 
  stat_boxplot(geom="errorbar", 
               width=0.1,
               position = position_dodge(0.9)) +
  geom_boxplot(width=0.4,
               position=position_dodge(0.9),
               outlier.shape = NA)+ #绘制箱线图，此处width=0.1控制小提琴图中箱线图的宽窄
  scale_fill_manual(values= c("#1177b0","#bf3553"))+
  scale_colour_manual(values = c("#1177b0","#bf3553"))+
  labs(title="Expression(Train)", x="", y = "",size=20) +
  stat_pvalue_manual(stat.test,
                     y.position = c(10,15,16,12,13,13),
                     size = 3.2,
                     family = "Times",
                     label = "p",
                     #parse = T,
                     face = "bold")+
  theme_bw()+
  theme(plot.title = element_text(hjust =0.5,colour="black",face="bold",size=15),
        axis.text.x=element_text(angle=0,hjust=0.5,colour="black",face="bold",size=12), 
        axis.text.y=element_text(hjust=0.5,colour="black",face="bold",size=12), 
        axis.title.x=element_text(size=16,face="bold"),
        axis.title.y=element_text(size=16,face="bold"),
        legend.text=element_text(face = "bold", hjust = 0.5,colour="black", size=12),
        legend.title = element_text(face = "bold", size = 12),
        legend.position = "top",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  guides(fill='none')+
  facet_wrap(~Symbol,scales = "free",nrow = 2) 
exp_plot
ggsave('01.expression(train).pdf',exp_plot,width = 8,height = 6)
ggsave('01.expression(train).png',exp_plot,width = 8,height = 6,units = 'in',dpi = 600)



## 验证集-----
## 读取数据
library(magrittr)
library(stringr)
library(limma)
library(magrittr)
library(stringr)
library(limma)
dat <- read.csv('../00_rawdata/dat.GSE75214.csv',row.names = 1)
group <- read.csv('../00_rawdata/group.GSE75214.csv')
table(group$group)
control.sample <- group$sample[which(group$group=='Normal')]
hub_exp<-dat[which(rownames(dat) %in% hubgene$symbol),]
hub_exp2 <- as.data.frame(hub_exp) %>% na.omit()
hub_exp2$Symbol<-rownames(hub_exp2)
hub_exp2<-gather(hub_exp2,
                 key = sample,
                 value = expr,
                 -c("Symbol"))

hub_exp2$Group<-ifelse(hub_exp2$sample%in%control.sample,'Normal','UC')

hub_exp2$Group <- factor(hub_exp2$Group,levels = c('Normal','UC'))
library(ggplot2)
library(tidyverse)
library(ggpubr)
library(ggsci)
library(rstatix)
stat.test<-hub_exp2%>%
  group_by(Symbol)%>%
  wilcox_test(expr ~ Group)%>%
  adjust_pvalue(method = 'fdr')
stat.test$p<-ifelse(stat.test$p<0.0001,"****",
                    ifelse(stat.test$p<0.001,"***",
                           ifelse(stat.test$p<0.01,"**",
                                  ifelse(stat.test$p<0.05,"*",'ns'))))
colnames(hub_exp2)
exp_plot <- ggplot(hub_exp2,aes(x = Group, y = expr, color = Group)) +
  geom_violin(trim=F,color="black",aes(fill=Group)) +  
  stat_boxplot(geom="errorbar", 
               width=0.1,
               position = position_dodge(0.9)) +
  geom_boxplot(width=0.4,
               position=position_dodge(0.9),
               outlier.shape = NA)+ #绘制箱线图，此处width=0.1控制小提琴图中箱线图的宽窄
  scale_fill_manual(values= c("#1177b0","#bf3553"))+
  scale_colour_manual(values = c("#1177b0","#bf3553"))+
  labs(title="Expression(Test)", x="", y = "",size=20) +
  stat_pvalue_manual(stat.test,
                     y.position = c(13,17,17,14,13,13),
                     size = 3.2,
                     family = "Times",
                     label = "p",
                     #parse = T,
                     face = "bold")+
  theme_bw()+
  theme(plot.title = element_text(hjust =0.5,colour="black",face="bold",size=15),
        axis.text.x=element_text(angle=0,hjust=0.5,colour="black",face="bold",size=12), 
        axis.text.y=element_text(hjust=0.5,colour="black",face="bold",size=12), 
        axis.title.x=element_text(size=16,face="bold"),
        axis.title.y=element_text(size=16,face="bold"),
        legend.text=element_text(face = "bold", hjust = 0.5,colour="black", size=12),
        legend.title = element_text(face = "bold", size = 12),
        legend.position = "top",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  guides(fill='none')+
  facet_wrap(~Symbol,scales = "free",nrow = 2) 
exp_plot
ggsave('02.expression(test).pdf',exp_plot,width = 8,height = 6)
ggsave('02.expression(test).png',exp_plot,width = 8,height = 6,units = 'in',dpi = 600)

