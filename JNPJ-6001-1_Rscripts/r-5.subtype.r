rm(list = ls())
setwd('/home/data/project/30.JNPJ-6001-1/')
file.dir = '05_subtype'
if (! dir.exists(file.dir)){
  dir.create(file.dir)
}
setwd(file.dir)

dat <- read.csv('../00_rawdata/dat.GSE87466.csv',row.names = 1,check.names = F)
group <- read.csv('../00_rawdata/group.GSE87466.csv')
table(group$group)
case.sample <- group$sample[which(group$group == 'UC')]
dat <- dat[,case.sample]

gene <- read.csv('../04_candigene/candigene.csv')
dat.cluster<-dat[gene$symbol,]
cluster_exp<-as.matrix(dat.cluster)
library(ConsensusClusterPlus)
df<-sweep(cluster_exp,1, apply(cluster_exp,1,median,na.rm=T))
maxK <-  9 #最多分成几组
results <-  ConsensusClusterPlus(df,
                                 maxK = maxK,
                                 reps = 1000,              # 抽样次数(一般1000或更多)
                                 pItem = 0.8,             # 抽样比例
                                 pFeature = 1,
                                 clusterAlg = "pam",      # 聚类方法
                                 seed = 123456,
                                 title="consensus",
                                 innerLinkage="complete",
                                 plot="pdf")
results <-  ConsensusClusterPlus(df,
                                 maxK = maxK,
                                 reps = 1000,              # 抽样次数(一般1000或更多)
                                 pItem = 0.8,             # 抽样比例
                                 pFeature = 1,
                                 clusterAlg = "pam",      # 聚类方法
                                 seed = 123456,
                                 title="consensus",
                                 innerLinkage="complete",
                                 plot="png")

icl = calcICL(results,
              title="consensus",
              plot="pdf")
cluster<-results[[2]]$consensusClass
cluster
library(tidyverse)
group <- data.frame(cluster)%>%rownames_to_column(var = 'sample')
group$cluster <- ifelse(group$cluster=='1','cluster 1','cluster 2')
write.csv(group,file = 'cluster.csv',row.names = F)
library(readxl)
library(GSVA)
gene_set <- read.table("../../../files/TIIC/mmc3.txt",
                       header = T,
                       sep ="\t")
dat.final2 <- as.matrix(dat)
gene_list <- split(as.matrix(gene_set)[,1],
                   gene_set[,2])
ssgsea.dat <- ssgseaParam(
  dat.final2,
  gene_list,
  assay = NA_character_,
  annotation = NA_character_,
  minSize = 1,
  maxSize = Inf,
  alpha = 0.25,
  normalize = TRUE
)
ssgsea_score = gsva(ssgsea.dat)
write.csv(ssgsea_score,file = "tiics_result.csv",row.names = T)

##加小提琴图
violin.dat <- ssgsea_score %>% as.data.frame()
violin.dat$Type<-rownames(violin.dat)
violin.dat<-gather(violin.dat,
                   key = sample,
                   value = score,
                   -c("Type"))

## 样本分组
violin.dat$Group <- group$cluster[match(violin.dat$sample,group$sample)]
table(violin.dat$Group)
##分面图形
library(ggplot2)
library(tidyverse)
library(ggpubr)
library(ggsci)
library(rstatix)
colnames(violin.dat)
stat.test<-violin.dat%>%
  group_by(Type) %>%
  wilcox_test(score ~ Group)%>%
  adjust_pvalue(method = 'fdr')
write.csv(stat.test,file = 'res.tiics.wilcox.csv',row.names = F,quote = F)
colnames(violin.dat)
de.cell <- stat.test[which(stat.test$p < 0.05),]
length(unique(de.cell$Type))
##23
violin.dat <- violin.dat[which(violin.dat$Type %in% unique(de.cell$Type)),]
boxplot_diff <- ggplot(violin.dat, aes(x=Type, 
                                       y=score,
                                       fill=Group)) +
  stat_boxplot(geom="errorbar",
               width=0.1,
               position = position_dodge(0.9)) +
  geom_boxplot(width=0.7,
               position=position_dodge(0.9),
               outlier.shape = NA)+ #绘制箱线图，此处width=0.1控制小提琴图中箱线图的宽窄
  geom_point(aes(fill = Group),
             size = 0.05,
             position = position_dodge(0.9))+
  scale_fill_manual(values= c("#519CBA", "#E05B5B"))+ #设置填充的颜色
  labs(title="", x="", y = "Score",size=20) +
  stat_compare_means(data = violin.dat,
                     mapping = aes(group = Group),
                     label ="p.signif",
                     hide.ns = F) +
  theme_bw()+#把背景设置为白底
  theme(plot.title = element_text(hjust =0.5,colour="black",face="bold",size=18), # 将图表标题居中
        axis.text.x=element_text(angle=65,hjust=1,colour="black",face="bold",size=10), #设置x轴刻度标签的字体显示倾斜角度为45度，并向下调整1(hjust = 1)，字体大小为14
        axis.text.y=element_text(hjust=0.5,colour="black",size=12), #设置y轴刻度标签的字体簇，字体大小，字体样式为plain
        axis.title.x=element_text(size=16,face="bold"),#设置x轴标题的字体属性
        axis.title.y=element_text(size=14,face="bold"), #设置y轴标题的字体属性
        legend.text=element_text(face="bold", hjust = 0.5,colour="black", size=11), #设置图例的子标题的字体属性
        legend.title=element_text(face="bold", colour="black", size=11),#设置图例的总标题的字体属性
        #legend.justification=c(-0.1,1.2), #可调整图例的位置。##(1,1)第一个1是调整图例在图的内外(左右移动)，第二个1是在图中上下移动。
        #legend.position=c(0, 1.04), #legend.position=c(0,1)左上角，(1,1)是在右上角。
        panel.grid.major = element_blank(), #不显示网格线
        panel.grid.minor = element_blank()) #不显示网格线
boxplot_diff
ggsave('01.TIIC.pdf',boxplot_diff,width = 10,height = 6)
ggsave('01.TIIC.png',boxplot_diff,width = 10,height = 6,units = 'in',dpi = 600,bg = 'white')

#代谢通路-------
geneset <- read.delim2('../../24.JNPJ-3020-12/00_rawdata/12894_2023_1317_MOESM1_ESM.xls',skip = 1,header = F,row.names = 1)
library(tidyr)

geneset$Pathway <- rownames(geneset)

# 将矩阵转成长表格
long_data <- geneset %>%
  pivot_longer(
    cols = -Pathway,      # 排除 Pathway 列，转长表格
    names_to = "seq",    # 新列保存基因名
    values_to = "Symbol" # 新列保存值（1或0）
  ) %>%
  filter(!Symbol == "") %>%  # 筛选属于该通路的基因
  dplyr::select(Symbol, Pathway)     # 仅保留基因名和通路名

# 查看结果
print(long_data)
gene_set <- long_data
dat.final2 <- as.matrix(dat)
gene_list <- split(as.matrix(gene_set)[,1],
                   gene_set[,2])
gseaparam <- ssgseaParam(
  dat.final2,
  gene_list,
  assay = NA_character_,
  annotation = NA_character_,
  minSize = 1,
  maxSize = Inf,
  alpha = 0.25,
  normalize = TRUE
)

intrest_score = gsva(gseaparam)
write.csv(intrest_score,file = "Metabolism.score.csv",row.names = T)

library(tidyverse)
violin.dat <- intrest_score %>% as.data.frame()
violin.dat$Type<-rownames(violin.dat)
violin.dat<-gather(violin.dat,
                   key = sample,
                   value = score,
                   -c("Type"))

## 样本分组
violin.dat$Group <- group$cluster[match(violin.dat$sample,group$sample)]
violin.dat$sample <- rownames(violin.dat)
table(violin.dat$Group)
violin.dat$Group <- factor(violin.dat$Group,levels = c('cluster 1','cluster 2'))
table(violin.dat$Group)
##分面图形
library(ggplot2)
library(tidyverse)
library(ggpubr)
library(ggsci)
library(rstatix)
colnames(violin.dat)
stat.test<-violin.dat%>%
  group_by(Type) %>%
  wilcox_test(score ~ Group)%>%
  adjust_pvalue(method = 'fdr')
write.csv(stat.test,file = 'metabolism.res.wilcox.csv',row.names = F,quote = F)

de.res <- stat.test[which(stat.test$p < 0.05),]
##37
de.score <- intrest_score[de.res$Type,]
group <- group[order(group$cluster),]
de.score <- de.score[,group$sample]
rownames(group) <- NULL
group.anno <- group %>% column_to_rownames(var = 'sample')
group_colors <- list("cluster" = c("cluster 1" = "#104680", "cluster 2" = "#FF7256"))
library(pheatmap)
pdf(file = '02.metabolism.heatmap.pdf',w=12,h=7)
pheatmap(de.score, 
         scale = "row",  # 行标准化
         clustering_distance_rows = "euclidean",  # 行聚类方法
         clustering_distance_cols = "euclidean",  # 列聚类方法
         show_rownames = TRUE,  # 显示行名
         show_colnames = F,  # 显示列名
         annotation_col = group.anno,  # 根据group列注释列
         annotation_colors = group_colors,
         cluster_cols = F,
         color = colorRampPalette(c("#104680", "white", "#FF7256"))(100))  # 热图配色
dev.off()

png(file = '02.metabolism.heatmap.png',w=12,h=7,units = 'in',res = 600,bg = 'white')
pheatmap(de.score, 
         scale = "row",  # 行标准化
         clustering_distance_rows = "euclidean",  # 行聚类方法
         clustering_distance_cols = "euclidean",  # 列聚类方法
         show_rownames = TRUE,  # 显示行名
         show_colnames = F,  # 显示列名
         annotation_col = group.anno,  # 根据group列注释列
         annotation_colors = group_colors,
         cluster_cols = F,
         color = colorRampPalette(c("#104680", "white", "#FF7256"))(100))  # 热图配色
dev.off()
## PCD-------
geneset <- readxl::read_xlsx('../../../files/Geneset/pcd18.xlsx')
library(GSVA)
dat.final2 <- as.matrix(dat)
gene_list <- split(as.matrix(geneset)[,1],
                   geneset[,2])
param <- ssgseaParam(
  dat.final2,
  gene_list,
  assay = NA_character_,
  annotation = NA_character_,
  minSize = 1,
  maxSize = Inf,
  alpha = 0.25,
  normalize = TRUE
)
interest_score <- gsva(param)

write.csv(interest_score,file = "PCDscore.csv",row.names = T)
library(tidyverse)
##加小提琴图
violin.dat <- interest_score %>% as.data.frame()
violin.dat$Type<-rownames(violin.dat)
violin.dat<-gather(violin.dat,
                   key = sample,
                   value = score,
                   -c("Type"))

## 样本分组
violin.dat$Group <- group$cluster[match(violin.dat$sample,group$sample)]
table(violin.dat$Group)
violin.dat$sample <- rownames(violin.dat)
violin.dat$Group <- factor(violin.dat$Group,levels = c('cluster 1','cluster 2'))
table(violin.dat$Group)
##分面图形
library(ggplot2)
library(tidyverse)
library(ggpubr)
library(ggsci)
library(rstatix)
colnames(violin.dat)
stat.test<-violin.dat%>%
  group_by(Type) %>%
  wilcox_test(score ~ Group)%>%
  adjust_pvalue(method = 'fdr')
write.csv(stat.test,file = 'pcd.res.wilcox.csv',row.names = F,quote = F)

de.res <- stat.test[which(stat.test$p < 0.05),]
##12
violin.dat <- violin.dat[which(violin.dat$Type%in%de.res$Type),]
boxplot_diff <- ggplot(violin.dat, aes(x=Type, 
                                       y=score,
                                       fill=Group)) +
  stat_boxplot(geom="errorbar",
               width=0.1,
               position = position_dodge(0.9)) +
  geom_boxplot(width=0.7,
               position=position_dodge(0.9),
               outlier.shape = NA)+ #绘制箱线图，此处width=0.1控制小提琴图中箱线图的宽窄
  geom_point(aes(fill = Group),
             size = 0.05,
             position = position_dodge(0.9))+
  scale_fill_manual(values= c("#519CBA", "#E05B5B"))+ #设置填充的颜色
  labs(title="", x="", y = "Score",size=20) +
  stat_compare_means(data = violin.dat,
                     mapping = aes(group = Group),
                     label ="p.signif",
                     hide.ns = F) +
  theme_bw()+#把背景设置为白底
  theme(plot.title = element_text(hjust =0.5,colour="black",face="bold",size=18), # 将图表标题居中
        axis.text.x=element_text(angle=65,hjust=1,colour="black",face="bold",size=10), #设置x轴刻度标签的字体显示倾斜角度为45度，并向下调整1(hjust = 1)，字体大小为14
        axis.text.y=element_text(hjust=0.5,colour="black",size=12), #设置y轴刻度标签的字体簇，字体大小，字体样式为plain
        axis.title.x=element_text(size=16,face="bold"),#设置x轴标题的字体属性
        axis.title.y=element_text(size=14,face="bold"), #设置y轴标题的字体属性
        legend.text=element_text(face="bold", hjust = 0.5,colour="black", size=11), #设置图例的子标题的字体属性
        legend.title=element_text(face="bold", colour="black", size=11),#设置图例的总标题的字体属性
        panel.grid.major = element_blank(), #不显示网格线
        panel.grid.minor = element_blank()) #不显示网格线
boxplot_diff
ggsave('03.PCD.pdf',boxplot_diff,w=8,h=6)
ggsave('03.PCD.png',boxplot_diff,w=8,h=6,units = 'in',dpi = 600,bg = 'white')
