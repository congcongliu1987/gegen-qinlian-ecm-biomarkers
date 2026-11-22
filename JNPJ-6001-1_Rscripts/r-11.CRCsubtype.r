rm(list = ls()) 
setwd('/home/data/project/30.JNPJ-6001-1/')
file.dir = '11_CRCsubtype'
if (! dir.exists(file.dir)){
  dir.create(file.dir)
}
setwd(file.dir)

dat <- read.csv('../00_rawdata/dat.fpkm.tcga.csv',row.names = 1,check.names = F)
group <- read.csv('../00_rawdata/group.tcga.csv')
table(group$group)
case.sample <- group$sample[which(group$group == 'Tumor')]
dat <- dat[,case.sample]
dat <- log2(dat + 1)
gene <- read.csv('../06_113model/bestModel.gene.csv')
dat.cluster<-dat[gene$symbol,]
cluster_exp<-as.matrix(dat.cluster)
library(ConsensusClusterPlus)
df<-sweep(cluster_exp,1, apply(cluster_exp,1,median,na.rm=T))
maxK <-  5 
results <-  ConsensusClusterPlus(df,
                                 maxK = maxK,
                                 reps = 1000,              
                                 pItem = 0.8,             
                                 pFeature = 1,
                                 clusterAlg = "pam",      
                                 seed = 123456,
                                 title="consensus",
                                 innerLinkage="complete",
                                 plot="pdf")
results <-  ConsensusClusterPlus(df,
                                 maxK = maxK,
                                 reps = 1000,              
                                 pItem = 0.8,           
                                 pFeature = 1,
                                 clusterAlg = "pam",      
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


##PCA------
pca_exp<-t(cluster_exp)
pca1<-prcomp(pca_exp,center = TRUE,scale. = TRUE)
df1<-pca1$x  
df1 <- as.data.frame(df1)
summ1 <- summary(pca1)
summ1
summ1 <- summary(pca1)
summ1
xlab1 <- paste0("PC1(",round(summ1$importance[2,1]*100,2),"%)")
ylab1 <- paste0("PC2(",round(summ1$importance[2,2]*100,2),"%)")
cluster<-as.data.frame(cluster)
cluster$cluster<-ifelse(cluster$cluster=='1','cluster 1',
                        ifelse(cluster$cluster=='2','cluster 2','cluster 3'))
library(ggplot2)
##survival---------
survival1 <- read_tsv('/home/data/t210508/files/database/TCGAsurvival/TCGA-COAD.survival.tsv.gz')
survival2 <- read_tsv('/home/data/t210508/files/database/TCGAsurvival/TCGA-READ.survival.tsv.gz')
survival <- rbind(survival1,survival2)
survival <- survival[,-4]
survival <- survival[which(survival$sample %in% colnames(dat)),]
write.csv(survival,file = 'survival.csv',row.names = F)
cluster_dat<-as.data.frame(t(dat.cluster))
cluster_dat$cluster<-as.vector(cluster$cluster,)
cluster_dat$sample<-rownames(cluster_dat)
cluster_dat<-merge(survival,cluster_dat,by='sample')
rownames(cluster_dat)<-cluster_dat$sample
cluster_dat<-cluster_dat[,-1]
library(survival)
library(survminer)
table(cluster_dat$cluster)
kmfit<-survfit(Surv(OS.time, OS) ~ cluster, data =  cluster_dat)
cluster_survival_median <- ggsurvplot(kmfit,
                                      pval = TRUE, 
                                      conf.int = F,
                                      legend.labs=c("cluster 1","cluster 2"),
                                      legend.title="cluster",
                                      title="Train KM",
                                      font.main = c(15,"bold"),
                                      risk.table = TRUE, 
                                      risk.table.col = "strata", 
                                      linetype = "strata", 
                                      surv.median.line = "hv", 
                                      ggtheme = theme_bw(), 
                                      palette = c("#A73030FF", "#0073C2FF"))
cluster_survival_median
print(cluster_survival_median)
pdf(file = '01.subtype.km.pdf',w = 8,h = 6,onefile = F)
print(cluster_survival_median) 
dev.off()


png(file = '01.subtype.km.png',w = 8,h = 6,units = 'in',res = 600)
print(cluster_survival_median) 
dev.off()

#B ESTIMATE--------
file.dir = '01_estimate'
if (! dir.exists(file.dir)){
  dir.create(file.dir)
}
setwd(file.dir)

table(group$cluster)
colnames(group)[2] <- 'group'
Low.sample<-group$sample[which(group$group=='cluster 1')]
High.sample<-group$sample[which(group$group=='cluster 2')]
dat<-dat[,group$sample]
group_estimate<-group$group%>%as.factor()
design<-model.matrix(~0 + group_estimate)
rownames(design)<-group$sample
colnames(design)<-levels(group_estimate)
design<-as.data.frame(design)
Low<-rownames(design)[which(design$`cluster 1`==1)]
High<-rownames(design)[which(design$`cluster 2`==1)]
length(Low)
length(High)
library(estimate)
expr_train <- dat
write.table(expr_train, 
            'expr.txt', 
            col.names = T, 
            row.names = T, 
            quote = F, sep="\t")
filterCommonGenes(input.f = './expr.txt', 
                  output.f = 'expr_train.gct', 
                  id = 'GeneSymbol')
estimateScore('expr_train.gct', 'train_purity.gct', platform="affymetrix")
es_score <- read.table('train_purity.gct', skip = 2, header = T, check.names = F)
immu_score <- es_score[,3:length(es_score)]
rownames(immu_score) <- es_score$NAME
write.table(es_score,
            file = "es_score.xls",
            sep = "\t",
            quote = F,
            row.names = F)

## 小提琴图--------
violin_dat <- data.frame(t(immu_score))
rownames(violin_dat)<-group$sample
violin_dat$sample <- rownames(violin_dat)
violin_dat$group <- ifelse(violin_dat$sample %in% Low.sample,
                           "Cluster 1", "Cluster 2")
table(violin_dat$group)
library(tidyr)
library(ggplot2)
library(ggpubr)
{
  p1 <- ggplot(violin_dat, aes(x=group,y=StromalScore, fill=group))+
    geom_violin() + #绘制小提琴图
    stat_boxplot(geom="errorbar",
                 width=0.1,
                 position = position_dodge(0.9)) +
    geom_boxplot(width=0.4,
                 position=position_dodge(0.9),
                 outlier.shape = NA,
                 fill='white')+ #绘制箱线图，此处width=0.1控制小提琴图中箱线图的宽窄
    geom_point(aes(fill = group),
               size = 0.5,
               position = position_dodge(0.9))+
    scale_fill_manual(values = c("#A28CC2","seagreen"), name = "Group") + 
    stat_compare_means(aes(group = group),method = 'wilcox')+ 
    theme_bw()+ #背景变为白色
    ylim(min(violin_dat$StromalScore-100),max(violin_dat$StromalScore +100)) +
    theme(panel.background = element_rect(colour = 'black', size = 1, fill = 'white'),
          panel.grid = element_blank()) + 
    theme(axis.text.x=element_text(family="Times",size=12,face="bold"), #设置x轴刻度标签的字体显示倾斜角度为15度，并向下调整1(hjust = 1)，字体簇为Times大小为20
          axis.text.y=element_text(family="Times",size=12,face="bold"), #设置y轴刻度标签的字体簇，字体大小，字体样式为plain
          axis.title.y=element_text(family="Times",size = 15,face="bold"), #设置y轴标题的字体属性
          axis.title.x=element_text(family="Times",size = 10,face="bold"),
          plot.title = element_text(hjust = 0.5,family="Times",size = 15,face="bold"),
          legend.text = element_text(face = "bold"),
          legend.title = element_text(face = "bold"))+
    labs(title="", x="", y="Stromal Score")+
    guides(fill='none')
  p1
  p2 <- ggplot(violin_dat, aes(x=group,y=ImmuneScore, fill=group))+
    geom_violin() + #绘制小提琴图
    stat_boxplot(geom="errorbar",
                 width=0.1,
                 position = position_dodge(0.9)) +
    geom_boxplot(width=0.4,
                 position=position_dodge(0.9),
                 outlier.shape = NA,
                 fill='white')+ #绘制箱线图，此处width=0.1控制小提琴图中箱线图的宽窄
    geom_point(aes(fill = group),
               size = 0.5,
               position = position_dodge(0.9))+
    scale_fill_manual(values = c("#A28CC2","seagreen"), name = "Group") +
    stat_compare_means(aes(group = group),method = 'wilcox')+ 
    theme_bw()+ #背景变为白色
    ylim(min(violin_dat$ImmuneScore-100),max(violin_dat$ImmuneScore +100)) +
    theme(panel.background = element_rect(colour = 'black', size = 1, fill = 'white'),
          panel.grid = element_blank()) + 
    theme(axis.text.x=element_text(family="Times",size=12,face="bold"), #设置x轴刻度标签的字体显示倾斜角度为15度，并向下调整1(hjust = 1)，字体簇为Times大小为20
          axis.text.y=element_text(family="Times",size=12,face="bold"), #设置y轴刻度标签的字体簇，字体大小，字体样式为plain
          axis.title.y=element_text(family="Times",size = 15,face="bold"), #设置y轴标题的字体属性
          axis.title.x=element_text(family="Times",size = 10,face="bold"),
          plot.title = element_text(hjust = 0.5,family="Times",size = 15,face="bold"),     
          legend.text = element_text(face = "bold"),
          legend.title = element_text(face = "bold"))+
    labs(title="", x="", y="Immune Score")+
    guides(fill='none')
  p2
  p3 <- ggplot(violin_dat, aes(x=group, y=ESTIMATEScore, fill=group))+
    geom_violin() + #绘制小提琴图
    stat_boxplot(geom="errorbar",
                 width=0.1,
                 position = position_dodge(0.9)) +
    geom_boxplot(width=0.4,
                 position=position_dodge(0.9),
                 outlier.shape = NA,
                 fill='white')+ #绘制箱线图，此处width=0.1控制小提琴图中箱线图的宽窄
    geom_point(aes(fill = group),
               size = 0.5,
               position = position_dodge(0.9))+
    scale_fill_manual(values = c("#A28CC2","seagreen"), name = "Group") +
    stat_compare_means(aes(group = group),method = 'wilcox')+ 
    theme_bw()+ #背景变为白色
    theme(panel.background = element_rect(colour = 'black', size = 1, fill = 'white'),
          panel.grid = element_blank()) + 
    ylim(min(violin_dat$ESTIMATEScore-100),max(violin_dat$ESTIMATEScore +100)) +
    theme(axis.text.x=element_text(family="Times",size=12,face="bold"), #设置x轴刻度标签的字体显示倾斜角度为15度，并向下调整1(hjust = 1)，字体簇为Times大小为20
          axis.text.y=element_text(family="Times",size=12,face="bold"), #设置y轴刻度标签的字体簇，字体大小，字体样式为plain
          axis.title.y=element_text(family="Times",size = 15,face="bold"), #设置y轴标题的字体属性
          axis.title.x=element_text(family="Times",size = 10,face="bold"),
          plot.title = element_text(hjust = 0.5,family="Times",size = 15,face="bold"),
          legend.text = element_text(face = "bold"),
          legend.title = element_text(face = "bold"))+
    labs(title="", x="", y="ESTIMATE Score")+
    guides(fill='none')
  p3
  p4 <- ggplot(violin_dat, aes(x=group, y=TumorPurity, fill=group))+
    geom_violin() + #绘制小提琴图
    stat_boxplot(geom="errorbar",
                 width=0.1,
                 position = position_dodge(0.9)) +
    geom_boxplot(width=0.4,
                 position=position_dodge(0.9),
                 outlier.shape = NA,
                 fill='white')+ #绘制箱线图，此处width=0.1控制小提琴图中箱线图的宽窄
    geom_point(aes(fill = group),
               size = 0.5,
               position = position_dodge(0.9))+
    scale_fill_manual(values = c("#A28CC2","seagreen"), name = "Group") +
    stat_compare_means(aes(group = group),method = 'wilcox')+
    theme_bw()+ #背景变为白色
    theme(panel.background = element_rect(colour = 'black', size = 1, fill = 'white'),
          panel.grid = element_blank()) +
    ylim(min(violin_dat$TumorPurity-0.05),max(violin_dat$TumorPurity +0.05)) +
    theme(axis.text.x=element_text(family="Times",size=12,face="bold"), #设置x轴刻度标签的字体显示倾斜角度为15度，并向下调整1(hjust = 1)，字体簇为Times大小为20
          axis.text.y=element_text(family="Times",size=12,face="bold"), #设置y轴刻度标签的字体簇，字体大小，字体样式为plain
          axis.title.y=element_text(family="Times",size = 15,face="bold"), #设置y轴标题的字体属性
          axis.title.x=element_text(family="Times",size = 10,face="bold"),
          plot.title = element_text(hjust = 0.5,family="Times",size = 15,face="bold"),
          legend.text = element_text(face = "bold"),
          legend.title = element_text(face = "bold"))+
    labs(title="", x="", y="Tumor Purity")+
    guides(fill='none')
  p4
  library(cowplot)
  p5 <- cowplot::plot_grid(p1,p2,p3,p4,
                           nrow = 1, 
                           align = 'h', 
                           vjust = -0.3) 
  p5
  
  
  
}
ggsave(filename = '01.estimate.all.pdf',p5,w=12,h=4)
ggsave(filename = '01.estimate.all.png',p5,w=12,h=4,units = 'in',dpi = 600,bg = 'white')
setwd('../')

#C TIDE--------
file.dir = '02_TIDE'
if (! dir.exists(file.dir)){
  dir.create(file.dir)
}
setwd(file.dir)

dat.tcga<-read.csv("../../00_rawdata/dat.fpkm.tcga.csv", row.names = 1,check.names = F)

tide_dat <- dat
rownmean <- apply(tide_dat,1,mean)
tide_dat2 <- sweep(tide_dat, 1, rownmean)
dim(tide_dat2)
write.table(tide_dat2,
            file ="tide_dat.txt",
            sep = "\t",
            quote = F,
            row.names = T)


tide_result <- read.csv("tide_results.csv",header = T)
colnames(tide_result)
tide_result2 <- subset(tide_result, select = c("Patient", "TIDE","Dysfunction","Exclusion"))
tide_plot_dat <- data.frame(Patient=group$sample,
                            group=group$group)
tide_plot_dat <- merge(tide_plot_dat, tide_result2, by = "Patient")

dim(tide_plot_dat)
colnames(tide_plot_dat)
library(ggsci)
table(tide_plot_dat$group)
my_comparisons = list( c("cluster 1", "cluster 2"))
p1<-ggplot(tide_plot_dat, aes(x=group, y=Exclusion,fill=group)) + 
  geom_boxplot()+
  stat_compare_means(comparisons = my_comparisons,
                     label = "p.signif",
                     method = "wilcox.test",
                     color='black',family="Times",face = "bold")+
  scale_fill_npg()+
  scale_fill_manual(values=c("#E05B5B","#519CBA")) + 
  labs(x = "", y = "Exclusion Score", title = "") + 
  theme_bw() +
  guides(fill=F) 
p1

p2<-ggplot(tide_plot_dat, aes(x=group, y=Dysfunction,fill=group)) + 
  geom_boxplot()+
  stat_compare_means(comparisons = my_comparisons,
                     label = "p.signif",
                     method = "wilcox.test",
                     color='black',family="Times",face = "bold")+
  scale_fill_npg()+
  scale_fill_manual(values=c("#E05B5B","#519CBA")) + 
  labs(x = "", y = "Dysfunction Score", title = "") + 
  theme_bw() +
  guides(fill=F) 
p2

p3<-ggplot(tide_plot_dat, aes(x=group, y=TIDE,fill=group)) + 
  geom_boxplot()+
  stat_compare_means(comparisons = my_comparisons,
                     label = "p.signif",
                     method = "wilcox.test",
                     color='black',family="Times",face = "bold")+
  scale_fill_npg()+
  scale_fill_manual(values=c("#E05B5B","#519CBA")) + 
  labs(x = "", y = "TIDE score", title = "") + 
  theme_bw() +
  guides(fill=F) 
p3
library(patchwork)
p<-p1+p2+p3+plot_layout(ncol = 3)
p
pdf('01.TIDE_boxplot.pdf',width=7,height=4,family='Times')
print(p)
dev.off()

png('01.TIDE_boxplot.png',width=7,height=4,family='Times',units='in',res=600)
print(p)
dev.off()

setwd('../')
#D IC50--------
file.dir = '03_IC50'
if (! dir.exists(file.dir)){
  dir.create(file.dir)
}
setwd(file.dir)
high.sample<-group$sample[which(group$group=='cluster 1')]
low.sample<-group$sample[which(group$group=='cluster 2')]

library(oncoPredict)
library(data.table)
library(gtools)
library(reshape2)
library(ggpubr)
library(reshape2)
library(ggpubr)

th=theme(axis.text.x = element_text(angle = 45,vjust = 0.5))
dir='../../../../files/IC50/DataFiles/Training Data/'
GDSC2_Expr = readRDS(file=file.path(dir,'GDSC2_Expr (RMA Normalized and Log Transformed).rds'))
GDSC2_Res = readRDS(file = file.path(dir,"GDSC2_Res.rds"))
GDSC2_Res <- exp(GDSC2_Res) 
model_expr<-dat[,group$sample]
testExpr<- as.matrix(model_expr)
dim(testExpr)  

calcPhenotype(trainingExprData = GDSC2_Expr,
              trainingPtype = GDSC2_Res,
              testExprData = testExpr,
              batchCorrect = 'eb',  #   "eb" for ComBat  
              powerTransformPhenotype = TRUE,
              removeLowVaryingGenes = 0.2,
              minNumSamples = 5, 
              printOutput = TRUE, 
              removeLowVaringGenesFrom = 'rawData' )
result <- read.csv('calcPhenotype_Output/DrugPredictions.csv',row.names = 1)
b <- result
removeRowsAllNa  <- function(x){x[apply(x, 1, function(y) any(!is.na(y))),]}
removeColsAllNa  <- function(x){x[, apply(x, 2, function(y) any(!is.na(y)))]}
c<-removeColsAllNa(b)
na_flag <- apply(is.na(c), 2, sum)
x <- c[, which(na_flag == 0)]
dim(x)
medicinal_result <- t(x) 
high_group <- high.sample
low_group <- low.sample
pvalue = padj = log2FoldChange <- matrix(0, nrow(medicinal_result), 1)
for (i in 1:nrow(medicinal_result)){
  pvalue[i, 1] = p.value = wilcox.test(medicinal_result[i, high_group],
                                       medicinal_result[i, low_group])$p.value
  log2FoldChange[i, 1] = mean(medicinal_result[i, high_group]) - 
    mean(medicinal_result[i, low_group])
}
padj <- p.adjust(as.vector(pvalue), "fdr", n = length(pvalue))
rTable <- data.frame(log2FoldChange, 
                     pvalue, 
                     padj,
                     row.names = rownames(medicinal_result))
high_group_res <- signif(apply(medicinal_result[rownames(rTable), high_group], 
                               1,
                               median), 4)
low_group_res <- signif(apply(medicinal_result[rownames(rTable), low_group], 
                              1, 
                              median), 4)
rTable <- data.frame(high_group_res, 
                     low_group_res,
                     rTable[, c("padj", "pvalue", "log2FoldChange")])
rTable$drugs <- rownames(rTable)
rTable$sig <- ifelse(rTable$pvalue < 0.05,
                     ifelse(rTable$pvalue < 0.01, 
                            ifelse(rTable$pvalue < 0.001,
                                   ifelse(rTable$pvalue < 0.0001,
                                          paste(rTable$drugs, "****",  sep = ""),
                                          paste(rTable$drugs, "***", sep = "")),
                                   paste(rTable$drugs, "**", sep = "")),
                            paste(rTable$drugs, "*",  sep = "")), 
                     rTable$drugs)

write.table(rTable,
            file = "drugs_wilcox_test.xls",
            quote = F,
            row.names = F)
DE.drug<-rTable[which(rTable$pvalue<0.05),]
druglist <- data.frame(rownames(rTable))
drug <- c('X5.Fluorouracil_1073','Oxaliplatin_1089','Irinotecan_1088','Lapatinib_1558')

library(tidyr)
library(ggplot2)
library(ggpubr) 
all(rownames(rTable) == rownames(medicinal_result))
drugs_res <- data.frame(drugs=rownames(medicinal_result), medicinal_result, pvalue=rTable$pvalue)
drugs_res <- drugs_res[rownames(drugs_res)%in%drug,]
violin_dat <- gather(drugs_res, key=indivs, value=score, -c("drugs","pvalue"))
violin_dat$indivs <- ifelse(gsub("\\.","-",violin_dat$indivs) %in% high_group,
                            "Cluster 1", "Cluster 2") 
violin_dat$indivs <- factor(violin_dat$indivs, levels = c("Cluster 1", "Cluster 2"))
head(violin_dat)
drugs_hub_boxplot1 <- ggboxplot(violin_dat, x = "indivs", y = "score",
                                color = "indivs", palette = c("#A73030FF", "#0073C2FF"),
                                add = "jitter",
                                short.panel.labs = T,
                                ggtheme = theme_bw()) +
  stat_compare_means(label = "p.signif", label.x = 1.4, vjust = 0.5)
drugs_hub_boxplot <- facet(drugs_hub_boxplot1,
                           facet.by = "drugs",
                           short.panel.labs = T,
                           panel.labs.background = list(fill = "white"),
                           ncol = 4,
                           scales = "free_y") + xlab("") + ylab("IC(50)") +
  theme(panel.grid = element_blank(),
        legend.position = "none",
        axis.title.y = element_text(size = 18, face = "bold"),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10, face = "bold"),
        axis.text.x = element_text(size = 10, face = "bold",angle = 45,hjust = 1),
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
        text = element_text(size = 13, face = "bold"))
drugs_hub_boxplot
ggsave(filename = "01.drugs.plot.pdf", height = 4, width = 8)
ggsave(filename = "01.drugs.plot.png", height = 4, width = 8,units = 'in',dpi = 600)

setwd('../')
#E clinicl--------
file.dir = '04_clinical'
if (! dir.exists(file.dir)){
  dir.create(file.dir)
}
setwd(file.dir)
library(readr)
phenotype1<-read_tsv(file ='/home/data/t210508/files/database/TCGAphenotype/TCGA-COAD.clinical.tsv.gz')
phenotype2<-read_tsv(file ='/home/data/t210508/files/database/TCGAphenotype/TCGA-READ.clinical.tsv.gz')
phenotype <- rbind(phenotype1,phenotype2)
train_phenotype<-data.frame(sample=phenotype$sample,
                            Age=phenotype$age_at_index.demographic,
                            Gender=phenotype$gender.demographic,
                            # Grade=phenotype$tumor_grade.diagnoses,
                            Stage=phenotype$ajcc_pathologic_stage.diagnoses,
                            T.stage=phenotype$ajcc_pathologic_t.diagnoses,
                            N.stage=phenotype$ajcc_pathologic_n.diagnoses,
                            M.stage=phenotype$ajcc_pathologic_m.diagnoses)

train_phenotype<-merge(train_phenotype,survival,by='sample')
train_phenotype2<-train_phenotype
train_phenotype2$OS<-as.numeric(train_phenotype2$OS)
train_phenotype2$OS.time<-as.numeric(train_phenotype2$OS.time)

table(train_phenotype2$T.stage)
train_phenotype2$T.stage<-gsub('a','',train_phenotype2$T.stage)
train_phenotype2$T.stage<-gsub('b','',train_phenotype2$T.stage)
train_phenotype2$T.stage<-gsub('Tis',NA,train_phenotype2$T.stage)

table(train_phenotype2$N.stage)
train_phenotype2$N.stage<-gsub('NX',NA,train_phenotype2$N.stage)
train_phenotype2$N.stage<-gsub('a','',train_phenotype2$N.stage)
train_phenotype2$N.stage<-gsub('b','',train_phenotype2$N.stage)
train_phenotype2$N.stage<-gsub('c','',train_phenotype2$N.stage)

table(train_phenotype2$M.stage)
train_phenotype2$M.stage<-gsub('MX',NA,train_phenotype2$M.stage)
train_phenotype2$M.stage<-gsub('a','',train_phenotype2$M.stage)
train_phenotype2$M.stage<-gsub('b','',train_phenotype2$M.stage)

table(train_phenotype2$Stage)
train_phenotype2$Stage<-gsub('A','',train_phenotype2$Stage)
train_phenotype2$Stage<-gsub('B','',train_phenotype2$Stage)
train_phenotype2$Stage<-gsub('C','',train_phenotype2$Stage)
train_phenotype2$Stage<-gsub('Stage IV','Stage4',train_phenotype2$Stage)
train_phenotype2$Stage<-gsub('Stage III','Stage3',train_phenotype2$Stage)
train_phenotype2$Stage<-gsub('Stage II','Stage2',train_phenotype2$Stage)
train_phenotype2$Stage<-gsub('Stage I','Stage1',train_phenotype2$Stage)


table(train_phenotype2$Age)
median(train_phenotype2$Age)
train_phenotype2$Age <- ifelse(train_phenotype2$Age>65,'>65','<=65')

table(train_phenotype2$Gender)

library(tidyverse)
colnames(train_phenotype2)
sub_risk <- group
train_phenotype3 <- merge(sub_risk,
                          train_phenotype2,
                          by = "sample")

write.csv(train_phenotype3,file = "clinical_subtype.csv",row.names = F,quote = F)

library(ggplot2)
library(dplyr)
library(gridExtra)
colnames(train_phenotype3) <- c('Sample','Group','Age','Gender','Stage','Tstage','Nstage','Mstage','OS.time','OS')
colors <- c('#A28CC2','#4484b1',"#4ea660",'#f7a27f',"#79cafb")

p1 <- ggplot(na.omit(train_phenotype3), aes(x = Group, fill = Gender)) +
  geom_bar(position = "fill") +
  theme_minimal() +
  labs(title = "Gender Distribution by Subtype", x = "Subtype", y = "Proportion", fill = "Stage") +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = colors) +
  theme(legend.position = "top") 

# 绘制Age的堆叠图
p2 <- ggplot(na.omit(train_phenotype3), aes(x = Group, fill = Age)) +
  geom_bar(position = "fill") +
  theme_minimal() +
  labs(title = "Age Distribution by Subtype", x = "Subtype", y = "Proportion", fill = "Age") +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = colors)+
  theme(legend.position = "top") 

p3 <- ggplot(na.omit(train_phenotype3), aes(x = Group, fill = Stage)) +
  geom_bar(position = "fill") +
  theme_minimal() +
  labs(title = "Stage Distribution by Subtype", x = "Subtype", y = "Proportion", fill = "Stage") +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = colors)+
  theme(legend.position = "top")

p4 <- ggplot(na.omit(train_phenotype3), aes(x = Group, fill = Tstage)) +
  geom_bar(position = "fill") +
  theme_minimal() +
  labs(title = "T.stage Distribution by Subtype", x = "Subtype", y = "Proportion", fill = "Tstage") +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = colors)+
  theme(legend.position = "top") 


p5 <- ggplot(na.omit(train_phenotype3), aes(x = Group, fill = Nstage)) +
  geom_bar(position = "fill") +
  theme_minimal() +
  labs(title = "N.stage Distribution by Subtype", x = "Subtype", y = "Proportion", fill = "Nstage") +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = colors)+
  theme(legend.position = "top") 

p6 <- ggplot(na.omit(train_phenotype3), aes(x = Group, fill = Mstage)) +
  geom_bar(position = "fill") +
  theme_minimal() +
  labs(title = "M.stage Distribution by Subtype", x = "Subtype", y = "Proportion", fill = "Mstage") +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = colors)+
  theme(legend.position = "top") 

p7 <- ggplot(na.omit(train_phenotype3), aes(x = Group, fill = OS)) +
  geom_bar(position = "fill") +
  theme_minimal() +
  labs(title = "OS Distribution by Subtype", x = "Subtype", y = "Proportion", fill = "OS") +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = colors)+
  theme(legend.position = "top") 

chisq_gender <- chisq.test(table(train_phenotype3$Group, train_phenotype3$Gender))
chisq_stage <- chisq.test(table(train_phenotype3$Group, train_phenotype3$Stage))
chisq_age <- chisq.test(table(train_phenotype3$Group, train_phenotype3$Age))
chisq_tstage <- chisq.test(table(train_phenotype3$Group, train_phenotype3$Tstage))
chisq_nstage <- chisq.test(table(train_phenotype3$Group, train_phenotype3$Nstage))
chisq_mstage <- chisq.test(table(train_phenotype3$Group, train_phenotype3$Mstage))
chisq_os <- chisq.test(table(train_phenotype3$Group, train_phenotype3$OS))

##打印结果
chisq_gender
chisq_stage
chisq_age
chisq_tstage
chisq_nstage
chisq_mstage
chisq_os

fisher_gender <- fisher.test(table(train_phenotype3$Group, train_phenotype3$Gender))
fisher_age <- fisher.test(table(train_phenotype3$Group, train_phenotype3$Age))
fisher_tstage <- fisher.test(table(train_phenotype3$Group, train_phenotype3$Tstage))
fisher_nstage <- fisher.test(table(train_phenotype3$Group, train_phenotype3$Nstage))
fisher_mstage <- fisher.test(table(train_phenotype3$Group, train_phenotype3$Mstage))
fisher_os <- fisher.test(table(train_phenotype3$Group, train_phenotype3$OS))
fisher_gender
fisher_age
fisher_tstage
fisher_nstage
fisher_mstage
fisher_os


p1 <- p1 + annotate("text", x = 1.5, y = 1.05, label = paste("Chi-sq p-value:", round(chisq_gender$p.value, 4)), size = 4, color = "black")
p2 <- p2 + annotate("text", x = 1.5, y = 1.05, label = paste("Chi-sq p-value:", round(chisq_age$p.value, 4)), size = 4, color = "black")
p3 <- p3 + annotate("text", x = 1.5, y = 1.05, label = paste("Chi-sq p-value:", round(chisq_stage$p.value, 4)), size = 4, color = "black")
p4 <- p4 + annotate("text", x = 1.5, y = 1.05, label = paste("Chi-sq p-value:", round(chisq_tstage$p.value, 4)), size = 4, color = "black")
p5 <- p5 + annotate("text", x = 1.5, y = 1.05, label = paste("Chi-sq p-value:", round(chisq_nstage$p.value, 4)), size = 4, color = "black")
p6 <- p6 + annotate("text", x = 1.5, y = 1.05, label = paste("Chi-sq p-value:", round(chisq_mstage$p.value, 4)), size = 4, color = "black")

p <- grid.arrange(p1, p2,p3,p4,p5,p6,ncol = 3)
print(p)
ggsave(filename = '01.clinical.pdf',p,w = 13, h =8)
ggsave(filename = '01.clinical.png',p,units = 'in',dpi = 600,w = 13, h = 8)

setwd('../')
#F TIICs--------
file.dir = '05_TIICs'
if (! dir.exists(file.dir)){
  dir.create(file.dir)
}
setwd(file.dir)

library(readxl)
library(GSVA)
gene_set <- read.table("../../../../files/TIIC/mmc3.txt",
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

violin.dat <- ssgsea_score %>% as.data.frame()
violin.dat$Type<-rownames(violin.dat)
violin.dat<-gather(violin.dat,
                   key = sample,
                   value = score,
                   -c("Type"))

violin.dat$Group <- group$group[match(violin.dat$sample,group$sample)]
table(violin.dat$Group)
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
        panel.grid.major = element_blank(), #不显示网格线
        panel.grid.minor = element_blank()) #不显示网格线
boxplot_diff
ggsave('01.TIIC.pdf',boxplot_diff,width = 10,height = 6)
ggsave('01.TIIC.png',boxplot_diff,width = 10,height = 6,units = 'in',dpi = 600,bg = 'white')
