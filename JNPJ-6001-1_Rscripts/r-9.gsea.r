rm(list = ls()) 
setwd('/home/data/project/30.JNPJ-6001-1/')
file.dir = '09_GSEA'
if (! dir.exists(file.dir)){
  dir.create(file.dir)
}
setwd(file.dir)
library(data.table)
library(org.Hs.eg.db)
library(clusterProfiler)
library(biomaRt)
library(enrichplot)
library(limma)
library(tidyverse)

dat <- read.csv('../00_rawdata/dat.GSE87466.csv',row.names = 1)
group <- read.csv('../00_rawdata/group.GSE87466.csv')
hubgene <- read.csv('../06_113model/bestModel.gene.csv')
hub.exp<-dat[hubgene$symbol,]
hub.exp<-hub.exp[order(rownames(hub.exp),decreasing = F),]  
kegggmt<-read.gmt("../../../files/GSVA/c2.cp.kegg.v7.5.1.symbols.gmt") 

# 每个hub基因均分析一遍 ------------------------------------------------------------
for (i in c(1:6)){
  train<-t(hub.exp)
  group=as.vector(ifelse(train[,i]>median(train[,i]),'high','low'))  
  group <- factor(group,levels = c("high","low")) 
  condition<-data.frame(group)
  rownames(condition)<-rownames(train)
  design<-model.matrix(~0+group)   
  contrast.matrix<-makeContrasts(contrasts = "grouphigh-grouplow", levels = design) 
  fit <- lmFit(dat,design)     
  fit1 <- contrasts.fit(fit, contrast.matrix)    
  fit2 <- eBayes(fit1)  
  tempOutput <- topTable(fit2, coef=1, n=nrow(fit2),adjust="fdr") 
  genelist<-data.frame(SYMBOL=rownames(tempOutput ),logFC=tempOutput $logFC)
  #开始ID转换
  gene=bitr(genelist$SYMBOL,fromType="SYMBOL",toType="ENTREZID",OrgDb="org.Hs.eg.db") 
  gene <- dplyr::distinct(gene,SYMBOL,.keep_all=TRUE)
  gene_df <- data.frame(ENTREZID=gene$ENTREZID , 
                        SYMBOL = gene$SYMBOL)  
  gene_df <- merge(gene_df,genelist,by="SYMBOL")
  geneList<-gene_df$logFC  
  names(geneList)=gene_df $SYMBOL 
  geneList=sort(geneList,decreasing = T)  
  set.seed(123)
  KEGG<-GSEA(geneList,TERM2GENE = kegggmt,pvalueCutoff = 1)  
  sortKEGG<-data.frame(KEGG)
  sortKEGG <- sortKEGG[which(sortKEGG$pvalue < 0.05),]
  sortKEGG<-sortKEGG[order(sortKEGG$p.adjust, decreasing = F),]
  write.table(sortKEGG,file=paste('0',i,'.',rownames(hub.exp)[i],'.xls',sep=''),sep="\t",quote=F,row.names=F)
  paths <- rownames(sortKEGG[c(1:5),])
  p = gseaplot2(KEGG,paths,base_size =10,color = c('#A28CC2','#104680','#5CB0C3','#F5BE8F','#E05B5B'),rel_heights = c(1.5, 0.3, 0.5),title=rownames(hub.exp)[i])
  fn1 = paste0("0", i, ".", rownames(hub.exp)[i], ".png")
  fn2 = paste0("0", i, ".", rownames(hub.exp)[i],".pdf")
  png(fn1,w=8,h=6,units = "in",res = 600)
  print(p)
  dev.off()
  pdf(fn2,w=8,h=6)
  print(p)
  dev.off()
}
