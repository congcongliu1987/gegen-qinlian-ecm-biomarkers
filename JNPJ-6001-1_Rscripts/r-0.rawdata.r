rm(list = ls())
setwd('/home/data/project/30.JNPJ-6001-1/')
file.dir = '00_rawdata'
if (! dir.exists(file.dir)){
  dir.create(file.dir)
}
setwd(file.dir)
library(data.table)
library(tidyverse)
library(GEOquery)
library(Biobase)
# GSE87466-----
gset<-getGEO("GSE87466",
             destdir = '.',
             GSEMatrix = T,
             getGPL = F)
expr<-as.data.frame(exprs(gset[[1]]))
gpl<-getGEO("GPL13158",destdir = '.')
a=gset[[1]]
gpl<-Table(gpl)
colnames(gpl)
probe2symobl<-gpl %>%
  dplyr::select('ID','Gene Symbol')%>%
  filter('Gene Symbol'!='')%>%
  separate('Gene Symbol',c('symbol','drop'),sep = '///')%>%
  dplyr::select(-drop)
probe2symobl=probe2symobl[probe2symobl$symbol!='',]
dat<-expr
dat$ID<-rownames(dat)
dat$ID<-as.character(dat$ID)
probe2symobl$ID<-as.character(probe2symobl$ID)
dat<-dat %>%
  inner_join(probe2symobl,by='ID')%>%
  dplyr::select(-ID)%>%     
  dplyr::select(symbol,everything())%>%     
  mutate(rowMean=rowMeans(.[grep('GSM',names(.))]))%>%    
  arrange(desc(rowMean))%>%       
  distinct(symbol,.keep_all = T)%>%     
  dplyr::select(-rowMean)%>%    
  tibble::column_to_rownames(colnames(.)[1])   
pd<-pData(a)
boxplot(dat,las = 2)
class(dat)
range(dat)
group <- data.frame(sample=pd$geo_accession,group=pd$`disease:ch1`)
table(group$group)
group$group <- ifelse(group$group == 'Ulcerative Colitis (UC)','UC','Normal')
dat <- dat[,group$sample]

write.csv(dat,file = 'dat.GSE87466.csv')
write.csv(group,file = 'group.GSE87466.csv',row.names = F)


# GSE75214-----
gset<-getGEO("GSE75214",
             destdir = '.',
             GSEMatrix = T,
             getGPL = F)
expr<-as.data.frame(exprs(gset[[1]]))
library(AnnoProbe)
probe2symobl <- idmap(gpl = 'GPL6244')
colnames(probe2symobl) <- c('ID','symbol')
dat<-expr
dat$ID<-rownames(dat)
dat$ID<-as.character(dat$ID)
probe2symobl$ID<-as.character(probe2symobl$ID)
dat<-dat %>%
  inner_join(probe2symobl,by='ID')%>%
  dplyr::select(-ID)%>%     
  dplyr::select(symbol,everything())%>%    
  mutate(rowMean=rowMeans(.[grep('GSM',names(.))]))%>%   
  arrange(desc(rowMean))%>%      
  distinct(symbol,.keep_all = T)%>%     
  dplyr::select(-rowMean)%>%   
  tibble::column_to_rownames(colnames(.)[1])  
pd<-pData(a)
boxplot(dat,las = 2)
class(dat)
range(dat)
group <- data.frame(sample=pd$geo_accession,group=pd$`disease:ch1`)
table(group$group)
group <- subset(group,group %in% c('control','ulcerative colitis'))
group$group <- ifelse(group$group == 'ulcerative colitis','UC','Normal')
dat <- dat[,group$sample]

write.csv(dat,file = 'dat.GSE75214.csv')
write.csv(group,file = 'group.GSE75214.csv',row.names = F)


#TCGA COAD +READ ------
library(TCGAbiolinks)
library(readr)
library(readxl)
library(tidyverse)
tcga.expr1<-read_tsv(file = '/home/data/t210508/files/database/TCGAmatrix/TCGA-COAD.star_counts.tsv.gz')%>%as.data.frame()%>%
  column_to_rownames(var = 'Ensembl_ID')
tcga.expr2<-read_tsv(file = '/home/data/t210508/files/database/TCGAmatrix/TCGA-READ.star_counts.tsv.gz')%>%as.data.frame()%>%
  column_to_rownames(var = 'Ensembl_ID')
tcga.expr1 <- tcga.expr1[rownames(tcga.expr2),]
tcga.expr <- cbind(tcga.expr1,tcga.expr2)

tcga.expr<-2^tcga.expr-1
genecode<-read.table(file = '/home/data/files/GENEANNO/gencode.v36.annotation.gtf.gene.probemap')
probe2symbol<-genecode[,(1:2)]
colnames(probe2symbol)<-c('ID','symbol')
probe2symbol<-probe2symbol[-1,]
dat.tcga<-tcga.expr
dat.tcga$ID <- rownames(dat.tcga)
dat.tcga$ID<-as.character(dat.tcga$ID)
probe2symbol$ID<-as.character(probe2symbol$ID)
dat.tcga<-dat.tcga %>%
  inner_join(probe2symbol,by='ID')%>% 
  select(-ID)%>%     
  select(symbol,everything())%>%    
  mutate(rowMean=rowMeans(.[grep('TCGA',names(.))]))%>%  
  arrange(desc(rowMean))%>%       
  distinct(symbol,.keep_all = T)%>%      
  select(-rowMean)%>% 
  tibble::column_to_rownames(colnames(.)[1])  
dim(dat.tcga)

mete=data.frame(colnames(dat.tcga))
for (i in 1:length(mete[,1])) {
  num=as.numeric(as.character(substring(mete[i,1],14,15)))
  if(num %in% seq(1,9)){mete[i,2]="T"}
  if(num %in% seq(10,29)){mete[i,2]="N"}
}
names(mete)=c("id","group")
table(mete$group)
mete$group=as.factor(mete$group)
mete=subset(mete,mete$group=="T")
exp_tumor<-dat.tcga[,which(colnames(dat.tcga)%in%mete$id)] %>% as.data.frame()
exp_control <- dat.tcga[,which(!colnames(dat.tcga) %in% mete$id)] %>% as.data.frame()
survival1 <- read_tsv('/home/data/files/database/TCGAsurvival/TCGA-COAD.survival.tsv.gz')
survival2 <- read_tsv('/home/data/files/database/TCGAsurvival/TCGA-READ.survival.tsv.gz')
survival <- rbind(survival1,survival2)
exp_tumor <- exp_tumor[,which(colnames(exp_tumor)%in%survival$sample)]
dat.final <- cbind(exp_control,exp_tumor)
pcg <- read.csv('../../../files/GENEANNO/pcg.v36.csv')
length(unique(pcg$gene_name))
colnames(pcg)
dat.final <- dat.final[pcg$gene_name,] %>% na.omit()
write.csv(dat.final,file = 'dat.count.tcga.csv',row.names = T,quote = F)
group <- data.frame(sample = colnames(dat.final), group = c(rep('Normal',51),rep('Tumor',607)))
write.csv(group,file = 'group.tcga.csv',row.names = F,quote = F)

##fpkm-----
expr_fpkm1<-read_tsv(file = '/home/data/files/database/TCGAmatrix/TCGA-COAD.star_fpkm.tsv.gz')%>%as.data.frame()%>%
  column_to_rownames(var = 'Ensembl_ID')
expr_fpkm2<-read_tsv(file = '/home/data/files/database/TCGAmatrix/TCGA-READ.star_fpkm.tsv.gz')%>%as.data.frame()%>%
  column_to_rownames(var = 'Ensembl_ID')
expr_fpkm1 <- expr_fpkm1[rownames(expr_fpkm2),]
expr_fpkm <- cbind(expr_fpkm1,expr_fpkm2)

expr_fpkm<-2^expr_fpkm-1
dat_fpkm<-expr_fpkm
dat_fpkm$ID <- rownames(dat_fpkm)
dat_fpkm$ID<-as.character(dat_fpkm$ID)
probe2symbol$ID<-as.character(probe2symbol$ID)
dat_fpkm<-dat_fpkm %>%
  inner_join(probe2symbol,by='ID')%>% 
  select(-ID)%>%    
  select(symbol,everything())%>%     
  mutate(rowMean=rowMeans(.[grep('TCGA',names(.))]))%>%    
  arrange(desc(rowMean))%>%      
  distinct(symbol,.keep_all = T)%>%     
  select(-rowMean)%>%     
  tibble::column_to_rownames(colnames(.)[1])  
dat_fpkm<-dat_fpkm[,colnames(dat.final)]
dat_fpkm<-dat_fpkm[rownames(dat.final),] %>% na.omit()

write.csv(dat_fpkm,file = 'dat.fpkm.tcga.csv',row.names = T)
FPKM2TPM <- function(fpkm){
  exp(log(fpkm) - log(sum(fpkm)) + log(1e6))
}

dat_tpm <- apply(dat_fpkm,2,FPKM2TPM)
write.csv(dat_tpm,file = 'dat.tpm.tcga.csv',row.names = T)
