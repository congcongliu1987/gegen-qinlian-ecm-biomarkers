rm(list = ls()) 
setwd('/home/data/project/30.JNPJ-6001-1/')
file.dir = '16_scExpression'
if (! dir.exists(file.dir)){
  dir.create(file.dir)
}
setwd(file.dir)
load('../02_scRNA/CRC/03_celltype/UMAP.celltype.Rdata')
pacman::p_load(reshape2,tidyverse,Seurat,rstatix,ggpubr,gridExtra)
UMAP.som <- UMAP
table(UMAP.som$celltype)
som.scobj <- UMAP.som
table(som.scobj$celltype)
gene <- read.csv('../06_113model/bestModel.gene.csv')
p1 = FeaturePlot(som.scobj, features = gene$symbol,repel = T, label = , slot = "data", 
                 cols = c("lightgrey","blue"), ncol =3, order = T)
p1
ggsave(filename = '01.feature.pdf',p1,w=10,h=6)
ggsave(filename = '01.feature.png',p1,w=10,h=6,units = 'in',dpi = 600)
