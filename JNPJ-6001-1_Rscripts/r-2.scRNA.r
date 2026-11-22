rm(list = ls())
setwd('/home/data/project/30.JNPJ-6001-1/')
file.dir = '02_scRNA'
if (! dir.exists(file.dir)){
  dir.create(file.dir)
}
setwd(file.dir)
##两个数据
##UC GSE214695--------
rm(list = ls())
setwd('/home/data/project/30.JNPJ-6001-1/02_scRNA/')
file.dir = 'UC'
if (! dir.exists(file.dir)){
  dir.create(file.dir)
}
setwd(file.dir)
file.dir = '00_rawdata'
if (! dir.exists(file.dir)){
  dir.create(file.dir)
}
setwd(file.dir)
library(Seurat)
library(SeuratDisk)
library(pacman)
pacman::p_load(data.table, tidyverse,Seurat, patchwork)

## PART A导入数据 --------------------------------------------------------------

folders <- list.files('GSE214695_raw/')
folders
group <- data.frame(sample=folders,group=c(rep('Health',6),rep('UC',6)))
table(group$group)
dir <- paste0('GSE214695_raw/',folders,'/')
dir
scList <- list()
for (i in c(1:12)) {
  count =  CreateSeuratObject(counts = Read10X(dir[i]),
                              min.cells = 3, min.features = 200)
  scList[[i]] <- count
}

scRNA <- merge(scList[[1]],
               y = c(scList[2:12]), 
               add.cell.ids = folders)

table(scRNA$orig.ident)
head(scRNA@meta.data)
##添加分组
smps <- Cells(scRNA) %>% str_sub(1,10)
table(smps)
grps <- group
group.label <- grps$group[match(smps,grps$sample)]
names(group.label) <-  Cells(scRNA)
group.label
table(group.label)
# Health     UC 
# 20351  17082 
scRNA <- AddMetaData(scRNA,group.label,col.name = 'group')

scRNA[["percent.mt"]] <- PercentageFeatureSet(scRNA, pattern = "^MT-")  
###计算核糖体基因比例
scRNA[["percent.rb"]] <- PercentageFeatureSet(scRNA, pattern = "^RP[SL]")

head(scRNA@meta.data)
system.time(save(scRNA, file = "01.scRNA_orig.Rdata"))  

#PART B QC---------
rm(list = ls())
setwd('/home/data/project/30.JNPJ-6001-1/02_scRNA/UC/')
if (! dir.exists("./01_QC")){
  dir.create("./01_QC")
}
setwd("./01_QC")
load('../00_rawdata/01.scRNA_orig.Rdata')

plot.featrures = c("nFeature_RNA", "nCount_RNA", "percent.mt")   
# 质控前小提琴图
library(patchwork)
library(ggplot2)
library(RColorBrewer)
library(IOBR)

plots <- list()
for (i in seq_along(plot.featrures)) {
  plots[[i]] <- VlnPlot(scRNA, 
                        pt.size = 0,
                        features = plot.featrures[i], 
                        cols = palettes(category = "random",22,show_col=F)) +
    NoLegend() +
    ggtitle(plot.featrures[i]) +  # 增加每个图的标题
    theme_minimal() +             # 使用极简主题
    theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5), # 标题字体调整
      axis.title = element_text(size = 12),                              # 坐标轴标题字体大小
      axis.text = element_text(size = 10),                               # 坐标轴标签字体大小
      panel.grid.major = element_line(color = "gray", size = 0.5),       # 主格网线
      panel.grid.minor = element_blank(),                                 # 隐藏次格网线                              # 坐标轴标签字体大小
      axis.text.x = element_text(angle = 45, hjust = 1)
    ) + NoLegend()
}

# 使用 wrap_plots() 组合图形
violin <- wrap_plots(plots = plots, nrow = 1, ncol = 3) & 
  plot_annotation(
    #  title = "Violin Plots of Features",
    theme = theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5)
    )
  )

print(violin)

ggsave("01.vlnplot_before_qc.pdf", plot = violin, width = 8, height = 4)

ggsave("01.vlnplot_before_qc.png", plot = violin, width = 8, height = 4,dpi = 600)

### 质控后 ###

pctMT=15
minGene=200
maxGene=4000
maxUMI=10000

scRNA <- subset(scRNA, subset =
                  nCount_RNA <maxUMI &
                  nFeature_RNA > minGene & 
                  nFeature_RNA < maxGene& 
                  percent.mt < pctMT)

# 设置绘图元素
plot.featrures = c("nFeature_RNA", "nCount_RNA", "percent.mt")
#GSE = "Patients_ID"

# 质控后小提琴图
plots <- list()
for (i in seq_along(plot.featrures)) {
  plots[[i]] <- VlnPlot(scRNA, 
                        pt.size = 0,
                        features = plot.featrures[i], 
                        cols = palettes(category = "random",22,show_col=F)) +
    NoLegend() +
    ggtitle(plot.featrures[i]) +  # 增加每个图的标题
    theme_minimal() +             # 使用极简主题
    theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5), # 标题字体调整
      axis.title = element_text(size = 12),                              # 坐标轴标题字体大小
      axis.text = element_text(size = 10),                               # 坐标轴标签字体大小
      panel.grid.major = element_line(color = "gray", size = 0.5),       # 主格网线
      panel.grid.minor = element_blank(),                                 # 隐藏次格网线                              # 坐标轴标签字体大小
      axis.text.x = element_text(angle = 45, hjust = 1)
    ) + NoLegend()
}

# 使用 wrap_plots() 组合图形
violin <- wrap_plots(plots = plots, nrow = 1, ncol = 3) & 
  plot_annotation(
    #  title = "Violin Plots of Features",
    theme = theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5)
    )
  )

print(violin)
ggsave("02.vlnplot_after_qc.pdf", plot = violin, width = 8, height = 4)
ggsave("02.vlnplot_after_qc.png", plot = violin, width = 8, height = 4,dpi = 600)

system.time(save(scRNA, file = "02.scRNA_qc.Rdata"))  

# 03 PCA UMAP降维聚类-------
rm(list = ls())
setwd('/home/data/project/30.JNPJ-6001-1/02_scRNA/UC/')
if (! dir.exists("./02_PCA")){
  dir.create("./02_PCA")
}
setwd("./02_PCA")
load('../01_QC/02.scRNA_qc.Rdata')

scRNA <- JoinLayers(scRNA)     #将count进行合并
library(Seurat)
scRNA.norm<-NormalizeData(scRNA,normalization.method = "LogNormalize",scale.factor = 10000)
all.genes<-rownames(scRNA.norm)
scRNA.norm<-FindVariableFeatures(scRNA.norm,selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(scRNA.norm), 10)
p <- VariableFeaturePlot(scRNA.norm)
p <- LabelPoints(plot = p, 
                 points = top10, 
                 repel = T)
p
ggsave(filename = "01.feature_selection.pdf",p,w=7,h=5)
ggsave(filename = "01.feature_selection.png",p,w=7,h=5,units = 'in',dpi = 600,bg = 'white')

scRNA.nor.sca<-ScaleData(scRNA.norm)

## PCA降维
scRNA.norm.pca<-RunPCA(scRNA.nor.sca,features = VariableFeatures(object = scRNA.nor.sca),npcs = 50)
print(scRNA.norm.pca[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(scRNA.norm.pca, dims = 1:2, reduction = "pca")
scRNA.norm.pca <- JackStraw(scRNA.norm.pca, num.replicate = 100, dims = 50)
scRNA.norm.pca <- ScoreJackStraw(scRNA.norm.pca, dims = 1:50)
plot_pca <- JackStrawPlot(scRNA.norm.pca, dims = 1:50)
plot_pca
ggsave(filename = '02.pca_cluster.pdf',plot_pca,w=13,h=8)
ggsave(filename = '02.pca_cluster.png',plot_pca,w=13,h=8,units = 'in',dpi = 600,bg = 'white')
plot_elbow <- ElbowPlot(scRNA.norm.pca, ndims = 50)
plot_elbow
ggsave(filename = '03.pca_sd.pdf',plot_elbow,w=6,h=5)
ggsave(filename = '03.pca_sd.png',plot_elbow,w=6,h=5,units = 'in',dpi = 600,bg = 'white')

scRNA.norm.pca.c<-FindNeighbors(scRNA.norm.pca,dims = 1:20)
scRNA.norm.pca.c<-FindClusters(scRNA.norm.pca.c,resolution = 0.4)
UMAP<-RunUMAP(scRNA.norm.pca.c,dims = 1:20)

DimPlot(UMAP,reduction = 'umap')
pdf('04.UMAP.pdf',w=8,h=7)
DimPlot(UMAP,reduction = 'umap',cols = palettes(category = "random",22,show_col=F))
dev.off()
png('04.UMAP.png',w=8,h=7,units = 'in',res = 600)
DimPlot(UMAP,reduction = 'umap',cols = palettes(category = "random",22,show_col=F))
dev.off()

system.time(save(UMAP, file = "UMAP.Rdata"))

# 04.celltype -----------------------------------------------------
rm(list = ls())
setwd('/home/data/project/30.JNPJ-6001-1/02_scRNA/UC/')
if (! dir.exists("./03_celltype")){
  dir.create("./03_celltype")
}
setwd("./03_celltype")
load('../02_PCA/UMAP.Rdata')
length(colnames(UMAP))

all.markers<-FindAllMarkers(UMAP,only.pos = TRUE,min.pct = 0.1,logfc.threshold = 1,test.use = 'wilcox',return.thresh = 0.01)
human.marker <- readxl::read_xlsx('../../../../../files/cellmarker/Cell_marker_Human.xlsx')
table(human.marker$tissue_type)
table(human.marker$cell_name)
map.marker <- human.marker[which(human.marker$marker %in% all.markers$gene),]
colnames(map.marker)
map.marker <- map.marker[,c(6,7,9)]
colnames(map.marker)[3] <- 'gene'
map.marker <- merge(map.marker,all.markers,by = 'gene')

library(Seurat)
library(SingleR)
scRNA_singleR <- GetAssayData(UMAP, slot = "data")
library(celldex)
library(ggplot2)
hpca.se <- HumanPrimaryCellAtlasData()
clusters=UMAP@meta.data$seurat_clusters
pred.hesc <- SingleR(test = scRNA_singleR,
                     ref = hpca.se,
                     labels = hpca.se$label.main,
                     method = "cluster",
                     clusters = clusters,
                     quantile = 0.8 ,
                     assay.type.test = "logcounts",
                     assay.type.ref = "logcounts")
table(pred.hesc$labels)

celltype = data.frame(ClusterID=rownames(pred.hesc), 
                      celltype=pred.hesc$labels, 
                      stringsAsFactors = F) 



library(Seurat)
library(SingleR)

list_gene <- list(
  'Epithelium' = c('EPCAM','AQP8','BEST4','MUC2','OLFM4','PLCG2','TRPM5','ZG16'),
  'B and Plasma' = c('CD79A', 'BANK1', 'CD19', 'DERL3', 'MS4A1', 'MZB1'),
  'T cells' = c('CD3D', 'CD3E', 'CD3G', 'CD8A', 'FOXP3', 'GZMA', 'GZMB', 'IL17A', 'NKG7', 'TRBC1'),
  'Stromal'=c('ACTA2', 'ADAMDEC1', 'CHI3L1', 'COL3A1', 'NRXN1', 'PVALP', 'SOX6', 'VWF'),
  'Myeloid' = c('AIF1', 'C1QA', 'C1QB', 'CD14', 'CMTM2', 'FCGR3B', 'LYZ', 'MS4A2', 'TPSAB1','TPSAB2')
)

allgene <- data.frame(rownames(UMAP))

gene<-c('EPCAM','AQP8','BEST4','MUC2','OLFM4','PLCG2','TRPM5','ZG16','CD79A', 'BANK1', 'CD19', 'DERL3', 'MS4A1', 'MZB1',
        'CD3D', 'CD3E', 'CD3G', 'CD8A', 'FOXP3', 'GZMA', 'GZMB', 'IL17A', 'NKG7', 'TRBC1','ACTA2', 'ADAMDEC1', 'CHI3L1', 'COL3A1', 'NRXN1', 'PVALP', 'SOX6', 'VWF',
        'AIF1', 'C1QA', 'C1QB', 'CD14', 'CMTM2', 'FCGR3B', 'LYZ', 'MS4A2', 'TPSAB1','TPSAB2')

theme.set = theme(
  axis.title = element_text(size = 18, face = "bold", family = "Times"),
  axis.text.x = element_text(angle = 45, vjust = 1,size = 14, family = "Times"),
  axis.text.y = element_text(size = 14,  face = "bold", family = "Times"),
  legend.text = element_text(size = 18, face = "bold", family = "Times"),
  legend.title = element_blank(),
  text = element_text(family = "Times"))

DotPlot(UMAP,features = list_gene)+RotatedAxis()+theme.set+
  scale_color_gradientn(colours = c('#619ac3','#e9ccd3','#ea7293'))

cluster.ids = c(
  "0" = "T cells",  
  "1" = "T cells",  
  "2" = "B and Plasma",  
  "3" = "T cells", 
  "4" = "Epithelium",   
  "5" = "B and Plasma",
  "6" = "Stromal",
  "7" = "T cells",
  "8" = "Myeloid",
  "9" = "Stromal",
  "10" = "Epithelium",
  "11" = "B and Plasma",
  "12" = "Epithelium",
  "13" = "Myeloid",
  "14" = "Myeloid",
  "15" = "T cells",
  "16" = "Myeloid",
  "17" = "Epithelium",
  "18" = "Stromal",
  "19" = "Stromal"
)

UMAP <- RenameIdents(UMAP, cluster.ids)
UMAP$celltype <- Idents(UMAP)
table(UMAP@meta.data$celltype)
library(IOBR)
### 绘制注释umap图
theme.set2 = theme(
  axis.title = element_text(size = 16, face = "bold", family = "Times"),
  axis.text.x = element_text(size = 12, face = "bold", family = "Times"),
  axis.text.y = element_text(size = 12, face = "bold", family = "Times"),
  legend.text = element_text(size = 14, face = "bold", family = "Times"),
  legend.title = element_blank(),
  text = element_text(family = "Times"))

pdf(file="01.celltype.pdf",width=9,height=7,family='Times')
DimPlot(UMAP ,reduction = "umap",
        label.size = 4,label.color = 'black',label = T,
        cols = palettes(category = "random",20,show_col=F))+theme.set2
dev.off()

png(file="01.celltype.png",width=9,height=7,family='Times',units='in',res=600)
DimPlot(UMAP ,reduction = "umap",
        label.size = 4,label.color = 'black',label =T,
        cols = palettes(category = "random",20,show_col=F))+theme.set2
dev.off()

phe=UMAP@meta.data
table(phe$celltype)

cell_type_stat <- as.data.frame(sort(table(phe$celltype)))
write.table(cell_type_stat,
            file = "cell_type_stat.txt",
            quote = F)
system.time(save(UMAP, file = "UMAP.celltype.Rdata"))  


library(stringr)
# genes_to_check=str_to_upper(unique(gene))
genes_to_check = gene

th=theme(axis.text.x = element_text(vjust = 0.5, hjust=0.5))
p_all_markers <- DotPlot(UMAP, features = genes_to_check,
                         assay='RNA' ,group.by = 'celltype' )  + coord_flip()+th
p_all_markers
data<-p_all_markers$data
data<-p_all_markers$data

colnames(data)
colnames(data)<-c("AverageExpression_unscaled","Precent Expressed","Features","celltype","Average Expression")
unique(data$`Precent Expressed`)

####用ggplot画图####
p = ggplot(data,aes(celltype,Features,size = `Precent Expressed` ))+
  geom_point(shape=21,aes(fill= `Average Expression`),position =position_dodge(0))+
  theme_minimal()+xlab(NULL)+ylab(NULL) +
  scale_size_continuous(range=c(1,10))+theme_bw()+
  scale_fill_gradient(low = "grey", high = "#E54924")+
  theme(legend.position = "right",legend.box = "vertical", #图例位置
        legend.margin=margin(t= 0.5, unit='cm'),
        legend.spacing = unit(0,"in"),
        axis.text.x  = element_text(color="black",size=15,angle = 90,
                                    vjust = 0.25, hjust=0.25),#x轴
        axis.text.y  = element_text(color="black",size=12),#y轴
        legend.text = element_text(size =10,color="black"),#图例
        legend.title = element_text(size =10,color="black"),#图例
        axis.title.y=element_text(vjust=1,
                                  size=16)
  )+labs(x=" ",y = "Features")
p
p1 = ggplot(data,aes(celltype,Features,size = `Precent Expressed` ))+
  geom_point(shape=21,aes(fill= `Average Expression`),position =position_dodge(0))+
  theme_minimal()+xlab(NULL)+ylab(NULL) +
  scale_size_continuous(range=c(1,10))+theme_bw()+
  scale_fill_gradient(low = "#E54924", high = "#498EA4")+
  theme(legend.position = "right",legend.box = "vertical",
        legend.margin=margin(t= 0, unit='cm'),
        legend.spacing = unit(0,"in"),
        axis.text.x  = element_text(color="black",size=15,angle = 90,
                                    vjust = 0.5, hjust=0.5),
        axis.text.y  = element_text(color="black",size=12),
        legend.text = element_text(size =10,color="black"),
        legend.title = element_text(size =10,color="black"),
        axis.title.y=element_text(vjust=1,
                                  size=16)
  )+labs(x=" ",y = "Features")

# 保存
library(patchwork)
p+p1
ggsave(file="02.dotplot_celltype.pdf",p,height= 10,width = 6)
ggsave(file="02.dotplot_celltype.png",p,height= 10,width = 6,units = 'in',dpi = 600)

table(UMAP$celltype)


color2=palettes(category = "random",25,show_col=F)

df=UMAP@reductions$umap@cell.embeddings %>%
  as.data.frame() %>%
  cbind(celltype=UMAP@meta.data$celltype)

head(df)
library(ggunchull)
library(Cairo)
p_tsne1=ggplot(df, aes(umap_1, umap_2, color=celltype))+
  geom_point(size = 0.5) + #细胞亚群点的大小
  scale_color_manual(values=color2)+
  theme(panel.border = element_blank(), #边框
        axis.title = element_blank(),  #轴标题
        axis.text = element_blank(), # 文本
        axis.ticks = element_blank(),
        panel.background = element_rect(fill = 'white'), #背景色
        plot.background=element_rect(fill="white"),
        legend.title = element_blank(), #去掉legend.title 
        legend.key=element_rect(fill='white'), 
        legend.text = element_text(size=14), #设置legend标签的大小
        legend.key.size=unit(0.6,'cm') )+
  stat_ellipse(aes(x = umap_1, y = umap_2, fill = celltype),
               geom = "polygon",#多边形
               linetype=2, # 椭圆，虚线
               alpha = 0.05,  #椭圆不透明度
               show.legend = FALSE, #去掉椭圆对应的图例
               level = 0.93)+ #置信区间
  guides(fill= guide_legend(override.aes = list(size = 2)))+ #图例点的大小
  scale_color_manual(values=color2)

celltype_position <- df %>%
  group_by(celltype) %>%
  dplyr::summarise(
    umap_1 = median(umap_1),
    umap_2 = median(umap_2));p_tsne1

library(ggrepel)
#install.packages('tidydr')
library(tidydr)
p_tsne2 <- p_tsne1+theme_dr(xlength=0.22,ylength=0.22,
                            arrow=grid::arrow(length=unit(0.15,'inches'),type='closed'))+
  theme(panel.grid = element_blank());p_tsne2



pdf(file="03.celltype_UMAP.pdf",width=6,height=6,family='Times')
p_tsne3<-p_tsne2+
  geom_label_repel(data = celltype_position,aes(label=celltype), 
                   fontface="bold",point.padding=unit(0.1, "lines"))+theme(legend.position = "none");p_tsne3

dev.off()


png(file="03.celltype_UMAP.png",width=6,height=6,family='Times',units='in',res=600)
p_tsne3<-p_tsne2+
  geom_label_repel(data = celltype_position,aes(label=celltype), 
                   fontface="bold",point.padding=unit(0.1, "lines"))+theme(legend.position = "none");p_tsne3

dev.off()

# 05.AUCell -----------------------------------------------------
rm(list = ls())
setwd('/home/data/project/30.JNPJ-6001-1/02_scRNA/UC/')
if (! dir.exists("./04_AUCell")){
  dir.create("./04_AUCell")
}
setwd("./04_AUCell")
load('../03_celltype/UMAP.celltype.Rdata')

library(AUCell)
library(ggplot2)
library(patchwork)
library(dplyr)
library(Seurat)

cells_rankings <- AUCell_buildRankings(UMAP@assays$RNA$data) 
geneSets <- read.csv('../../../01_Target/ecm.geneset.csv')
geneSets <- geneSets$symbol
names(geneSets) <- 'ECM.score'
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings, aucMaxRank=nrow(cells_rankings)*0.1)

aucs <- as.numeric(getAUC(cells_AUC))
UMAP$AUCell  <- aucs
library(ggraph)
p <- ggplot(data.frame(UMAP@meta.data, UMAP@reductions$umap@cell.embeddings), aes(umap_1, umap_2, color=AUCell)
) + geom_point( size=1.5
) + scale_color_viridis(option="A")  + theme_light(base_size = 15)+labs(title = "AUCell score")+
  theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"))+
  theme(plot.title = element_text(hjust = 0.5))
p
ggsave(filename = '01.AUcell.pdf',p,w=6,h=5)
ggsave(filename = '01.AUcell.png',p,w=6,h=5,units = 'in',dpi = 600,bg = 'white')

##高低评分组
aucs.group <- data.frame(cell = colnames(UMAP),auc.score = aucs)
aucs.group$group <- ifelse(aucs.group$auc.score > median(aucs.group$auc.score),'High.score','Low.score')
UMAP <- AddMetaData(UMAP,aucs.group$group,col.name = 'Score.group')
table(UMAP$Score.group)

p1 <- DimPlot(UMAP,reduction = 'umap',cols = c('red','blue'),group.by = 'Score.group')
p1
ggsave(filename = '02.AUcell.group.pdf',p1,w=6,h=5)
ggsave(filename = '02.AUcell.group.png',p1,w=6,h=5,units = 'in',dpi = 600,bg = 'white')

all.count <- data.frame(celltype=UMAP$celltype, sample=UMAP$Score.group)
library(reshape2)
plot.celltype = table(all.count$celltype, all.count$sample) %>% as.data.frame() %>% recast(Var1 ~ Var2)
colnames(plot.celltype)[1] = "Cell"
plot.celltype <- plot.celltype[,c(1,3,2)]
table(UMAP$Score.group)
plot.celltype = cbind(plot.celltype, apply(plot.celltype[-1], 2, proportions))
colnames(plot.celltype)[2:3] <- paste0('Count_',colnames(plot.celltype)[2:3])
colnames(plot.celltype)[4:5] = paste0('Proportion_',colnames(plot.celltype)[4:5])
plot.celltype <- plot.celltype[order(plot.celltype$Proportion_High.score,decreasing = T),]
write.table(plot.celltype, "CellCount.xls", row.names = F, col.names = T, quote = F, sep = "\t")

dat.celltype = plot.celltype[c(1,4:5)]
colnames(dat.celltype)[2:3] = gsub('Proprotion','',colnames(dat.celltype)[2:3])
dat.celltype[2:3] = 100 * dat.celltype[2:3]

dat.plot = melt(dat.celltype, id.vars = "Cell", variable.name = "Sample", value.name = "value")
dat.plot$Cell = factor(dat.plot$Cell, levels = plot.celltype$Cell)
library(IOBR)
library(ggplot2)
library(IOBR)
dat.plot <- dat.plot[order(dat.plot$value),]
dat.plot$Cell = factor(dat.plot$Cell, levels = plot.celltype$Cell)
library(IOBR)
dat.plot$label = round(dat.plot$value, 2) %>% paste0(.,"%")
str(dat.plot)
dat.plot$Sample <- factor(dat.plot$Sample,levels = c('Proportion_High.score','Proportion_Low.score'))
ggplot(data = dat.plot, mapping = aes(x = value, y = Cell)) + 
  geom_text(aes(x = value + 1, label = label), hjust = 0) +
  geom_bar(aes(fill = Cell), stat = "identity") + 
  facet_wrap(~Sample) +
  scale_fill_manual(values = palettes(category = "random",22,show_col=F)) + 
  xlim(c(0,70)) +
  xlab("Cell Fraction") + ylab("Cell Type") +
  theme_minimal(base_size = 16) + 
  theme(panel.border = element_rect(color = "grey60", fill = "transparent"), legend.position = "none", 
        panel.grid = element_blank())
system.time(save(UMAP, file = "UMAP.celltype.ecm.Rdata"))  

# 06.Function enrichment -----------------------------------------------------
rm(list = ls())
setwd('/home/data/project/30.JNPJ-6001-1/02_scRNA/UC/')
if (! dir.exists("./05_Function")){
  dir.create("./05_Function")
}
setwd("./05_Function")
load('../04_AUCell/UMAP.celltype.ecm.Rdata')
table(UMAP$Score.group)
##GSVA------
group <- as.data.frame(UMAP$Score.group) %>% rownames_to_column(var = 'sample')
colnames(group)[2] <- 'group'
control.sample <- group$sample[which(group$group == 'Low.score')]
group2 <- group
colnames(group2)
colnames(group2)<-c('id','label')
gsva_exp <- as.data.frame(UMAP@assays$RNA$data) 
all(colnames(gsva_exp) == group2$id)
dim(gsva_exp)
# 分组
library(limma)
library(GSEABase)
library(GSVA)
group<- group2$label %>% as.factor()
design <- model.matrix(~0 + group)
rownames(design) <- colnames(gsva_exp)
colnames(design) <- levels(group)
compare<- makeContrasts("High.score-Low.score", levels = design)
KEGG_ref <- getGmt("../../../../../files/GSVA/c2.cp.kegg.v7.5.1.symbols.gmt")
gsvapar <- gsvaParam(
  as.matrix(gsva_exp),
  KEGG_ref,
  assay = NA_character_,
  annotation = NA_character_,
  minSize = 10,
  maxSize = Inf,
  kcdf = c("Gaussian", "Poisson", "none"),
  tau = 1,
  maxDiff = TRUE,
  absRanking = FALSE
)

es_KEGG <- gsva(gsvapar)
es_KEGG <- as.data.frame(es_KEGG)
saveRDS(es_KEGG,file = 'es_KEGG.rds')

fit <- lmFit(es_KEGG, design)
fit2 <- contrasts.fit(fit ,compare)
fit3 <- eBayes(fit2)
allGeneSets <- topTable(fit3, coef = 1, number = Inf)
logFCcutoff <- 0
allGeneSets$change = as.factor(
  ifelse(allGeneSets$adj.P.Val < 0.05 & abs(allGeneSets$logFC) > logFCcutoff,
         ifelse(allGeneSets$logFC > logFCcutoff,'UP','DOWN'),'NOT')
)
table(allGeneSets$change)
DEGeneSets <- subset(allGeneSets,
                     allGeneSets$adj.P.Val < 0.05 & abs(allGeneSets$logFC) > logFCcutoff )

write.table(allGeneSets, 
            file = "KEGG.xls",
            quote = F,
            sep = "\t",
            row.names = T)
DEGeneSets <- DEGeneSets[order(DEGeneSets$adj.P.Val),]
dim(DEGeneSets)
write.table(DEGeneSets,
            file = 'diff_KEGG.xls',
            quote = F,
            sep = "\t",
            row.names = T)
library(ggprism)
## barplot
dat_plot<-data.frame(id=rownames(DEGeneSets),
                     t=DEGeneSets$t)
allGeneSets <- allGeneSets[dat_plot$id,]
dat_plot$threshold = factor(allGeneSets$change)
# dat_plot$threshold = factor(ifelse(dat_plot$t  >-1, ifelse(dat_plot$t >= 1 ,'Up','NoSignifi'),'Down'),levels=c('Up','Down','NoSignifi'))
table(dat_plot$threshold)
dat_plot$threshold <- ifelse(dat_plot$threshold=='UP','Up',
                             ifelse(dat_plot$threshold=='DOWN','Down','NoSignifi'))
dat_plot<-dat_plot%>%arrange(t)
dat_plot <- dat_plot[which(dat_plot$id %in% rownames(DEGeneSets)),]
dat_plot$id<-factor(dat_plot$id,levels=dat_plot$id)

p <- ggplot(data = dat_plot,aes(x = id,y = t,fill = threshold)) +
  geom_col()+
  coord_flip() +
  scale_fill_manual(values = c('Up'= '#FF0088','NoSignifi'='#cccccc','Down'='#7bcd7b')) +
  #  geom_hline(yintercept = c(-1,1),color = 'white',size = 0.5,lty='dashed') +
  scale_y_continuous(breaks = c(seq(-20,20,5))) +
  xlab('') + 
  ylab('t value of GSVA score, High versus Low') + #注意坐标轴旋转了
  guides(fill=F)+ # 不显示图例
  theme_prism(border = T) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )
p
low1 <- dat_plot %>% filter(t < -1) %>% nrow()
low0 <- dat_plot %>% filter( t < 0) %>% nrow()
high0 <- dat_plot %>% filter(t < 1) %>% nrow()
high1 <- nrow(dat_plot)

# 依次从下到上添加标签
p <- p + geom_text(data = dat_plot[1:low1,],aes(x = id,y = 0.1,label = id),
                   hjust = 0,color = 'black')  + # 小于-1的为黑色标签
  geom_text(data = dat_plot[(high0 +1):high1,],aes(x = id,y = -0.1,label = id),
            hjust = 1,color = 'black') # 大于1的为黑色标签
p
ggsave(filename = '01.GSVA.pdf',p,w=13,h=25)
ggsave(filename = '01.GSVA.png',p,w=13,h=25,units = 'in',dpi = 600)


##GO--------
library(tidyverse)
library(Seurat)
UMAP.som <- UMAP
table(UMAP.som$Score.group)
all.markers<-FindMarkers(UMAP.som,only.pos = TRUE,min.pct = 0.1,logfc.threshold = 1,test.use = 'wilcox',
                         return.thresh = 1,ident.1 = 'High.score',ident.2 = 'Low.score',group.by = 'Score.group')
select.markers <- all.markers[which(all.markers$p_val_adj < 0.05 & abs(all.markers$avg_log2FC) > 1),]
select.markers$gene <- rownames(select.markers)
library(Seurat)
library(gplots)
library(ggplot2)
library(clusterProfiler)
library(org.Hs.eg.db)
ids=bitr(select.markers$gene,'SYMBOL','ENTREZID','org.Hs.eg.db') ## 将SYMBOL转成ENTREZID
ego <- enrichGO(gene = ids$ENTREZID,
                OrgDb = org.Hs.eg.db,
                keyType = "ENTREZID",
                ont = "ALL",
                pAdjustMethod = "BH",
                pvalueCutoff = 0.05,
                qvalueCutoff = 1,
                readable = TRUE)
write.csv(ego,file = "GO.csv",row.names = F)
go_bar<-barplot(ego, showCategory=5, split="ONTOLOGY",label_format = 70) +
  facet_grid(ONTOLOGY ~ ., scales = "free",space = 'free')
go_bar
ggsave('02.GO_bar.png',go_bar,width =8,height = 6,units = 'in',dpi = 600)
ggsave('02.GO_bar.pdf',go_bar,width =8,height = 6)

library(DOSE)
erich.do<-DOSE::enrichDO(gene=ids$ENTREZID,
                         ont = "HDO",organism = 'hsa',
                         pvalueCutoff = 0.05,
                         qvalueCutoff = 0.05,
                         readable = T)
res.do <- setReadable(erich.do,OrgDb = 'org.Hs.eg.db',keyType = 'SYMBOL')
write.csv(ego,file = "DO.csv",row.names = F)
# 展示富集最显著的 GO term
do_dot<-dotplot(erich.do, showCategory=15,label_format = 70)
do_dot
ggsave('03.DO_dot.png',do_dot,width =6,height = 6,units = 'in',dpi = 600)
ggsave('03.DO_dot.pdf',do_dot,width =6,height = 6)


# 07.DEcell -----------------------------------------------------
rm(list = ls())
setwd('/home/data/project/30.JNPJ-6001-1/02_scRNA/UC/')
if (! dir.exists("./06_DEcell")){
  dir.create("./06_DEcell")
}
setwd("./06_DEcell")

pacman::p_load(tidyverse, reshape2, ggpubr, ggplot2, ggsci)
load('../03_celltype/UMAP.celltype.Rdata')
colnames(UMAP)
all.count <- data.frame(celltype = UMAP$celltype,
                        group = UMAP$group,
                        sample = stringr::str_extract(colnames(UMAP),"^GSM\\d+"))

table(all.count$celltype)
table(all.count$group)
table(all.count$sample)
plot.celltype <- table(all.count$celltype, all.count$sample ) %>% as.data.frame() %>% recast(Var1 ~ Var2)
colnames(plot.celltype)[1] = "Cell"
plot.celltype <- cbind(plot.celltype, apply(plot.celltype[-1], 2, proportions))
colnames(plot.celltype)
percent.celltype <- plot.celltype[,c(1,14:25)]
percent.celltype[2:13] = 100 *percent.celltype[2:13] 

dat.celltype <- melt(percent.celltype)
colnames(dat.celltype) <- c('Cell','Sample','Percentage')

##添加一下分组
dat.celltype$Group <- all.count$group[match(dat.celltype$Sample,all.count$sample)]

##wilcox检验
vln.df <- dat.celltype
table(vln.df$Group)
vln.df$Group <- factor(vln.df$Group,levels = c('Health','UC'))
vln.df$Percentage
exp_plot <- ggplot(vln.df,aes(x = Group, y = Percentage, fill = Group)) +
  geom_violin(trim=F,color="black",aes(fill=Group)) + #绘制小提琴图, “color=”设置小提琴图的轮廓线的颜色(#不要轮廓可以设为white以下设为背景为白色，其实表示不要轮廓线)
  stat_boxplot(geom="errorbar", 
               width=0.1,
               position = position_dodge(0.9)) +
  geom_boxplot(width=0.4,
               position=position_dodge(0.9),
               outlier.shape = NA,fill = 'white')+ #绘制箱线图，此处width=0.1控制小提琴图中箱线图的宽窄
  scale_fill_manual(values= c("#5cb3cc","#ed556a"), name = "Group")+
  scale_colour_manual(values = c("#5cb3cc","#ed556a"))+
  labs(title="", x="", y = "Percentage",size=20) +
  scale_color_brewer(palette = 'Set2')+
  stat_compare_means(data = vln.df,
                     mapping = aes(group = Group),
                     label ="p.signif",
                     method = 'wilcox.test',
                     paired = F,
                     label.x = 1.4
  ) +
  theme_bw()+
  theme(plot.title = element_text(hjust =0.5,colour="black",face="bold",size=15),
        axis.text.x=element_text(angle=45,hjust=1,colour="black",face="bold",size=12), 
        axis.text.y=element_text(hjust=0.5,colour="black",face="bold",size=12), 
        axis.title.x=element_text(size=16,face="bold"),
        axis.title.y=element_text(size=16,face="bold"),
        legend.text=element_text(face = "bold", hjust = 0.5,colour="black", size=12),
        legend.title = element_text(face = "bold", size = 12),
        legend.position = "top",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  guides(fill='none')+facet_wrap(~Cell,ncol = 5,scales = 'free')
exp_plot

ggsave(filename = '01.DEcell.pdf',exp_plot,w=10,h=4)
ggsave(filename = '01.DEcell.png',exp_plot,w=10,h=4,dpi = 600)

# 循环生成每个细胞类型的图
cell_types <- unique(dat.celltype$Cell)
for (cell in cell_types) {
  vln.df <- dat.celltype[dat.celltype$Cell == cell,]
  
  exp_plot <- ggplot(vln.df, aes(x = Group, y = Percentage, fill = Group)) +
    geom_violin(trim = FALSE, color = "black", aes(fill = Group)) +
    stat_boxplot(geom = "errorbar", width = 0.1, position = position_dodge(0.9)) +
    geom_boxplot(width = 0.4, position = position_dodge(0.9), outlier.shape = NA, fill = 'white') +
    scale_fill_manual(values = c("#5cb3cc", "#ed556a"), name = "Group") +
    scale_colour_manual(values = c("#5cb3cc", "#ed556a")) +
    labs(title = "", x = "", y = "Percentage", size = 20) +
    scale_color_brewer(palette = 'Set2') +
    stat_compare_means(data = vln.df, 
                       mapping = aes(group = Group), 
                       label = "p.signif", 
                       method = 'wilcox.test', 
                       paired = FALSE, 
                       label.x = 1.4) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, colour = "black", face = "bold", size = 15),
          axis.text.x = element_text(angle = 45, hjust = 1, colour = "black", face = "bold", size = 12),
          axis.text.y = element_text(hjust = 0.5, colour = "black", face = "bold", size = 12),
          axis.title.x = element_text(size = 16, face = "bold"),
          axis.title.y = element_text(size = 16, face = "bold"),
          legend.text = element_text(face = "bold", hjust = 0.5, colour = "black", size = 12),
          legend.title = element_text(face = "bold", size = 12),
          legend.position = "top", 
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank()) +
    guides(fill = 'none') +
    facet_wrap(~Cell, ncol = 5, scales = 'free')
  
  # 保存每个细胞类型的图
  ggsave(filename = paste0(cell, "_DEcell.pdf"), exp_plot, width = 3, height = 4)
  ggsave(filename = paste0(cell, "_DEcell.png"), exp_plot, width = 3, height = 4, dpi = 600)
}



## 差异cell
table(UMAP$celltype)
UMAP.som <- subset(UMAP, subset = celltype %in% c('Myeloid'))

system.time(save(UMAP.som, file = "UMAP.main.celltype.Rdata"))  


##CRC GSE132257--------
rm(list = ls())
setwd('/home/data/project/30.JNPJ-6001-1/')
file.dir = '02_scRNA'
if (! dir.exists(file.dir)){
  dir.create(file.dir)
}
setwd(file.dir)
rm(list = ls())
setwd('/home/data/project/30.JNPJ-6001-1/02_scRNA/')
file.dir = 'CRC'
if (! dir.exists(file.dir)){
  dir.create(file.dir)
}
setwd(file.dir)
file.dir = '00_rawdata'
if (! dir.exists(file.dir)){
  dir.create(file.dir)
}
setwd(file.dir)
library(Seurat)
library(SeuratDisk)
library(pacman)
pacman::p_load(data.table, tidyverse,Seurat, patchwork)

## PART A导入数据 --------------------------------------------------------------

dat <- fread('GSE132257_GEO_processed_protocol_and_fresh_frozen_raw_UMI_count_matrix.txt.gz',check.names = F) %>% column_to_rownames(var = 'Index')
colnames(dat) <- gsub('.','-',colnames(dat),fixed = T)
anno <- read.delim2('GSE132257_processed_protocol_and_fresh_frozen_cell_annotation.txt')
table(anno$Cell_type)
group <- data.frame(sample = anno$Index,group = anno$Class)
table(group$group)
scRNA <- CreateSeuratObject(counts = dat,min.cells = 3,min.features = 200)

table(scRNA$orig.ident)
head(scRNA@meta.data)
length(colnames(scRNA))
# 18409 细胞数量

length(rownames(scRNA)) 
# 21977 基因数量

##添加分组
smps <- Cells(scRNA)
table(smps)
grps <- group
group.label <- grps$group[match(smps,grps$sample)]
names(group.label) <-  Cells(scRNA)
group.label
table(group.label)
# Normal  Tumor 
# 7686  10723 
scRNA <- AddMetaData(scRNA,group.label,col.name = 'group')

scRNA[["percent.mt"]] <- PercentageFeatureSet(scRNA, pattern = "^MT-")  
###计算核糖体基因比例
scRNA[["percent.rb"]] <- PercentageFeatureSet(scRNA, pattern = "^RP[SL]")

head(scRNA@meta.data)
system.time(save(scRNA, file = "01.scRNA_orig.Rdata"))  

#PART B QC---------
rm(list = ls())
setwd('/home/data/project/30.JNPJ-6001-1/02_scRNA/CRC/')
if (! dir.exists("./01_QC")){
  dir.create("./01_QC")
}
setwd("./01_QC")
load('../00_rawdata/01.scRNA_orig.Rdata')
# load('../00_rawdata/GSE214295_RAW/01.scRNA_orig.Rdata')
### 质控前 ###

# 设置绘图元素
plot.featrures = c("nFeature_RNA", "nCount_RNA", "percent.mt")   
#GSE = "Patients_ID"

# 质控前小提琴图
library(patchwork)
library(ggplot2)
library(RColorBrewer)
library(IOBR)

plots <- list()
for (i in seq_along(plot.featrures)) {
  plots[[i]] <- VlnPlot(scRNA, 
                        pt.size = 0,
                        features = plot.featrures[i], 
                        cols = palettes(category = "random",22,show_col=F)) +
    NoLegend() +
    ggtitle(plot.featrures[i]) +  # 增加每个图的标题
    theme_minimal() +             # 使用极简主题
    theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5), # 标题字体调整
      axis.title = element_text(size = 12),                              # 坐标轴标题字体大小
      axis.text = element_text(size = 10),                               # 坐标轴标签字体大小
      panel.grid.major = element_line(color = "gray", size = 0.5),       # 主格网线
      panel.grid.minor = element_blank(),                                 # 隐藏次格网线                              # 坐标轴标签字体大小
      axis.text.x = element_text(angle = 45, hjust = 1)
    ) + NoLegend()
}

# 使用 wrap_plots() 组合图形
violin <- wrap_plots(plots = plots, nrow = 1, ncol = 3) & 
  plot_annotation(
    #  title = "Violin Plots of Features",
    theme = theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5)
    )
  )

print(violin)
ggsave("01.vlnplot_before_qc.pdf", plot = violin, width = 12, height = 4)
ggsave("01.vlnplot_before_qc.png", plot = violin, width = 12, height = 4,dpi = 600)

### 质控后 ###

pctMT=15
minGene=200
maxGene=4000
maxUMI=10000

scRNA <- subset(scRNA, subset =
                  nCount_RNA <maxUMI &
                  nFeature_RNA > minGene & 
                  nFeature_RNA < maxGene& 
                  percent.mt < pctMT)
plot.featrures = c("nFeature_RNA", "nCount_RNA", "percent.mt")
plots <- list()
for (i in seq_along(plot.featrures)) {
  plots[[i]] <- VlnPlot(scRNA, 
                        pt.size = 0,
                        features = plot.featrures[i], 
                        cols = palettes(category = "random",22,show_col=F)) +
    NoLegend() +
    ggtitle(plot.featrures[i]) +  # 增加每个图的标题
    theme_minimal() +             # 使用极简主题
    theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5), # 标题字体调整
      axis.title = element_text(size = 12),                              # 坐标轴标题字体大小
      axis.text = element_text(size = 10),                               # 坐标轴标签字体大小
      panel.grid.major = element_line(color = "gray", size = 0.5),       # 主格网线
      panel.grid.minor = element_blank(),                                 # 隐藏次格网线                              # 坐标轴标签字体大小
      axis.text.x = element_text(angle = 45, hjust = 1)
    ) + NoLegend()
}

# 使用 wrap_plots() 组合图形
violin <- wrap_plots(plots = plots, nrow = 1, ncol = 3) & 
  plot_annotation(
    #  title = "Violin Plots of Features",
    theme = theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5)
    )
  )

print(violin)
ggsave("02.vlnplot_after_qc.pdf", plot = violin, width = 12, height = 4)
ggsave("02.vlnplot_after_qc.png", plot = violin, width = 12, height = 4,dpi = 600)

system.time(save(scRNA, file = "02.scRNA_qc.Rdata"))  

# 03 PCA UMAP降维聚类-------
rm(list = ls())
setwd('/home/data/project/30.JNPJ-6001-1/02_scRNA/CRC/')
if (! dir.exists("./02_PCA")){
  dir.create("./02_PCA")
}
setwd("./02_PCA")
load('../01_QC/02.scRNA_qc.Rdata')

scRNA <- JoinLayers(scRNA)     #将count进行合并
library(Seurat)
scRNA.norm<-NormalizeData(scRNA,normalization.method = "LogNormalize",scale.factor = 10000)
all.genes<-rownames(scRNA.norm)
scRNA.norm<-FindVariableFeatures(scRNA.norm,selection.method = "vst", nfeatures = 2000)
# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(scRNA.norm), 10)
# plot variable features with and without labels
p <- VariableFeaturePlot(scRNA.norm)
p <- LabelPoints(plot = p, 
                 points = top10, 
                 repel = T)
p
ggsave(filename = "01.feature_selection.pdf",p,w=7,h=5)
ggsave(filename = "01.feature_selection.png",p,w=7,h=5,units = 'in',dpi = 600,bg = 'white')

scRNA.nor.sca<-ScaleData(scRNA.norm)

## PCA降维
scRNA.norm.pca<-RunPCA(scRNA.nor.sca,features = VariableFeatures(object = scRNA.nor.sca),npcs = 50)
print(scRNA.norm.pca[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(scRNA.norm.pca, dims = 1:2, reduction = "pca")
scRNA.norm.pca <- JackStraw(scRNA.norm.pca, num.replicate = 100, dims = 50)
scRNA.norm.pca <- ScoreJackStraw(scRNA.norm.pca, dims = 1:50)
plot_pca <- JackStrawPlot(scRNA.norm.pca, dims = 1:50)
plot_pca
ggsave(filename = '02.pca_cluster.pdf',plot_pca,w=13,h=8)
ggsave(filename = '02.pca_cluster.png',plot_pca,w=13,h=8,units = 'in',dpi = 600,bg = 'white')
plot_elbow <- ElbowPlot(scRNA.norm.pca, ndims = 50)
plot_elbow
ggsave(filename = '03.pca_sd.pdf',plot_elbow,w=6,h=5)
ggsave(filename = '03.pca_sd.png',plot_elbow,w=6,h=5,units = 'in',dpi = 600,bg = 'white')

scRNA.norm.pca.c<-FindNeighbors(scRNA.norm.pca,dims = 1:20)
scRNA.norm.pca.c<-FindClusters(scRNA.norm.pca.c,resolution = 0.4)
UMAP<-RunUMAP(scRNA.norm.pca.c,dims = 1:20)

DimPlot(UMAP,reduction = 'umap')
pdf('04.UMAP.pdf',w=8,h=7)
DimPlot(UMAP,reduction = 'umap',cols = palettes(category = "random",22,show_col=F))
dev.off()
png('04.UMAP.png',w=8,h=7,units = 'in',res = 600)
DimPlot(UMAP,reduction = 'umap',cols = palettes(category = "random",22,show_col=F))
dev.off()

system.time(save(UMAP, file = "UMAP.Rdata"))

# 04.celltype -----------------------------------------------------
rm(list = ls())
setwd('/home/data/project/30.JNPJ-6001-1/02_scRNA/CRC/')
if (! dir.exists("./03_celltype")){
  dir.create("./03_celltype")
}
setwd("./03_celltype")
load('../02_PCA/UMAP.Rdata')
length(colnames(UMAP))

all.markers<-FindAllMarkers(UMAP,only.pos = TRUE,min.pct = 0.1,logfc.threshold = 1,test.use = 'wilcox',return.thresh = 0.01)

library(Seurat)
library(SingleR)
scRNA_singleR <- GetAssayData(UMAP, slot = "data")
library(celldex)
library(ggplot2)
hpca.se <- HumanPrimaryCellAtlasData()
#hpca.se<-BlueprintEncodeData()
#hpca.se<-DatabaseImmuneCellExpressionData()
clusters=UMAP@meta.data$seurat_clusters
pred.hesc <- SingleR(test = scRNA_singleR,
                     ref = hpca.se,
                     labels = hpca.se$label.main,
                     method = "cluster",
                     clusters = clusters,
                     quantile = 0.8 ,
                     assay.type.test = "logcounts",
                     assay.type.ref = "logcounts")
table(pred.hesc$labels)

celltype = data.frame(ClusterID=rownames(pred.hesc), 
                      celltype=pred.hesc$labels, 
                      stringsAsFactors = F) 

anno <- read.delim2('../00_rawdata/GSE132257_processed_protocol_and_fresh_frozen_cell_annotation.txt')
table(anno$Cell_type)

library(Seurat)

list_gene <- list(
  'Epithelium' = c('EPCAM','AQP8','BEST4','MUC2','OLFM4','PLCG2','TRPM5','ZG16'),
  'B cells' = c('CD79A', 'BANK1', 'CD19', 'DERL3', 'MS4A1', 'MZB1'),
  'T cells' = c('CD3D', 'CD3E', 'CD3G', 'CD8A', 'FOXP3', 'GZMA', 'GZMB', 'IL17A', 'NKG7', 'TRBC1'),
  'Stromal'=c('ACTA2', 'ADAMDEC1', 'CHI3L1', 'COL3A1', 'NRXN1', 'PVALP', 'SOX6', 'VWF'),
  'Myeloid' = c('AIF1', 'C1QA', 'C1QB', 'CD14', 'CMTM2', 'FCGR3B', 'LYZ', 'MS4A2', 'TPSAB1','TPSAB2'),
  'Mast cells' =c('TPSB2','KIT','CPA3','KITLG')
)

allgene <- data.frame(rownames(UMAP))

gene<-c('EPCAM','AQP8','BEST4','MUC2','OLFM4','PLCG2','TRPM5','ZG16','CD79A', 'BANK1', 'CD19', 'DERL3', 'MS4A1', 'MZB1',
        'CD3D', 'CD3E', 'CD3G', 'CD8A', 'FOXP3', 'GZMA', 'GZMB', 'IL17A', 'NKG7', 'TRBC1','ACTA2', 'ADAMDEC1', 'CHI3L1', 'COL3A1', 'NRXN1', 'PVALP', 'SOX6', 'VWF',
        'AIF1', 'C1QA', 'C1QB', 'CD14', 'CMTM2', 'FCGR3B', 'LYZ', 'MS4A2', 'TPSAB1','TPSAB2','TPSB2','KIT','CPA3','KITLG')

theme.set = theme(
  axis.title = element_text(size = 18, face = "bold", family = "Times"),
  axis.text.x = element_text(angle = 45, vjust = 1,size = 14, family = "Times"),
  axis.text.y = element_text(size = 14,  face = "bold", family = "Times"),
  legend.text = element_text(size = 18, face = "bold", family = "Times"),
  legend.title = element_blank(),
  text = element_text(family = "Times"))

# FeaturePlot(UMAP,features = gene,label = T,ncol = 2)

DotPlot(UMAP,features = list_gene)+RotatedAxis()+theme.set+
  scale_color_gradientn(colours = c('#619ac3','#e9ccd3','#ea7293'))


cluster.ids = c(
  "0" = "T cells",  
  "1" = "T cells",  
  "2" = "T cells",  
  "3" = "B cells", 
  "4" = "Epithelium",   
  "5" = "B cells",
  "6" = "Myeloid",
  "7" = "Epithelium",
  "8" = "Epithelium",
  "9" = "Stromal",
  "10" = "T cells",
  "11" = "Mast cells",
  "12" = "Stromal",
  "13" = "Stromal",
  "14" = "Stromal",
  "15" = "Stromal"
)

UMAP <- RenameIdents(UMAP, cluster.ids)
UMAP$celltype <- Idents(UMAP)
table(UMAP@meta.data$celltype)
library(IOBR)
### 绘制注释umap图
theme.set2 = theme(
  axis.title = element_text(size = 16, face = "bold", family = "Times"),
  axis.text.x = element_text(size = 12, face = "bold", family = "Times"),
  axis.text.y = element_text(size = 12, face = "bold", family = "Times"),
  legend.text = element_text(size = 14, face = "bold", family = "Times"),
  legend.title = element_blank(),
  text = element_text(family = "Times"))

pdf(file="01.celltype.pdf",width=9,height=7,family='Times')
DimPlot(UMAP ,reduction = "umap",
        label.size = 4,label.color = 'black',label = T,
        cols = palettes(category = "random",20,show_col=F))+theme.set2
dev.off()

png(file="01.celltype.png",width=9,height=7,family='Times',units='in',res=600)
DimPlot(UMAP ,reduction = "umap",
        label.size = 4,label.color = 'black',label =T,
        cols = palettes(category = "random",20,show_col=F))+theme.set2
dev.off()

phe=UMAP@meta.data
table(phe$celltype)

cell_type_stat <- as.data.frame(sort(table(phe$celltype)))
write.table(cell_type_stat,
            file = "cell_type_stat.txt",
            quote = F)
system.time(save(UMAP, file = "UMAP.celltype.Rdata"))  


library(stringr)
# genes_to_check=str_to_upper(unique(gene))
genes_to_check = gene

th=theme(axis.text.x = element_text(vjust = 0.5, hjust=0.5))
p_all_markers <- DotPlot(UMAP, features = genes_to_check,
                         assay='RNA' ,group.by = 'celltype' )  + coord_flip()+th
p_all_markers
data<-p_all_markers$data
data<-p_all_markers$data

colnames(data)
colnames(data)<-c("AverageExpression_unscaled","Precent Expressed","Features","celltype","Average Expression")
unique(data$`Precent Expressed`)

####用ggplot画图####
p = ggplot(data,aes(celltype,Features,size = `Precent Expressed` ))+
  geom_point(shape=21,aes(fill= `Average Expression`),position =position_dodge(0))+
  theme_minimal()+xlab(NULL)+ylab(NULL) +
  scale_size_continuous(range=c(1,10))+theme_bw()+
  scale_fill_gradient(low = "grey", high = "#E54924")+
  theme(legend.position = "right",legend.box = "vertical", #图例位置
        legend.margin=margin(t= 0.5, unit='cm'),
        legend.spacing = unit(0,"in"),
        axis.text.x  = element_text(color="black",size=15,angle = 90,
                                    vjust = 0.25, hjust=0.25),#x轴
        axis.text.y  = element_text(color="black",size=12),#y轴
        legend.text = element_text(size =10,color="black"),#图例
        legend.title = element_text(size =10,color="black"),#图例
        axis.title.y=element_text(vjust=1,
                                  size=16)
  )+labs(x=" ",y = "Features")
p
p1 = ggplot(data,aes(celltype,Features,size = `Precent Expressed` ))+
  geom_point(shape=21,aes(fill= `Average Expression`),position =position_dodge(0))+
  theme_minimal()+xlab(NULL)+ylab(NULL) +
  scale_size_continuous(range=c(1,10))+theme_bw()+
  scale_fill_gradient(low = "#E54924", high = "#498EA4")+
  theme(legend.position = "right",legend.box = "vertical",
        legend.margin=margin(t= 0, unit='cm'),
        legend.spacing = unit(0,"in"),
        axis.text.x  = element_text(color="black",size=15,angle = 90,
                                    vjust = 0.5, hjust=0.5),
        axis.text.y  = element_text(color="black",size=12),
        legend.text = element_text(size =10,color="black"),
        legend.title = element_text(size =10,color="black"),
        axis.title.y=element_text(vjust=1,
                                  size=16)
  )+labs(x=" ",y = "Features")

# 保存
library(patchwork)
p+p1
ggsave(file="02.dotplot_celltype.pdf",p,height= 10,width = 6)
ggsave(file="02.dotplot_celltype.png",p,height= 10,width = 6,units = 'in',dpi = 600)

# library(tidyverse)
# UMAP@meta.data$celltype=celltype[match(clusters,celltype$ClusterID),'celltype']
table(UMAP$celltype)

color2=palettes(category = "random",25,show_col=F)

df=UMAP@reductions$umap@cell.embeddings %>%
  as.data.frame() %>%
  cbind(celltype=UMAP@meta.data$celltype)

head(df)
# devtools::install_github('junjunlab/scRNAtoolVis')
# devtools::install_github("sajuukLyu/ggunchull", type = "source")

# install.packages('alphahull')
library(ggunchull)
library(Cairo)
# CairoPDF("tsne_ellipse.pdf", width=12, height=9) 
#用stat_ellipse加置信区间
p_tsne1=ggplot(df, aes(umap_1, umap_2, color=celltype))+
  geom_point(size = 0.5) + #细胞亚群点的大小
  scale_color_manual(values=color2)+
  theme(panel.border = element_blank(), #边框
        axis.title = element_blank(),  #轴标题
        axis.text = element_blank(), # 文本
        axis.ticks = element_blank(),
        # panel.grid = element_blank(),
        # panel.grid.major = element_blank(), #网格线
        # panel.grid.minor = element_blank(), #网格线
        panel.background = element_rect(fill = 'white'), #背景色
        plot.background=element_rect(fill="white"),
        legend.title = element_blank(), #去掉legend.title 
        legend.key=element_rect(fill='white'), 
        legend.text = element_text(size=14), #设置legend标签的大小
        legend.key.size=unit(0.6,'cm') )+
  stat_ellipse(aes(x = umap_1, y = umap_2, fill = celltype),
               geom = "polygon",#多边形
               linetype=2, # 椭圆，虚线
               alpha = 0.05,  #椭圆不透明度
               show.legend = FALSE, #去掉椭圆对应的图例
               level = 0.93)+ #置信区间
  guides(fill= guide_legend(override.aes = list(size = 2)))+ #图例点的大小
  scale_color_manual(values=color2)

celltype_position <- df %>%
  group_by(celltype) %>%
  dplyr::summarise(
    umap_1 = median(umap_1),
    umap_2 = median(umap_2));p_tsne1

library(ggrepel)
#install.packages('tidydr')
library(tidydr)
p_tsne2 <- p_tsne1+theme_dr(xlength=0.22,ylength=0.22,
                            arrow=grid::arrow(length=unit(0.15,'inches'),type='closed'))+
  theme(panel.grid = element_blank());p_tsne2



pdf(file="03.celltype_UMAP.pdf",width=6,height=6,family='Times')
p_tsne3<-p_tsne2+
  geom_label_repel(data = celltype_position,aes(label=celltype), 
                   fontface="bold",point.padding=unit(0.1, "lines"))+theme(legend.position = "none");p_tsne3

dev.off()


png(file="03.celltype_UMAP.png",width=6,height=6,family='Times',units='in',res=600)
p_tsne3<-p_tsne2+
  geom_label_repel(data = celltype_position,aes(label=celltype), 
                   fontface="bold",point.padding=unit(0.1, "lines"))+theme(legend.position = "none");p_tsne3

dev.off()

# 05.AUCell -----------------------------------------------------
rm(list = ls())
setwd('/home/data/project/30.JNPJ-6001-1/02_scRNA/CRC/')
if (! dir.exists("./04_AUCell")){
  dir.create("./04_AUCell")
}
setwd("./04_AUCell")
load('../03_celltype/UMAP.celltype.Rdata')

library(AUCell)
library(ggplot2)
library(patchwork)
library(dplyr)
library(Seurat)

#获得原始表达量作为matrix
cells_rankings <- AUCell_buildRankings(UMAP@assays$RNA$data) 
geneSets <- read.csv('../../../01_Target/ecm.geneset.csv')
geneSets <- geneSets$symbol
names(geneSets) <- 'ECM.score'
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings, aucMaxRank=nrow(cells_rankings)*0.1)

##set gene set of interest here for plotting
aucs <- as.numeric(getAUC(cells_AUC))
UMAP$AUCell  <- aucs
library(ggraph)
p <- ggplot(data.frame(UMAP@meta.data, UMAP@reductions$umap@cell.embeddings), aes(umap_1, umap_2, color=AUCell)
) + geom_point( size=1.5
) + scale_color_viridis(option="A")  + theme_light(base_size = 15)+labs(title = "AUCell score")+
  # scale_color_gradient(high = 'firebrick',low = '#498EA4')+
  theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"))+
  theme(plot.title = element_text(hjust = 0.5))
p
ggsave(filename = '01.AUcell.pdf',p,w=6,h=5)
ggsave(filename = '01.AUcell.png',p,w=6,h=5,units = 'in',dpi = 600,bg = 'white')

##高低评分组
aucs.group <- data.frame(cell = colnames(UMAP),auc.score = aucs)
aucs.group$group <- ifelse(aucs.group$auc.score > median(aucs.group$auc.score),'High.score','Low.score')
UMAP <- AddMetaData(UMAP,aucs.group$group,col.name = 'Score.group')
table(UMAP$Score.group)

p1 <- DimPlot(UMAP,reduction = 'umap',cols = c('red','blue'),group.by = 'Score.group',label = T)
p1
ggsave(filename = '02.AUcell.group.pdf',p1,w=6,h=5)
ggsave(filename = '02.AUcell.group.png',p1,w=6,h=5,units = 'in',dpi = 600,bg = 'white')

all.count <- data.frame(celltype=UMAP$celltype, sample=UMAP$Score.group)
library(reshape2)
plot.celltype = table(all.count$celltype, all.count$sample) %>% as.data.frame() %>% recast(Var1 ~ Var2)
colnames(plot.celltype)[1] = "Cell"
plot.celltype <- plot.celltype[,c(1,3,2)]
table(UMAP$Score.group)
plot.celltype = cbind(plot.celltype, apply(plot.celltype[-1], 2, proportions))
colnames(plot.celltype)[2:3] <- paste0('Count_',colnames(plot.celltype)[2:3])
colnames(plot.celltype)[4:5] = paste0('Proportion_',colnames(plot.celltype)[4:5])
plot.celltype <- plot.celltype[order(plot.celltype$Proportion_High.score,decreasing = T),]
write.table(plot.celltype, "CellCount.xls", row.names = F, col.names = T, quote = F, sep = "\t")

dat.celltype = plot.celltype[c(1,4:5)]
colnames(dat.celltype)[2:3] = gsub('Proprotion','',colnames(dat.celltype)[2:3])
dat.celltype[2:3] = 100 * dat.celltype[2:3]

dat.plot = melt(dat.celltype, id.vars = "Cell", variable.name = "Sample", value.name = "value")
dat.plot$Cell = factor(dat.plot$Cell, levels = plot.celltype$Cell)
library(IOBR)
library(ggplot2)
library(IOBR)
dat.plot <- dat.plot[order(dat.plot$value),]
dat.plot$Cell = factor(dat.plot$Cell, levels = plot.celltype$Cell)
library(IOBR)
dat.plot$label = round(dat.plot$value, 2) %>% paste0(.,"%")
str(dat.plot)
dat.plot$Sample <- factor(dat.plot$Sample,levels = c('Proportion_High.score','Proportion_Low.score'))
ggplot(data = dat.plot, mapping = aes(x = value, y = Cell)) + 
  geom_text(aes(x = value + 1, label = label), hjust = 0) +
  geom_bar(aes(fill = Cell), stat = "identity") + 
  facet_wrap(~Sample) +
  scale_fill_manual(values = palettes(category = "random",22,show_col=F)) + 
  xlim(c(0,70)) +
  xlab("Cell Fraction") + ylab("Cell Type") +
  theme_minimal(base_size = 16) + 
  theme(panel.border = element_rect(color = "grey60", fill = "transparent"), legend.position = "none", 
        panel.grid = element_blank())
# ggsave("03.CellProportion.png", width = 8, height = 5, units = "in", dpi = 300)
# ggsave("03.CellProportion.pdf", width = 8, height = 5, units = "in", dpi = 300)
system.time(save(UMAP, file = "UMAP.celltype.ecm.Rdata"))  


# 06.DEcell -----------------------------------------------------
rm(list = ls())
setwd('/home/data/project/30.JNPJ-6001-1/02_scRNA/CRC/')
if (! dir.exists("./05_DEcell")){
  dir.create("./05_DEcell")
}
setwd("./05_DEcell")

pacman::p_load(tidyverse, reshape2, ggpubr, ggplot2, ggsci)
load('../03_celltype/UMAP.celltype.Rdata')
##绘制散点图 wilcox比较
### 不同细胞比例不同样本
colnames(UMAP)
# all.count <- data.frame(celltype=scRNA.norm.umap$celltype,
#                         group=scRNA.norm.umap$group, 
#                         sample = substr(colnames(scRNA.norm.umap),1,10))
all.count <- data.frame(celltype = UMAP$celltype,
                        group = UMAP$group,
                        sample = UMAP$orig.ident)

table(all.count$celltype)
table(all.count$group)
table(all.count$sample)
plot.celltype <- table(all.count$celltype, all.count$sample ) %>% as.data.frame() %>% recast(Var1 ~ Var2)
colnames(plot.celltype)[1] = "Cell"
plot.celltype <- cbind(plot.celltype, apply(plot.celltype[-1], 2, proportions))
colnames(plot.celltype)
percent.celltype <- plot.celltype[,c(1,12:21)]
percent.celltype[2:11] = 100 *percent.celltype[2:11] 

dat.celltype <- melt(percent.celltype)
colnames(dat.celltype) <- c('Cell','Sample','Percentage')

##添加一下分组
dat.celltype$Group <- all.count$group[match(dat.celltype$Sample,all.count$sample)]

##wilcox检验
vln.df <- dat.celltype
table(vln.df$Group)
vln.df$Group <- factor(vln.df$Group,levels = c('Normal','Tumor'))
# vln.df <- vln.df[vln.df$symbol%in%de.stat$symbol,]
vln.df$Percentage
exp_plot <- ggplot(vln.df,aes(x = Group, y = Percentage, fill = Group)) +
  geom_violin(trim=F,color="black",aes(fill=Group)) + #绘制小提琴图, “color=”设置小提琴图的轮廓线的颜色(#不要轮廓可以设为white以下设为背景为白色，其实表示不要轮廓线)
  #"trim"如果为TRUE(默认值),则将小提琴的尾部修剪到数据范围。如果为FALSE,不修剪尾部。
  stat_boxplot(geom="errorbar", 
               width=0.1,
               position = position_dodge(0.9)) +
  geom_boxplot(width=0.4,
               position=position_dodge(0.9),
               outlier.shape = NA,fill = 'white')+ #绘制箱线图，此处width=0.1控制小提琴图中箱线图的宽窄
  scale_fill_manual(values= c("#5cb3cc","#ed556a"), name = "Group")+
  scale_colour_manual(values = c("#5cb3cc","#ed556a"))+
  labs(title="", x="", y = "Percentage",size=20) +
  scale_color_brewer(palette = 'Set2')+
  stat_compare_means(data = vln.df,
                     mapping = aes(group = Group),
                     label ="p.signif",
                     method = 'wilcox.test',
                     paired = F,
                     label.x = 1.4
  ) +
  theme_bw()+
  theme(plot.title = element_text(hjust =0.5,colour="black",face="bold",size=15),
        axis.text.x=element_text(angle=45,hjust=1,colour="black",face="bold",size=12), 
        axis.text.y=element_text(hjust=0.5,colour="black",face="bold",size=12), 
        axis.title.x=element_text(size=16,face="bold"),
        axis.title.y=element_text(size=16,face="bold"),
        legend.text=element_text(face = "bold", hjust = 0.5,colour="black", size=12),
        legend.title = element_text(face = "bold", size = 12),
        legend.position = "top",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  guides(fill='none')+facet_wrap(~Cell,ncol = 4,scales = 'free')
exp_plot

ggsave(filename = '01.DEcell.pdf',exp_plot,w=8,h=6)
ggsave(filename = '01.DEcell.png',exp_plot,w=8,h=6,dpi = 600)

# 循环生成每个细胞类型的图
cell_types <- unique(dat.celltype$Cell)
for (cell in cell_types) {
  vln.df <- dat.celltype[dat.celltype$Cell == cell,]
  
  exp_plot <- ggplot(vln.df, aes(x = Group, y = Percentage, fill = Group)) +
    geom_violin(trim = FALSE, color = "black", aes(fill = Group)) +
    stat_boxplot(geom = "errorbar", width = 0.1, position = position_dodge(0.9)) +
    geom_boxplot(width = 0.4, position = position_dodge(0.9), outlier.shape = NA, fill = 'white') +
    scale_fill_manual(values = c("#5cb3cc", "#ed556a"), name = "Group") +
    scale_colour_manual(values = c("#5cb3cc", "#ed556a")) +
    labs(title = "", x = "", y = "Percentage", size = 20) +
    scale_color_brewer(palette = 'Set2') +
    stat_compare_means(data = vln.df, 
                       mapping = aes(group = Group), 
                       label = "p.signif", 
                       method = 'wilcox.test', 
                       paired = FALSE, 
                       label.x = 1.4) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, colour = "black", face = "bold", size = 15),
          axis.text.x = element_text(angle = 45, hjust = 1, colour = "black", face = "bold", size = 12),
          axis.text.y = element_text(hjust = 0.5, colour = "black", face = "bold", size = 12),
          axis.title.x = element_text(size = 16, face = "bold"),
          axis.title.y = element_text(size = 16, face = "bold"),
          legend.text = element_text(face = "bold", hjust = 0.5, colour = "black", size = 12),
          legend.title = element_text(face = "bold", size = 12),
          legend.position = "top", 
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank()) +
    guides(fill = 'none') +
    facet_wrap(~Cell, ncol = 5, scales = 'free')
  
  # 保存每个细胞类型的图
  ggsave(filename = paste0(cell, "_DEcell.pdf"), exp_plot, width = 3, height = 4)
  ggsave(filename = paste0(cell, "_DEcell.png"), exp_plot, width = 3, height = 4, dpi = 600)
}



## 差异cell
table(UMAP$celltype)
UMAP.som <- subset(UMAP, subset = celltype %in% c('Myeloid'))

system.time(save(UMAP.som, file = "UMAP.main.celltype.Rdata"))  


#PART3 关键细胞确定后的分析（UC数据）-----------
# 08.cellchat -----------------------------------------------------
rm(list = ls())
setwd('/home/data/project/30.JNPJ-6001-1/02_scRNA/UC/')
if (! dir.exists("./07_cellchat")){
  dir.create("./07_cellchat")
}
setwd("./07_cellchat")

library(CellChat) 
library(Seurat)
library(ggplot2)
load('../03_celltype/UMAP.celltype.Rdata')
table(UMAP$group)
##提取疾病组
UMAP.som <- subset(UMAP,group == 'UC')
table(UMAP.som$group)
#读入seurat处理后的rds文件
scRNA<-UMAP.som
#data.input = as.data.frame(GetAssayData(subset(scRNA), slot='counts'))# raw count
data.input<-scRNA@assays$RNA$data # normalized data matrix

meta = scRNA@meta.data # a dataframe with rownames containing cell mata data
table(meta$celltype)
unique(meta$celltype)
save(data.input,meta,file='UC/01.inputdata.rda')

##CellChat的输入需要的matrix和meta我们已经准备好，下面开始创建
cellchat <- createCellChat(object = data.input, meta = meta, group.by = "celltype",)
cellchat <- addMeta(cellchat, meta = meta)##增加其他meta信息9
cellchat <- setIdent(cellchat, ident.use = "celltype") # 将 "labels" 设为默认细胞标记类型，这个可以根据自己的数据自定义
cellchat@idents<-droplevels(cellchat@idents,exclude=setdiff(levels(cellchat@idents),unique(cellchat@idents)))

groupSize <- as.numeric(table(cellchat@idents)) # 每组细胞的数量

##基于配受体分析的数据库
CellChatDB <- CellChatDB.human # 包括人和老鼠的，这里我们用人
showDatabaseCategory(CellChatDB)###作者提供了可视化的代码，可以看到该数据库中“Secreted Signaling”占比过半
library(tidyverse)
dplyr::glimpse(CellChatDB$interaction)  ###看一下CellChatDB的基本结构
CellChatDB_interaction <- CellChatDB$interaction
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # 我们这里使用“Secreted Signaling”部分做后续的细胞通讯分析
# CellChatDB.use <- CellChatDB
cellchat@DB <- CellChatDB.use

##表达数据做进一步预处理,节省算力
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)###首先识别过表达基因（配体——受体）
cellchat <- identifyOverExpressedInteractions(cellchat)###然后识别过表达配受体之间过表达的相互作用   （绕绕绕绕绕....）
# project gene expression data onto PPI network (optional)
cellchat <- projectData(cellchat, PPI.human)

##计算胞间通讯概率，预测通讯网络
cellchat <- computeCommunProb(cellchat,type = 'truncatedMean')#
#cellchat <- filterCommunication(cellchat, min.cells = 10)##过滤掉小于10个细胞的胞间通讯网络
##胞间通讯网络的输出代码
df.net <- subsetCommunication(cellchat)  
#df.net<-na.omit(df.net)
write.csv(df.net,'UC/01.df.net.csv',row.names = F)
##信号通路的水平进一步推测胞间通讯，计算聚合网络
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
cellchat.UC <- cellchat
saveRDS(cellchat,'UC/cellchat.rds')
##可视化细胞互作的结果
groupSize <- as.numeric(table(cellchat@idents))

#### circle ####
pdf('UC/01.net_circle_number.pdf',w=9,h=9,family='serif')
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
dev.off()

png('UC/01.net_circle_number.png',w=9,h=9,units='in',res=300,family='serif')
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
dev.off()
# 
pdf('UC//02.net_circle_weight.pdf',w=9,h=9,family='serif')
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weight/strength")
dev.off()

png('UC/02.net_circle_weight.png',w=9,h=9,units='in',res=300,family='serif')
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weight/strength")
dev.off()

# show all the significant interactions (L-R pairs)
table(UMAP.som$celltype)
#需要指定受体细胞和配体细胞
cellchat@data.signaling
p1 <- netVisual_bubble(cellchat, remove.isolate = T,sources.use = c('Myeloid'))
p2 <- netVisual_bubble(cellchat, remove.isolate = T,targets.use = c('Myeloid'))
pdf('UC/03.bubble.pdf',w=10,h=7)
p1+p2
dev.off()
png('UC/03.bubble.png',w=10,h=7,units = 'in',res = 600)
p1+p2
dev.off()


##提取对照组----
UMAP.som <- subset(UMAP,group == 'Health')
table(UMAP.som$group)
#读入seurat处理后的rds文件
scRNA<-UMAP.som
#data.input = as.data.frame(GetAssayData(subset(scRNA), slot='counts'))# raw count
data.input<-scRNA@assays$RNA$data # normalized data matrix

meta = scRNA@meta.data # a dataframe with rownames containing cell mata data
save(data.input,meta,file='Health/01.inputdata.rda')
table(meta$celltype)
unique(meta$celltype)
##CellChat的输入需要的matrix和meta我们已经准备好，下面开始创建
cellchat <- createCellChat(object = data.input, meta = meta, group.by = "celltype",)
cellchat <- addMeta(cellchat, meta = meta)##增加其他meta信息9
cellchat <- setIdent(cellchat, ident.use = "celltype") # 将 "labels" 设为默认细胞标记类型，这个可以根据自己的数据自定义
cellchat@idents<-droplevels(cellchat@idents,exclude=setdiff(levels(cellchat@idents),unique(cellchat@idents)))

groupSize <- as.numeric(table(cellchat@idents)) # 每组细胞的数量

##基于配受体分析的数据库
CellChatDB <- CellChatDB.human # 包括人和老鼠的，这里我们用人
showDatabaseCategory(CellChatDB)###作者提供了可视化的代码，可以看到该数据库中“Secreted Signaling”占比过半
library(tidyverse)
dplyr::glimpse(CellChatDB$interaction)  ###看一下CellChatDB的基本结构
CellChatDB_interaction <- CellChatDB$interaction
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # 我们这里使用“Secreted Signaling”部分做后续的细胞通讯分析
# CellChatDB.use <- CellChatDB
cellchat@DB <- CellChatDB.use

##表达数据做进一步预处理,节省算力
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)###首先识别过表达基因（配体——受体）
cellchat <- identifyOverExpressedInteractions(cellchat)###然后识别过表达配受体之间过表达的相互作用   （绕绕绕绕绕....）
# project gene expression data onto PPI network (optional)
cellchat <- projectData(cellchat, PPI.human)

##计算胞间通讯概率，预测通讯网络
cellchat <- computeCommunProb(cellchat,type = 'truncatedMean')#
#cellchat <- filterCommunication(cellchat, min.cells = 10)##过滤掉小于10个细胞的胞间通讯网络
##胞间通讯网络的输出代码
df.net <- subsetCommunication(cellchat)  
#df.net<-na.omit(df.net)
write.csv(df.net,'Health/01.df.net.csv',row.names = F)
##信号通路的水平进一步推测胞间通讯，计算聚合网络
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
##可视化细胞互作的结果
groupSize <- as.numeric(table(cellchat@idents))
# Health组分析完成后保存
cellchat.Health <- cellchat

saveRDS(cellchat.Health,file = "Health/cellchat.Health.rds")
#### circle ####
pdf('Health/01.net_circle_number.pdf',w=9,h=9,family='serif')
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
dev.off()

png('Health/01.net_circle_number.png',w=9,h=9,units='in',res=300,family='serif')
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
dev.off()
# 
pdf('Health/02.net_circle_weight.pdf',w=9,h=9,family='serif')
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weight/strength")
dev.off()

png('Health/02.net_circle_weight.png',w=9,h=9,units='in',res=300,family='serif')
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weight/strength")
dev.off()

# show all the significant interactions (L-R pairs)
p1 <- netVisual_bubble(cellchat, remove.isolate = T,sources.use = c('Myeloid'))
p2 <- netVisual_bubble(cellchat, remove.isolate = T,targets.use = c('Myeloid'))
pdf('Health/03.bubble.pdf',w=9,h=6)
p1+p2
dev.off()
png('Health/03.bubble.png',w=9,h=6,units = 'in',res = 600)
p1+p2
dev.off()


##合并做比较分析
library(ggplot2)
library(patchwork)

# 读取两个组的CellChat对象
cellchat.Health <- readRDS("Health/cellchat.Health.rds")
cellchat.UC <- readRDS("UC/cellchat.UC.rds")

# 合并CellChat对象
object.list <- list(
  Health = cellchat.Health,
  UC = cellchat.UC
)

cellchat.merge <- mergeCellChat(
  object.list,
  add.names = names(object.list),
  cell.prefix = TRUE
)

saveRDS(
  cellchat.merge,
  file = "CellChat_merged.rds"
)
levels(cellchat.Health@idents)
levels(cellchat.UC@idents)
# 比较整体通讯数量
p1 <- compareInteractions(
  cellchat.merge,
  group = c(1, 2),
  show.legend = FALSE,
  measure = "count"
)

# 比较整体通讯强度
p2 <- compareInteractions(
  cellchat.merge,
  group = c(1, 2),
  show.legend = FALSE,
  measure = "weight"
)

pdf("03.compareInteractions.pdf", width = 8, height = 4)
p1 + p2
dev.off()

png(
  "03.compareInteractions.png",
  width = 8,
  height = 4,
  units = "in",
  res = 300
)
p1 + p2
dev.off()
library(ggplot2)
library(reshape2)
library(dplyr)

# 提取矩阵
count.Health <- cellchat.merge@net$Health$count
count.UC <- cellchat.merge@net$UC$count

weight.Health <- cellchat.merge@net$Health$weight
weight.UC <- cellchat.merge@net$UC$weight

# 确保行列顺序一致
celltype <- levels(cellchat.Health@idents)

count.Health <- count.Health[celltype, celltype]
count.UC <- count.UC[celltype, celltype]

weight.Health <- weight.Health[celltype, celltype]
weight.UC <- weight.UC[celltype, celltype]

# 计算UC相对于Health的差值
count.diff <- count.UC - count.Health
weight.diff <- weight.UC - weight.Health
df.count <- melt(count.diff)

colnames(df.count) <- c(
  "Source",
  "Target",
  "Difference"
)

p1 <- ggplot(
  df.count,
  aes(x = Target, y = Source, fill = Difference)
) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = Difference), size = 4) +
  scale_fill_gradient2(
    low = "#519CBA",
    mid = "white",
    high = "#E05B5B",
    midpoint = 0
  ) +
  labs(
    x = "Target",
    y = "Source",
    fill = "Difference",
    title = "Difference in interaction numbers (UC - Health)"
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(
      angle = 45, hjust = 1, color = "black"
    ),
    axis.text.y = element_text(color = "black"),
    plot.title = element_text(hjust = 0.5)
  )
p1
ggsave(
  "04.diffInteraction_count_heatmap.pdf",
  p1, width = 6, height = 4
)

ggsave(
  "04.diffInteraction_count_heatmap.png",
  p1, width = 6, height = 4, dpi = 600
)
df.weight <- melt(weight.diff)

colnames(df.weight) <- c(
  "Source",
  "Target",
  "Difference"
)

p2 <- ggplot(
  df.weight,
  aes(x = Target, y = Source, fill = Difference)
) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(
    aes(label = sprintf("%.3f", Difference)),
    size = 4
  ) +
  scale_fill_gradient2(
    low = "#519CBA",
    mid = "white",
    high = "#E05B5B",
    midpoint = 0
  ) +
  labs(
    x = "Target",
    y = "Source",
    fill = "Difference",
    title = "Difference in interaction strength (UC - Health)"
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(
      angle = 45, hjust = 1, color = "black"
    ),
    axis.text.y = element_text(color = "black"),
    plot.title = element_text(hjust = 0.5)
  )
p2
ggsave(
  "05.diffInteraction_weight_heatmap.pdf",
  p2, width = 6, height = 4
)

ggsave(
  "05.diffInteraction_weight_heatmap.png",
  p2, width = 6, height = 4, dpi = 600
)

# 比较 Myeloid 相关的配体-受体互作
# Myeloid作为发送细胞
p1 <- netVisual_bubble(
  cellchat.merge,
  comparison = c(1, 2),
  sources.use = "Myeloid",
  remove.isolate = TRUE,
  angle.x = 45
)

# Myeloid作为接收细胞
p2 <- netVisual_bubble(
  cellchat.merge,
  comparison = c(1, 2),
  targets.use = "Myeloid",
  remove.isolate = TRUE,
  angle.x = 45
)

pdf("06.Myeloid_comparison.pdf", width = 12, height = 8)

p1 + p2

dev.off()

png(
  "06.Myeloid_comparison.png",
  width = 12,
  height = 8,
  units = "in",
  res = 600
)
p1 + p2
dev.off()
# 6. 比较不同信号通路的信息流
# 比较信号通路信息流
pdf("07.pathway_rank.pdf", width = 10, height = 8)
rankNet(
  cellchat.merge,
  mode = "comparison",
  stacked = TRUE,
  do.stat = TRUE
)

dev.off()
png(
  "07.pathway_rank.png",
  width = 10,
  height = 8,
  units = "in",
  res = 600
)
rankNet(
  cellchat.merge,
  mode = "comparison",
  stacked = TRUE,
  do.stat = TRUE
)
dev.off()
# 09.Trajectory----------------------------------------------------
rm(list = ls())
setwd('/home/data/project/30.JNPJ-6001-1/02_scRNA/UC/')
if (! dir.exists("./08_Trajectory")){
  dir.create("./08_Trajectory")
}
setwd("./08_Trajectory")

# 利⽤ UMAP 将关键细胞降维，将关键细胞更精细的分为不同亚群
load('../06_DEcell/UMAP.main.celltype.Rdata')
sce.maincell=UMAP.som
sce.maincell <- NormalizeData(sce.maincell, normalization.method = "LogNormalize", scale.factor = 10000)
sce.maincell <- FindVariableFeatures(sce.maincell, selection.method = 'vst', nfeatures = 2000)
sce.maincell <- ScaleData(sce.maincell)
sce.maincell <- RunPCA(sce.maincell, features = VariableFeatures(object = sce.maincell))

sce.maincell <- JackStraw(sce.maincell, num.replicate = 100, dims = 50)
sce.maincell <- ScoreJackStraw(sce.maincell, dims = 1:50)
JackStrawPlot(sce.maincell, dims = 1:50)
ElbowPlot(sce.maincell, ndims = 50)

sce.maincell <- FindNeighbors(sce.maincell, dims = 1:20)
sce.maincell <- FindClusters(sce.maincell, resolution = 0.4)
sce.maincell<-RunUMAP(sce.maincell,dims = 1:20)
# sce.maincell$seurat_clusters <- factor(sce.maincell$seurat_clusters,levels = c('0','1','2','3','4','5'))
# remotes::install_github("IOBR/IOBR")
library(IOBR)
table(UMAP.som$celltype)
DimPlot(sce.maincell, reduction = "umap", group.by = "seurat_clusters", label = T,cols = palettes(category = "random",22,show_col=F))+ggtitle('Subtype of Myeloid')
pdf(file = '01.subtype.pdf',w=6,h=5)
DimPlot(sce.maincell, reduction = "umap", group.by = "seurat_clusters", label = T,cols = palettes(category = "random",22,show_col=F))+ggtitle('Subtype of Myeloid')
dev.off()
png(file = '01.subtype.png',w=6,h=5,units = 'in',res = 300)
DimPlot(sce.maincell, reduction = "umap", group.by = "seurat_clusters", label = T,cols = palettes(category = "random",22,show_col=F))+ggtitle('Subtype of Myeloid')
dev.off()
system.time(save(sce.maincell, file = "maincell.subtype.Rdata"))

library(monocle)
library(devtools)
devtools::load_all('/home/data/files/pkg/monocle')
library(Seurat)
library(tidyverse)
# load('../07_subtype/maincell.subtype.Rdata')
Mono_tj = sce.maincell
Mono_matrix = GetAssayData(Mono_tj, slot = "count", assay = "RNA")
feature_ann = data.frame(gene_id=rownames(Mono_matrix),gene_short_name=rownames(Mono_matrix))
rownames(feature_ann) = rownames(Mono_matrix)
Mono_fd = new("AnnotatedDataFrame", data = feature_ann)
sample_ann = Mono_tj@meta.data
cell.type = Idents(Mono_tj) %>% data.frame
sample_ann = cbind(sample_ann, cell.type)
colnames(sample_ann)[ncol(sample_ann)] = "Subtype"
Mono_pd = new("AnnotatedDataFrame", data =sample_ann)

Mono.cds = newCellDataSet(Mono_matrix, phenoData =Mono_pd, featureData=Mono_fd, expressionFamily=negbinomial.size())
rm(Mono_tj)

Mono.cds = estimateSizeFactors(Mono.cds)
Mono.cds = estimateDispersions(Mono.cds)
WGCNA::collectGarbage()

#plot_ordering_genes(Mono.cds)
#8 9+3 6 1+2
Mono.cds = reduceDimension(Mono.cds, max_components = 2, verbose = T,
                           norm_method = "log"#, residualModelFormulaStr = "~nFeature_RNA+group"
)

Mono.cds = orderCells(Mono.cds)
save(Mono.cds, file = "mycds_order.RData")

# load('mycds_order_enod.RData')
table(Mono.cds$seurat_clusters)
# Mono.cds$seurat_clusters = factor(Mono.cds$seurat_clusters, levels = c('0','1','2'))
Mono.cds$seurat_clusters
Mono.cds$Pseudotime
p = plot_cell_trajectory(Mono.cds,color_by="Pseudotime",cell_size = 1, theta = 180,
                         size=1, show_backbone=TRUE, show_branch_points = F) +
  #  scale_color_manual(values = c("red","orange","green","blue","yellow",'brown','darkgreen','midnightblue','pink','purple')) +
  theme(text = element_text(size = 14))
p
ggsave("02.PSD_Trajectory.png", width = 6, height = 6, units = "in", dpi = 300)
ggsave("02.PSD_Trajectory.pdf", width = 6, height = 6, units = "in", dpi = 300)
p3 = plot_cell_trajectory(Mono.cds,color_by="State",cell_size = 1, theta = 180,
                          size=1, show_backbone=TRUE, show_branch_points = F) +
  scale_color_manual(values = c("red","orange","blue",'darkgreen','purple',"firebrick",'skyblue','pink','cyan','lightseagreen')) +
  theme(text = element_text(size = 14))
p3
ggsave("03.State_Trajectory.png", width = 6, height = 6, units = "in", dpi = 300)
ggsave("03.State_Trajectory.pdf", width = 6, height = 6, units = "in", dpi = 300)

p4 = plot_cell_trajectory(Mono.cds,color_by="State",cell_size = 1, theta = 180,
                          size=1, show_backbone=TRUE, show_branch_points = F) +
  scale_color_manual(values = c("red","orange","blue",'darkgreen','purple',"firebrick",'skyblue','pink','cyan','lightseagreen')) +
  theme(text = element_text(size = 14))+
  facet_wrap(~group,nrow = 1)
p4
ggsave(filename = '04.State_Trajectory(group).pdf',p4,w=10,h=6, units = "in", dpi = 300)
ggsave(filename = '04.State_Trajectory(group).png',p4,w=10,h=6, units = "in", dpi = 300)


# 10.ECMscore----------------------------------------------------
rm(list = ls())
setwd('/home/data/project/30.JNPJ-6001-1/02_scRNA/UC/')
if (! dir.exists("./09_ECMscore")){
  dir.create("./09_ECMscore")
}
setwd("./09_ECMscore")

load('../08_Trajectory/maincell.subtype.Rdata')
table(sce.maincell$seurat_clusters)
dat <- as.data.frame(sce.maincell@assays$RNA$counts)

perform.corplot <- function(dat1,dat2,type1,type2){
  x <- as.numeric(dat1[,1])
  y <- as.data.frame(dat2)
  library(psych)
  d <- corr.test(y,x,use="complete",method = 'spearman')
  r <- data.frame(d$r)
  p <- data.frame(d$p)
  correlation<-data.frame(rownames(p),r$d.r,p$d.p)
  colnames(correlation) <- c(type2,"Correlation","p.value")
  write.csv(correlation,file = paste0('correlation.',type2,'.csv'),row.names = F,quote = F)
  correlation<-correlation[correlation$p.value<0.05,]
  correlation$sig[correlation$p.value>0.05] <- "ns"   
  correlation$sig[correlation$p.value<0.05&correlation$p.value>0.01] <- "*"   
  correlation$sig[correlation$p.value<0.01&correlation$p.value>0.001] <- "**"  
  correlation$sig[correlation$p.value<0.001&correlation$p.value>0.0001] <- "***"   
  correlation$sig[correlation$p.value<0.0001] <- "****"
  correlation$'correlation'<-abs(correlation$Correlation)
  #correlation_sig<-correlation[c(3,9,10,11,23,29,7,12,18,31,1,14,15,25,28),]
  library(ggpubr)
  library("viridis")
  #trace(ggdotchart,edit=T)
  p0<-ggdotchart(correlation, x = type2, y = "Correlation",
                 
                 dot.size ='correlation',
                 color ='p.value',
                 sorting = "descending",
                 add = "segments",                             # 添加棒子
                 rotate = TRUE,
                 ggtheme = theme_pubr(),                        # 改变主题
                 ylab="Correlation Coefficient (r)",
                 xlab='',
                 # dot.size = 6 ,
                 title=paste0('Correlation between ',type1,' and ',type2)
                 
  )+scale_colour_gradient( high = "#4682B4",low = "#E05B5B")
  
  p10<-p0+theme(legend.position = "right",
                panel.background = element_blank())+geom_hline(aes(yintercept=0),linetype="dashed",lwd = 0.2)+
    theme(axis.title.x =element_text(size=15,family = "Times", face = "bold"),
          axis.text.x =element_text(size=12,family = "Times", face = "bold"),
          axis.title.y =element_text(size=20,family = "Times", face = "bold"),
          axis.text.y=element_text(size=12,family = "Times", face = "bold"),
          plot.title=element_text(size=20,family = "Times", face = "bold",hjust=0.5),
          legend.text = element_text(size = 16, family = "Times"),
          legend.title = element_text(size = 18, family = "Times",face = "bold"))+
    theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())
  
  p10
  height <- length(rownames(correlation))
  
  ggsave(filename = paste0('cor_',type2,'.pdf'),p10,w=12,h=height *0.25+2,limitsize = F)  
  ggsave(filename = paste0('cor_',type2,'.png'),p10,w=12,h=height *0.25+2,dpi = 600,limitsize = F)
  
}

##每个亚群单独做，先评估ECM score，得到一列数，再评估通路score，得到N列数，再进行Spearman相关性，绘制棒棒糖图。

clusters <- seq(0,8,1)
for (i in c(1:length(clusters))) {
  filename = paste0('0',i,'_Cluster',clusters[i])
  if (! dir.exists(filename)){
    dir.create(filename)
  }
  setwd(filename)
  
  cluster.id <- colnames(sce.maincell)[which(sce.maincell$seurat_clusters == clusters[i])]
  cluster.dat <- dat[,cluster.id]
  geneset<- read.csv('../../../../01_Target/ecm.geneset.csv')
  geneset$type <- 'ECM'
  length(unique(geneset$symbol))
  gene_list <- split(as.matrix(geneset)[,1],
                     geneset[,2])
  
  library(GSVA)
  ssgsea.param <- ssgseaParam(
    as.matrix(cluster.dat),
    gene_list,
    assay = NA_character_,
    annotation = NA_character_,
    minSize = 1,
    maxSize = Inf,
    alpha = 0.25,
    normalize = TRUE
  )
  score = gsva(ssgsea.param)
  write.csv(score,file = paste0('Cluster',clusters[i],'_ECMscore.csv'),row.names = T)
  
  # Hallmark------
  geneset<- read.gmt('../../../../../../files/GSVA/h.all.v2024.1.Hs.symbols.gmt')
  length(unique(geneset$gene))
  gene_list <- split(as.matrix(geneset)[,2],
                     geneset[,1])
  
  library(GSVA)
  ssgsea.param <- ssgseaParam(
    as.matrix(cluster.dat),
    gene_list,
    assay = NA_character_,
    annotation = NA_character_,
    minSize = 1,
    maxSize = Inf,
    alpha = 0.25,
    normalize = TRUE
  )
  score2 = gsva(ssgsea.param)
  
  write.csv(score2,file = paste0('Cluster',clusters[i],'_Hallmark.csv'),row.names = T)
  
  perform.corplot(dat1 = t(score),dat2 = t(score2),type1 = 'ECM',type2 = 'Hallmark')
  
  ##BIOCARTA--------
  geneset<- read.gmt('../../../../../../files/GSVA/c2.cp.biocarta.v2024.1.Hs.symbols.gmt')
  length(unique(geneset$gene))
  gene_list <- split(as.matrix(geneset)[,2],
                     geneset[,1])
  
  library(GSVA)
  ssgsea.param <- ssgseaParam(
    as.matrix(cluster.dat),
    gene_list,
    assay = NA_character_,
    annotation = NA_character_,
    minSize = 1,
    maxSize = Inf,
    alpha = 0.25,
    normalize = TRUE
  )
  score2 = gsva(ssgsea.param)
  
  write.csv(score2,file = paste0('Cluster',clusters[i],'_BIOCARTA.csv'),row.names = T)
  
  perform.corplot(dat1 = t(score),dat2 = t(score2),type1 = 'ECM',type2 = 'BIOCARTA')
  
  
  ##KEGG--------
  geneset<- read.gmt('../../../../../../files/GSVA/c2.cp.kegg.v7.5.1.symbols.gmt')
  length(unique(geneset$gene))
  gene_list <- split(as.matrix(geneset)[,2],
                     geneset[,1])
  
  library(GSVA)
  ssgsea.param <- ssgseaParam(
    as.matrix(cluster.dat),
    gene_list,
    assay = NA_character_,
    annotation = NA_character_,
    minSize = 1,
    maxSize = Inf,
    alpha = 0.25,
    normalize = TRUE
  )
  score2 = gsva(ssgsea.param)
  
  write.csv(score2,file = paste0('Cluster',clusters[i],'_KEGG.csv'),row.names = T)
  
  perform.corplot(dat1 = t(score),dat2 = t(score2),type1 = 'ECM',type2 = 'KEGG')
  
  ##Wikipathway--------
  geneset<- read.gmt('../../../../../../files/GSVA/c2.cp.wikipathways.v2024.1.Hs.symbols.gmt')
  length(unique(geneset$gene))
  gene_list <- split(as.matrix(geneset)[,2],
                     geneset[,1])
  
  library(GSVA)
  ssgsea.param <- ssgseaParam(
    as.matrix(cluster.dat),
    gene_list,
    assay = NA_character_,
    annotation = NA_character_,
    minSize = 1,
    maxSize = Inf,
    alpha = 0.25,
    normalize = TRUE
  )
  score2 = gsva(ssgsea.param)
  
  write.csv(score2,file = paste0('Cluster',clusters[i],'_Wikipathway.csv'),row.names = T)
  
  # perform.corplot(dat1 = t(score),dat2 = t(score2),type1 = 'ECM',type2 = 'Wikipathway')
  
  ##GOBP--------
  geneset<- read.gmt('../../../../../../files/GSVA/c5.go.bp.v2024.1.Hs.symbols.gmt')
  length(unique(geneset$gene))
  gene_list <- split(as.matrix(geneset)[,2],
                     geneset[,1])
  
  library(GSVA)
  ssgsea.param <- ssgseaParam(
    as.matrix(cluster.dat),
    gene_list,
    assay = NA_character_,
    annotation = NA_character_,
    minSize = 1,
    maxSize = Inf,
    alpha = 0.25,
    normalize = TRUE
  )
  score2 = gsva(ssgsea.param)
  
  write.csv(score2,file = paste0('Cluster',clusters[i],'_GOBP.csv'),row.names = T)
  
  # perform.corplot(dat1 = t(score),dat2 = t(score2),type1 = 'ECM',type2 = 'GOBP')
  
  setwd('../')
}

#ALL-------
filename = 'All'
if (! dir.exists(filename)){
  dir.create(filename)
}
setwd(filename)
geneset<- read.csv('../../../../01_Target/ecm.geneset.csv')
geneset$type <- 'ECM'
length(unique(geneset$symbol))
gene_list <- split(as.matrix(geneset)[,1],
                   geneset[,2])

library(GSVA)
ssgsea.param <- ssgseaParam(
  as.matrix(dat),
  gene_list,
  assay = NA_character_,
  annotation = NA_character_,
  minSize = 1,
  maxSize = Inf,
  alpha = 0.25,
  normalize = TRUE
)
score = gsva(ssgsea.param)
write.csv(score,file = 'ECMscore.csv',row.names = T)

# Hallmark------
geneset<- read.gmt('../../../../../../files/GSVA/h.all.v2024.1.Hs.symbols.gmt')
length(unique(geneset$gene))
gene_list <- split(as.matrix(geneset)[,2],
                   geneset[,1])

library(GSVA)
ssgsea.param <- ssgseaParam(
  as.matrix(dat),
  gene_list,
  assay = NA_character_,
  annotation = NA_character_,
  minSize = 1,
  maxSize = Inf,
  alpha = 0.25,
  normalize = TRUE
)
score2 = gsva(ssgsea.param)

write.csv(score2,file = 'Hallmark.csv',row.names = T)

perform.corplot(dat1 = t(score),dat2 = t(score2),type1 = 'ECM',type2 = 'Hallmark')

##BIOCARTA--------
geneset<- read.gmt('../../../../../../files/GSVA/c2.cp.biocarta.v2024.1.Hs.symbols.gmt')
length(unique(geneset$gene))
gene_list <- split(as.matrix(geneset)[,2],
                   geneset[,1])

library(GSVA)
ssgsea.param <- ssgseaParam(
  as.matrix(dat),
  gene_list,
  assay = NA_character_,
  annotation = NA_character_,
  minSize = 1,
  maxSize = Inf,
  alpha = 0.25,
  normalize = TRUE
)
score2 = gsva(ssgsea.param)

write.csv(score2,file = 'BIOCARTA.csv',row.names = T)

perform.corplot(dat1 = t(score),dat2 = t(score2),type1 = 'ECM',type2 = 'BIOCARTA')

##KEGG--------
geneset<- read.gmt('../../../../../../files/GSVA/c2.cp.kegg.v7.5.1.symbols.gmt')
length(unique(geneset$gene))
gene_list <- split(as.matrix(geneset)[,2],
                   geneset[,1])

library(GSVA)
ssgsea.param <- ssgseaParam(
  as.matrix(dat),
  gene_list,
  assay = NA_character_,
  annotation = NA_character_,
  minSize = 1,
  maxSize = Inf,
  alpha = 0.25,
  normalize = TRUE
)
score2 = gsva(ssgsea.param)
write.csv(score2,file = 'KEGG.csv',row.names = T)
perform.corplot(dat1 = t(score),dat2 = t(score2),type1 = 'ECM',type2 = 'KEGG')

##Wikipathway--------
geneset<- read.gmt('../../../../../../files/GSVA/c2.cp.wikipathways.v2024.1.Hs.symbols.gmt')
length(unique(geneset$gene))
gene_list <- split(as.matrix(geneset)[,2],
                   geneset[,1])

library(GSVA)
ssgsea.param <- ssgseaParam(
  as.matrix(dat),
  gene_list,
  assay = NA_character_,
  annotation = NA_character_,
  minSize = 1,
  maxSize = Inf,
  alpha = 0.25,
  normalize = TRUE
)
score2 = gsva(ssgsea.param)

write.csv(score2,file = 'Wikipathway.csv',row.names = T)

# perform.corplot(dat1 = t(score),dat2 = t(score2),type1 = 'ECM',type2 = 'Wikipathway')

##GOBP--------
geneset<- read.gmt('../../../../../../files/GSVA/c5.go.bp.v2024.1.Hs.symbols.gmt')
length(unique(geneset$gene))
gene_list <- split(as.matrix(geneset)[,2],
                   geneset[,1])

library(GSVA)
ssgsea.param <- ssgseaParam(
  as.matrix(dat),
  gene_list,
  assay = NA_character_,
  annotation = NA_character_,
  minSize = 1,
  maxSize = Inf,
  alpha = 0.25,
  normalize = TRUE
)
score2 = gsva(ssgsea.param)

write.csv(score2,file = '_OBP.csv',row.names = T)

# perform.corplot(dat1 = t(score),dat2 = t(score2),type1 = 'ECM',type2 = 'GOBP')

setwd('../')

