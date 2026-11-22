rm(list = ls()) 
setwd('/home/data/project/30.JNPJ-6001-1/')
file.dir = '10_scRNA_part2'
if (! dir.exists(file.dir)){
  dir.create(file.dir)
}
setwd(file.dir)

#A 组间表达-------
setwd('/home/data/project/30.JNPJ-6001-1/10_scRNA_part2')
file.dir = '01_expression'
if (! dir.exists(file.dir)){
  dir.create(file.dir)
}
setwd(file.dir)

load('../../02_scRNA/UC/08_Trajectory/maincell.subtype.Rdata')
pacman::p_load(reshape2,tidyverse,Seurat,rstatix,ggpubr,gridExtra)
UMAP.som <- sce.maincell
table(UMAP.som$celltype)
som.scobj <- UMAP.som
table(som.scobj$celltype)
gene <- read.csv('../../06_113model/bestModel.gene.csv')
VlnPlot(UMAP.som,features = gene$symbol,group.by = 'group')
##加个差异分析----
library(reshape2)
library(tidyverse)
vln.df <- as.data.frame(UMAP.som@assays$RNA$data[gene$symbol[-6], ]) %>% 
  apply(2, function(x) {
    if (sum(x) == 0) {
      return(rep(0, length(x)))  # 如果总和为0，返回一个全0向量
    } else {
      return(x / sum(x) * 10000)  # 正常执行归一化操作
    }
  })
rownames(vln.df) <- gene$symbol[-6]
vln.df <- log2(vln.df + 1) %>% as.data.frame()

vln.df$symbol=rownames(vln.df)
vln.df=melt(vln.df,id="symbol")
colnames(vln.df)[c(2,3)]=c("sample","exp")
head(vln.df)

anno=UMAP.som@meta.data[,c("group","celltype")]
anno$sample <- rownames(anno)
vln.df=inner_join(vln.df,anno,by="sample")
colnames(vln.df)
library(rstatix)
stat.test<-vln.df%>%
  group_by(symbol)%>%
  wilcox_test(exp ~ group)%>%
  adjust_pvalue(method = 'fdr')

library(ggpubr)
library(ggplot2)
library(ggsci)
vln.df$group <- factor(vln.df$group,levels = c('UC','Health'))

exp_plot <- ggplot(vln.df,aes(x = group, y = exp, fill = group)) +
  geom_violin(trim=F,color="black") +
  stat_boxplot(geom="errorbar",
               width=0.1,
               position = position_dodge(0.9)) +
  geom_boxplot(width=0.3,
               position=position_dodge(0.9),
               outlier.shape = NA,fill = 'white')+ #绘制箱线图，此处width=0.1控制小提琴图中箱线图的宽窄
  scale_fill_manual(values= c("#CD3700","#4682B4"), name = "Group")+
  labs(title="", x="", y = "Expression level",size=20) +
  
  stat_compare_means(data = vln.df,
                     mapping = aes(group = group),
                     label ="p.signif",
                     method = 'wilcox.test',
                     paired = F,label.x = 1.4) +
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
  guides(fill='none')+facet_wrap(~symbol + ,ncol = 3,scales = 'free')
exp_plot
ggsave(filename = '01.expression(group).pdf',exp_plot,w=6,h=5)
ggsave(filename = '01.expression(group).png',exp_plot,w=6,h=5,units = 'in',dpi = 600,bg = 'white')

vln.df$seurat_clusters <- sce.maincell$seurat_clusters[match(vln.df$sample,colnames(sce.maincell))]
table(vln.df$seurat_clusters)
vln.df$seurat_clusters <- paste0('Cluster',vln.df$seurat_clusters)

exp_plot <- ggplot(vln.df,aes(x = group, y = exp, fill = group)) +
  geom_violin(trim=F,color="black") + 
  stat_boxplot(geom="errorbar",
               width=0.1,
               position = position_dodge(0.9)) +
  geom_boxplot(width=0.3,
               position=position_dodge(0.9),
               outlier.shape = NA,fill = 'white')+ #绘制箱线图，此处width=0.1控制小提琴图中箱线图的宽窄
  scale_fill_manual(values= c("#CD3700","#4682B4"), name = "Group")+
  labs(title="", x="", y = "Expression level",size=20) +
  
  stat_compare_means(data = vln.df,
                     mapping = aes(group = group),
                     label ="p.signif",
                     method = 'wilcox.test',
                     paired = F,label.x = 1.4) +
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
  guides(fill='none')+facet_wrap(~symbol + seurat_clusters ,ncol = 9,scales = 'free')
exp_plot
ggsave(filename = '02.expression(subtype).pdf',exp_plot,w=20,h=15)
ggsave(filename = '02.expression(subtype).png',exp_plot,w=20,h=15,units = 'in',dpi = 600,bg = 'white')


#B 拟时序动态-------
setwd('/home/data/project/30.JNPJ-6001-1/10_scRNA_part2')
file.dir = '02_Trajectory'
if (! dir.exists(file.dir)){
  dir.create(file.dir)
}
setwd(file.dir)

load('../../02_scRNA/UC/08_Trajectory/mycds_order.RData')
gene <- read.csv('../../06_113model/bestModel.gene.csv')
allgene <- data.frame(rownames(Mono.cds))
gene <- read.csv('../../01_Target/ecm.geneset.csv')

p5 <- plot_genes_in_pseudotime(Mono.cds[gene$symbol[-1],],
                               color_by = "State",
                               ncol = 3)+
  theme(axis.title = element_text(face = "bold"),
        strip.text = element_text(face = "bold"))
p5
ggsave(filename = '01.sdheat_hubgene.pdf',p5,w=6,h=4)
ggsave(filename = '01.sdheat_hubgene.png',p5,w=6,h=4,units = 'in',dpi = 300)
