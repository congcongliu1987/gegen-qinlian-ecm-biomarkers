rm(list = ls())
setwd('/home/data/project/30.JNPJ-6001-1/')
file.dir = '03_DEGs'
if (! dir.exists(file.dir)){
  dir.create(file.dir)
}
setwd(file.dir)
library(magrittr)
library(stringr)
library(limma)
df = read.csv('../00_rawdata/dat.GSE87466.csv',row.names = 1) 
df.group = read.csv('../00_rawdata/group.GSE87466.csv')
df = df[,df.group$sample]
table(df.group$group)
df.group$group <- ifelse(df.group$group == 'Normal','Normal','UC')

df.group$group = factor(df.group$group, levels = c("Normal", "UC"))
design.mat = cbind(Normal = ifelse(df.group$group == "Normal", 1, 0), 
                   UC = ifelse(df.group$group == "Normal", 0, 1))

contrast.mat = makeContrasts(contrasts="UC-Normal", levels=design.mat)
fit = lmFit(df, design.mat)
fit = contrasts.fit(fit, contrast.mat)
fit = eBayes(fit)
fit = topTable(fit, coef = 1, number = Inf, adjust.method = "fdr")
DEG=na.omit(fit)
logFC_cutoff <- 1
DEG$change = as.factor(
  ifelse(DEG$adj.P.Val <0.05 & abs(DEG$logFC) > logFC_cutoff,
         ifelse(DEG$logFC > logFC_cutoff ,'UP','DOWN'),'NOT')
)
sig_diff <- subset(DEG,
                   DEG$adj.P.Val < 0.05 & abs(DEG$logFC) > logFC_cutoff)

dim(DEG)
dim(sig_diff)
summary(sig_diff$change) 
write.csv(DEG,file = "DEG_all.csv",row.names = T)
write.csv(sig_diff,file = "DEG_sig.csv",row.names = T)




ecm <- read.csv('../01_Target/ecm.geneset.csv')

diff.ecm <- data.frame(symbol = intersect(ecm$symbol,rownames(sig_diff)))
## 18

library(ggplot2)
library(ggthemes)
library(RColorBrewer)
library(scales)
library(ggrepel)
dat_rep<-DEG[rownames(DEG)%in%
               rownames(rbind(head(sig_diff[order(sig_diff$logFC,decreasing = T),],10),
                              head(sig_diff[order(sig_diff$logFC,decreasing = F),],10))),]
volcano_plot<- ggplot(data = DEG, 
                      aes(x = logFC,
                          y = -log10(adj.P.Val), 
                          color =change)) +
  scale_color_manual(values = c("skyblue4", "darkgray","firebrick2")) +
  geom_point(size = 1.2, alpha = 0.4, na.rm=T) +
  theme_bw(base_size = 12, base_family = "Times") +
  geom_vline(xintercept = c(-0.5,0.5),
             lty = 4,
             col = "darkgray",
             lwd = 0.6)+
  geom_hline(yintercept = -log10(0.05),
             lty = 4,
             col = "darkgray",
             lwd = 0.6)+
  theme(legend.position = "right",
        panel.grid = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(face="bold",
                                   color="black",
                                   family = "Times",
                                   size=13),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(face = "bold",
                                   color = "black",
                                   size = 15),
        axis.text.y = element_text(face = "bold",
                                   color = "black",
                                   size = 15),
        axis.title.x = element_text(face = "bold",
                                    color = "black",
                                    size = 15),
        axis.title.y = element_text(face = "bold",
                                    color = "black",
                                    size = 15)) +
  geom_label_repel(
    data = dat_rep,
    aes(label = rownames(dat_rep)),
    max.overlaps = 20,
    size = 4,
    box.padding = unit(0.5, "lines"),
    min.segment.length = 0,
    point.padding = unit(0.8, "lines"), segment.color = "black", show.legend = FALSE )+
  labs(x = "log2(Fold Change)",
       y = "-log10 (adj.P.Val)")
volcano_plot
ggsave('01.volcano.png', volcano_plot,width = 8, height = 6,units = 'in',dpi = 600)
ggsave('01.volcano.pdf', volcano_plot,width = 8, height = 6)

# 密度热图------
library(ComplexHeatmap)

rt <- df
rt <- rt[rownames(dat_rep),]
group <- df.group

x<-rt
mat <- t(scale(t(x)))#归一化
df1 <- as.data.frame(mat)
range(mat)

pdf('02.heatmap.pdf',  w=6,h=6,family='Times')
densityHeatmap(mat ,title = "Distribution as heatmap", ylab = " ",height = unit(3, "cm")) %v%
  HeatmapAnnotation(Group = group$group, col = list(Group = c("UC" = "#B72230", "Normal" = "#104680"))) %v%
  Heatmap(mat, 
          row_names_gp = gpar(fontsize = 9),
          show_column_names = F,
          show_row_names = T,
          ###show_colnames = FALSE,
          name = "expression", 
          ###cluster_cols = F,
          cluster_rows = T,
          height = unit(6, "cm"),
          #cluster_columns = FALSE,
          ###cluster_rows = FALSE,
          use_raster = T,
          col = colorRampPalette(c("#104680", "white","#D80305"))(100))
dev.off()

png('02.heatmap.png',w=6,h=6,units='in',res=600,family='Times')
densityHeatmap(mat ,title = "Distribution as heatmap", ylab = " ",height = unit(3, "cm")) %v%
  HeatmapAnnotation(Group = group$group, col = list(Group = c("UC" = "#B72230", "Normal" = "#104680"))) %v%
  Heatmap(mat, 
          row_names_gp = gpar(fontsize = 9),
          show_column_names = F,
          show_row_names = T,
          ###show_colnames = FALSE,
          name = "expression", 
          ###cluster_cols = F,
          cluster_rows = T,
          height = unit(6, "cm"),
          #cluster_columns = FALSE,
          ###cluster_rows = FALSE,
          use_raster = T,
          col = colorRampPalette(c("#104680", "white","#D80305"))(100))
dev.off()
