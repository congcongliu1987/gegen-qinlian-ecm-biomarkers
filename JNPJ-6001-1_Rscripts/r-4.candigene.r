rm(list = ls())
setwd('/home/data/project/30.JNPJ-6001-1/')
file.dir = '04_candigene'
if (! dir.exists(file.dir)){
  dir.create(file.dir)
}
setwd(file.dir)

sig.diff <- read.csv('../03_DEGs/DEG_sig.csv',row.names = 1)
ecm <- read.csv('../01_Target/ecm.geneset.csv')

candigene <- data.frame(symbol = intersect(rownames(sig.diff),ecm$symbol))
write.csv(candigene,file = 'candigene.csv',row.names = F,quote = F)
library(ggvenn)

mydata<-list('DEGs'= rownames(sig.diff),'ECM'= ecm$symbol)
venn_plot <- ggvenn(mydata, c('DEGs','ECM'),
                    fill_color = c("#f8bc31", '#5CB0C3'),
                    show_percentage = F,
                    stroke_alpha = 0.5,
                    stroke_size = 0.3,
                    text_size = 4,
                    set_name_size = 5,
                    stroke_color = "black",
                    stroke_linetype = "solid",
                    set_name_color = c("#f8bc31", '#5CB0C3'),
                    text_color = 'black')

print(venn_plot)
# 保存为 PDF
pdf('01.venn.pdf', w = 8, h = 6)
print(venn_plot)
dev.off()

png('01.venn.png',w=8,h=6,units = 'in',res = 600)
print(venn_plot)
dev.off()
