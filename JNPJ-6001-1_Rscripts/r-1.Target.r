rm(list = ls())
setwd('/home/data/project/30.JNPJ-6001-1/')
file.dir = '01_Target'
if (! dir.exists(file.dir)){
  dir.create(file.dir)
}
setwd(file.dir)

# 黄芩、葛根、黄连、⽢草

#匹配一下herb的tcmsp id
drug2mol <- readxl::read_xlsx('drug2mol.xlsx')
colnames(drug2mol)[1] <- 'herbname'
table(drug2mol$herbname)

library(readxl)
mol2target.all <- readxl::read_xlsx('../../../files/TCMSP/mol2target.xlsx')
mol2target <- mol2target.all[which(mol2target.all$`Mol ID` %in% drug2mol$`Mol ID`),]
mol2target <- mol2target[,c(1:3)]
length(unique(mol2target$`Target name`))

target.name <- data.frame(protein = unique(mol2target$`Target name`)) 

unique.pro <- read_xlsx('../../../files/TCMSP/uniport.xlsx')
unique.pro <- unique.pro[!duplicated(unique.pro$protein),]
map <- unique.pro[which(unique.pro$protein %in% target.name$protein),]
library(tidyverse)
unmap <- data.frame(protein = target.name[!target.name$protein %in% unique.pro$protein,]) %>% na.omit()

unique(drug2mol$herbname)
target.map <- read_xlsx('unmap.pro.xlsx')
target.map <- rbind(target.map,map)
colnames(mol2target)[3] <- 'protein'
mol2target <- merge(mol2target,target.map,by = 'protein')
colnames(mol2target)
mol2target <- mol2target[,c(2,3,1,4)]

drug2mol <- drug2mol[,c(1:3)]
allnet <- merge(drug2mol,mol2target,by = 'Mol ID')

colnames(allnet)
allnet <- allnet[,-4]
allnet$link <- paste0(allnet$`Mol ID`,allnet$herbname,allnet$symbol)
allnet <- allnet[!duplicated(allnet$link),]
allnet <- allnet[,-6]
write.csv(allnet,file = 'allnet.csv')

##画图的
cordrug2mol <- allnet[,c(2,3)]
cordrug2mol$link <- paste0(cordrug2mol$herbname,cordrug2mol$`Molecule Name.x`)
cordrug2mol <- cordrug2mol[!duplicated(cordrug2mol$link),]
cormol2tar <- allnet[,c(3,5)]
cormol2tar$link <- paste0(cormol2tar$`Molecule Name.x`,cormol2tar$symbol)
cormol2tar <- cormol2tar[!duplicated(cormol2tar$link),]
write.csv(cordrug2mol,file = 'cordrug2mol.csv',row.names = F)
write.csv(cormol2tar,file = 'cormol2tar.csv',row.names = F)

ecm <- readxl::read_xlsx('../../../files/Geneset/ecm.xlsx') %>% as.data.frame()
colnames(ecm)
ecm <- ecm[,-c(1:4)]
class(ecm) 

list <- list()
for (i in c(1:length(colnames(ecm)))) {
  list[[i]] <- data.frame(Term = colnames(ecm)[i],na.omit(ecm[,i]))
}

ecm.all <- do.call(rbind,list)

table(ecm.all$Term)
colnames(ecm.all)[2] <- 'Symbol'

write.csv(ecm.all,file = 'ecm.all.csv',row.names = F)
length(unique(ecm.all$Symbol))

ecm.geneset <- data.frame(symbol = intersect(ecm.all$Symbol,unique(allnet$symbol)))

write.csv(ecm.geneset,file = 'ecm.geneset.csv',row.names = F,quote = F)

library(ggvenn)

mydata<-list('All ECM Genes'= ecm.all$Symbol,'Target'= unique(allnet$symbol))
venn_plot <- ggvenn(mydata, c('All ECM Genes','Target'),
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
