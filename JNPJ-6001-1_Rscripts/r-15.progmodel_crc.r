rm(list = ls())
setwd('/home/data/project/30.JNPJ-6001-1/')
file.dir = '15_progmodel_CRC'
if (! dir.exists(file.dir)){
  dir.create(file.dir)
}
setwd(file.dir)
train_phenotype <- read.csv('../11_CRCsubtype/04_clinical/clinical_subtype.csv') 

train_phenotype$OS<-as.numeric(train_phenotype$OS)
train_phenotype$OS.time<-as.numeric(train_phenotype$OS.time)
train_phenotype2<-train_phenotype
colnames(train_phenotype2)

table(train_phenotype2$T.stage)
train_phenotype2$T.stage <- gsub('T','',train_phenotype2$T.stage)

table(train_phenotype2$N.stage)
train_phenotype2$N.stage <- gsub('N','',train_phenotype2$N.stage)

table(train_phenotype2$M.stage)
train_phenotype2$M.stage <- gsub('M','',train_phenotype2$M.stage)

table(train_phenotype2$Age)
train_phenotype2$Age<-ifelse(train_phenotype2$Age=='>65','1','0')


table(train_phenotype2$Gender)
train_phenotype2$Gender<-ifelse(train_phenotype2$Gender=='female','1','0')


colnames(train_phenotype2)
dat<-read.csv('../00_rawdata/dat.fpkm.tcga.csv',row.names = 1,check.names = F)
dat <- log2(dat +1)
hubgene <- read.csv('../06_113model/bestModel.gene.csv')
sub_dat <- dat[hubgene$symbol,]
sub_dat <- t(sub_dat) %>% as.data.frame() %>% rownames_to_column(var = 'sample')
train_phenotype3 <- merge(sub_dat,
                          train_phenotype2,
                          by = "sample")

write.csv(train_phenotype3,file = "clinical_gene.csv",row.names = F,quote = F)

train_risk_clinical <- train_phenotype3
train_risk_clinical = subset(train_risk_clinical, select = -c(sample))
dim(train_risk_clinical)
colnames_train <- colnames(train_risk_clinical)
covariates_train <- colnames_train[-which(colnames_train %in% c("OS", "OS.time"))]

train_risk_clinical$T.stage <- factor(train_risk_clinical$T.stage)
train_risk_clinical$N.stage <- factor(train_risk_clinical$N.stage)
library(survival)

res.MMP3 = coxph(Surv(time = OS.time, event = OS) ~ MMP3, data = train_risk_clinical) %>% summary
res.MMP3 = c(res.MMP3$conf.int[-2], res.MMP3$coefficients[5])

res.VCAM1 = coxph(Surv(time = OS.time, event = OS) ~ VCAM1, data = train_risk_clinical) %>% summary
res.VCAM1 = c(res.VCAM1$conf.int[-2], res.VCAM1$coefficients[5])

res.MMP1 = coxph(Surv(time = OS.time, event = OS) ~ MMP1, data = train_risk_clinical) %>% summary
res.MMP1 = c(res.MMP1$conf.int[-2], res.MMP1$coefficients[5])

res.PCOLCE = coxph(Surv(time = OS.time, event = OS) ~ PCOLCE, data = train_risk_clinical) %>% summary
res.PCOLCE = c(res.PCOLCE$conf.int[-2], res.PCOLCE$coefficients[5])

res.PPARG = coxph(Surv(time = OS.time, event = OS) ~ PPARG, data = train_risk_clinical) %>% summary
res.PPARG = c(res.PPARG$conf.int[-2], res.PPARG$coefficients[5])

res.KDR = coxph(Surv(time = OS.time, event = OS) ~ KDR, data = train_risk_clinical) %>% summary
res.KDR = c(res.KDR$conf.int[-2], res.KDR$coefficients[5])

res.age = coxph(Surv(time = OS.time, event = OS) ~ Age, data = train_risk_clinical) %>% summary
res.age = c(res.age$conf.int[-2], res.age$coefficients[5])

res.gender = coxph(Surv(time = OS.time, event = OS) ~ Gender, data = train_risk_clinical) %>% summary
res.gender = c(res.gender$conf.int[-2], res.gender$coefficients[5])

res.tstage = coxph(Surv(time = OS.time, event = OS) ~ T.stage, data = train_risk_clinical) %>% summary
res.tstage = cbind(res.tstage$conf.int[,-2], res.tstage$coefficients[,5])

res.nstage = coxph(Surv(time = OS.time, event = OS) ~ N.stage, data = train_risk_clinical) %>% summary
res.nstage = cbind(res.nstage$conf.int[,-2], res.nstage$coefficients[,5])

res.mstage = coxph(Surv(time = OS.time, event = OS) ~ M.stage, data = train_risk_clinical) %>% summary
res.mstage = c(res.mstage$conf.int[-2], res.mstage$coefficients[5])

res.ref = c(1,1,1,NA)
res = rbind(res.MMP1,res.MMP3,res.KDR,res.PCOLCE,res.PPARG,res.VCAM1,res.age,res.gender,res.ref,res.tstage,res.ref,res.nstage,res.mstage) %>% as.data.frame()
rownames(res)
res$Indicators = c("MMP1","MMP3","KDR","PCOLCE","PPARG","VCAM1","Age",'Gender','T1(Reference)','T2',
                   'T3','T4','N0(Reference)',"N1",'N2','M1 vs.M0')
colnames(res) = c("hr","low","up","pv","Indicator")
res$p = signif(res$pv, 2) %>% paste0("p = ", .)
res$p[is.na(res$pv)] = NA
res$Indicator = factor(res$Indicator, levels = rev(res$Indicator))
rownames(res) <- res$Indicator
res2 <- data.frame(p.value=res$pv,
                   HR=res$hr,
                   HR.95L=res$low,
                   HR.95H=res$up,
                   Indicator=res$Indicator)
rownames(res2) <- res2$Indicator
write.table(res2, file = "univariate_cox_prog_forest.xls", sep = "\t", quote = F)
res2 <- subset(res2, select = -c(Indicator))
library(tidyr)
hz <- paste(round(res2$HR,3),
            "(",round(res2$HR.95L,3),
            "-",round(res2$HR.95H,3),")",sep = "")
hz
hz[c(9,13)] <- ""
tabletext <- cbind(c(NA,rownames(res2)),
                   c("P value",ifelse(res2$p.value<0.0001,
                                      "< 0.0001",
                                      round(res2$p.value,4))),
                   c("Hazard Ratio(95% CI)",hz))
nrow(tabletext) + 1
library(forestplot)
pdf(file = '01.uncoxforest.pdf',height = 7, width = 12, onefile = F)
forestplot(labeltext=tabletext, 
           graph.pos=4,  #为Pvalue箱线图所在的位置
           is.summary = c(TRUE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,TRUE,FALSE,FALSE,FALSE,TRUE,rep(FALSE, 7)),
           col=fpColors(box="red", lines="royalblue", zero = "gray50"),
           mean=c(NA,res2$HR),
           lower=c(NA,res2$HR.95L), #95%置信区间下限
           upper=c(NA,res2$HR.95H), #95%置信区间上限
           boxsize=0.1,lwd.ci=3,   #箱子大小，线的宽度
           ci.vertices.height = 0.08,ci.vertices=TRUE, #置信区间用线宽、高、型
           zero=1,lwd.zero=0.5,    #zero线宽 基准线的位置
           colgap=unit(5,"mm"),    #列间隙
           xticks = c(seq(0,26,2)), #横坐标刻度
           lwd.xaxis=2,            #X轴线宽
           lineheight = unit(0.8,"cm"), #固定行高
           graphwidth = unit(.6,"npc"), #图在表中的宽度比例
           cex=0.9, fn.ci_norm = fpDrawCircleCI, #误差条显示方式
           hrzl_lines = list("2" = gpar(col = "black", lty = 1, lwd = 2)),
           txt_gp=fpTxtGp(label=gpar(cex=1),
                          ticks=gpar(cex=0.8, fontface = "bold"),
                          xlab=gpar(cex = 1, fontface = "bold"),
                          title=gpar(cex = 1.25, fontface = "bold")),
           xlab="Hazard Ratio",
           grid = T,
           title = "Univariate") # 垂直于x轴的网格线，对应每个刻度
dev.off()

png(file = '01.uncoxforest.png',height = 7, width = 12,units = 'in',res = 600)
forestplot(labeltext=tabletext, 
           graph.pos=4,  #为Pvalue箱线图所在的位置
           is.summary = c(TRUE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,TRUE,FALSE,FALSE,FALSE,TRUE,rep(FALSE, 7)),
           col=fpColors(box="red", lines="royalblue", zero = "gray50"),
           mean=c(NA,res2$HR),
           lower=c(NA,res2$HR.95L), #95%置信区间下限
           upper=c(NA,res2$HR.95H), #95%置信区间上限
           boxsize=0.1,lwd.ci=3,   #箱子大小，线的宽度
           ci.vertices.height = 0.08,ci.vertices=TRUE, #置信区间用线宽、高、型
           zero=1,lwd.zero=0.5,    #zero线宽 基准线的位置
           colgap=unit(5,"mm"),    #列间隙
           xticks = c(seq(0,26,2)), #横坐标刻度
           lwd.xaxis=2,            #X轴线宽
           lineheight = unit(0.8,"cm"), #固定行高
           graphwidth = unit(.6,"npc"), #图在表中的宽度比例
           cex=0.9, fn.ci_norm = fpDrawCircleCI, #误差条显示方式
           hrzl_lines = list("2" = gpar(col = "black", lty = 1, lwd = 2)),
           txt_gp=fpTxtGp(label=gpar(cex=1),
                          ticks=gpar(cex=0.8, fontface = "bold"),
                          xlab=gpar(cex = 1, fontface = "bold"),
                          title=gpar(cex = 1.25, fontface = "bold")),
           xlab="Hazard Ratio",
           grid = T,
           title = "Univariate") # 垂直于x轴的网格线，对应每个刻度
dev.off()

cox_more <- coxph(Surv(time = OS.time, event = OS) ~ MMP1+MMP3+Age+T.stage+N.stage+M.stage, data = train_risk_clinical)
cox_zph <- cox.zph(cox_more)
cox_table <- cox_zph$table[-nrow(cox_zph$table),]      #PH假定检验
write.csv(cox_table,file = 'PH.tab.csv',row.names = T)
res.mul = coxph(Surv(time = OS.time, event = OS) ~ MMP1+MMP3+Age+T.stage+N.stage+M.stage, data = train_risk_clinical)%>% summary
res.mul = cbind(res.mul$conf.int[,-2], res.mul$coefficients[,5]) %>% as.data.frame()
rownames(res.mul)

res.mul = rbind(res.mul[c(1:3),], res.ref, res.mul[c(4:6),],res.ref,res.mul[c(7:9),])

res.mul$Indicators = c("MMP1","MMP3",'Age','T1(Reference)','T2','T3','T4','N0(Reference)','N1','N2','M1 vs.M0')
colnames(res.mul) = c("hr","low","up","pv","Indicator")
res.mul$p = signif(res.mul$pv, 2) %>% paste0("p = ", .)
res.mul$p[is.na(res.mul$pv)] = NA
res.mul$Indicator = factor(res.mul$Indicator, levels = rev(res.mul$Indicator))
rownames(res.mul) <- res.mul$Indicator

multi_res <- data.frame(p.value=res.mul$pv,
                        HR=res.mul$hr,
                        HR.95L=res.mul$low,
                        HR.95H=res.mul$up,
                        Indicator=res.mul$Indicator)
rownames(multi_res) <- multi_res$Indicator
multi_res
write.csv(multi_res,
          file = "multivariate_cox_prog_result.csv",
          quote = F,
          row.names = T)
multi_res <- subset(multi_res, select = -c(Indicator))

library(tidyr)
hz <- paste(round(multi_res$HR,3),
            "(",round(multi_res$HR.95L,3),
            "-",round(multi_res$HR.95H,3),")",sep = "")
hz
hz[c(4,8)] <- ""
tabletext <- cbind(c(NA,rownames(multi_res)),
                   c("P value",ifelse(multi_res$p.value<0.0001,
                                      "< 0.0001",
                                      round(multi_res$p.value,4))),
                   c("Hazard Ratio(95% CI)",hz))
nrow(tabletext) + 1
pdf(file = '02.mulcoxforest.pdf',height = 5, width = 12, onefile = F)
forestplot(labeltext=tabletext, 
           graph.pos=4,  #为Pvalue箱线图所在的位置
           is.summary = c(TRUE,FALSE,FALSE,FALSE,TRUE,FALSE,FALSE,FALSE,TRUE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE),
           col=fpColors(box="red", lines="royalblue", zero = "gray50"),
           mean=c(NA,multi_res$HR),
           lower=c(NA,multi_res$HR.95L), #95%置信区间下限
           upper=c(NA,multi_res$HR.95H), #95%置信区间上限
           boxsize=0.1,lwd.ci=3,   #箱子大小，线的宽度
           ci.vertices.height = 0.08,ci.vertices=TRUE, #置信区间用线宽、高、型
           zero=1,lwd.zero=0.5,    #zero线宽 基准线的位置
           colgap=unit(5,"mm"),    #列间隙
           xticks = c(seq(0,24,2)), #横坐标刻度
           lwd.xaxis=2,            #X轴线宽
           lineheight = unit(0.8,"cm"), #固定行高
           graphwidth = unit(.6,"npc"), #图在表中的宽度比例
           cex=0.9, fn.ci_norm = fpDrawCircleCI, #误差条显示方式
           hrzl_lines = list("2" = gpar(col = "black", lty = 1, lwd = 2)),
           txt_gp=fpTxtGp(label=gpar(cex=1),
                          ticks=gpar(cex=0.8, fontface = "bold"),
                          xlab=gpar(cex = 1, fontface = "bold"),
                          title=gpar(cex = 1.25, fontface = "bold")),
           xlab="Hazard Ratio",
           grid = T,
           title = "Multivariate") # 垂直于x轴的网格线，对应每个刻度
dev.off()

png(file = '02.mulcoxforest.png',height = 5, width = 12,units = 'in',res = 600)
forestplot(labeltext=tabletext, 
           graph.pos=4,  #为Pvalue箱线图所在的位置
           is.summary = c(TRUE,FALSE,FALSE,FALSE,TRUE,FALSE,FALSE,FALSE,TRUE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE),
           col=fpColors(box="red", lines="royalblue", zero = "gray50"),
           mean=c(NA,multi_res$HR),
           lower=c(NA,multi_res$HR.95L), #95%置信区间下限
           upper=c(NA,multi_res$HR.95H), #95%置信区间上限
           boxsize=0.1,lwd.ci=3,   #箱子大小，线的宽度
           ci.vertices.height = 0.08,ci.vertices=TRUE, #置信区间用线宽、高、型
           zero=1,lwd.zero=0.5,    #zero线宽 基准线的位置
           colgap=unit(5,"mm"),    #列间隙
           xticks = c(seq(0,24,2)), #横坐标刻度
           lwd.xaxis=2,            #X轴线宽
           lineheight = unit(0.8,"cm"), #固定行高
           graphwidth = unit(.6,"npc"), #图在表中的宽度比例
           cex=0.9, fn.ci_norm = fpDrawCircleCI, #误差条显示方式
           hrzl_lines = list("2" = gpar(col = "black", lty = 1, lwd = 2)),
           txt_gp=fpTxtGp(label=gpar(cex=1),
                          ticks=gpar(cex=0.8, fontface = "bold"),
                          xlab=gpar(cex = 1, fontface = "bold"),
                          title=gpar(cex = 1.25, fontface = "bold")),
           xlab="Hazard Ratio",
           grid = T,
           title = "Multivariate") # 垂直于x轴的网格线，对应每个刻度
dev.off()
multi_cov<-c('MMP1','MMP3','Age','T.stage','N.stage','M.stage')
cox_data_prog <- as.formula(paste0('Surv(OS.time, OS)~',
                                   paste(multi_cov,
                                         sep = '',
                                         collapse = '+')))


cox_more_prog <- coxph(cox_data_prog,
                       data = as.data.frame(train_risk_clinical))

# Nomogram
library(rms)
ddist <- datadist(train_risk_clinical)
options(datadist='ddist')

# 构建COX模型，绘制列线图
train_risk_clinical$Age <- ifelse(train_risk_clinical$Age==1,'>65','<=65')
res.cox <- psm(cox_data_prog,
               data = train_risk_clinical, dist = 'lognormal')
surv <- Survival(res.cox) 
function(x) surv(365, x) 
function(x) surv(1095, x) 
function(x) surv(1825, x) 

nom.cox <- nomogram(res.cox,
                    fun = list(function(x) surv(365, x),
                               function(x) surv(1095, x),
                               function(x) surv(1825, x)),
                    funlabel=c("1-year Survival Probability", "3-year Survival Probability", "5-year Survival Probability"),
                    maxscale = 10,
                    fun.at = c(0.01,seq(0.1,0.9,by=0.2),0.95,0.99),
                    lp=F)
plot(nom.cox, cex.axis  = 1.5, cex.var = 1.7)
png(filename = "03.nomogram_line_points.png", height = 8, width = 17,units = 'in',res = 600)
plot(nom.cox, cex.axis  = 1.5, cex.var = 1.7)
dev.off()
pdf(file = "03.nomogram_line_points.pdf", height = 8, width = 17)
plot(nom.cox, cex.axis  = 1.5, cex.var = 1.7)
dev.off()

train_risk_clinical$Age <- ifelse(train_risk_clinical$Age == '>65',1,0)
coxm_1 <- cph(cox_data_prog,
              data=train_risk_clinical,
              surv=T,
              x=T,
              y=T,
              time.inc = 365*1)
cal_1<-calibrate(coxm_1,u=365*1,cmethod='KM',m=100,B=100)

coxm_3 <- cph(cox_data_prog,
              data=train_risk_clinical,
              surv=T,
              x=T,y=T,
              time.inc = 3*365)
cal_3 <-calibrate(coxm_3,u=3*365,cmethod='KM',m=100,B=100)

coxm_5 <- cph(cox_data_prog,
              data=train_risk_clinical,
              surv=T,
              x=T,y=T,
              time.inc = 5*365)
cal_5 <-calibrate(coxm_5,u=5*365,cmethod='KM',m=100,B=100)
pdf(file = '04.calibrate.pdf',w=8,h=7)
par(mar=c(7,4,4,3),cex=1.5)

plot(cal_1,
     subtitles = F,
     lwd=2,lty=1, ##设置线条形状和尺寸
     errbar.col=c(rgb(0,118,192,maxColorValue = 255)), ##设置一个颜色
     xlab='Nomogram-Predicted Probability of 1-5 year Survival Probability)',#便签
     ylab='Actual 1-5 year Survival Probability',#标签
     col="#519CBA",#设置一个颜色
     xlim = c(0,1),ylim = c(0,1)) ##x轴和y轴范围
plot(cal_3,
     add = T,
     subtitles = F,
     lwd=2,lty=1,  ##设置线条宽度和线条类型
     errbar.col=c(rgb(0,118,192,maxColorValue = 255)), ##设置一个颜色
     xlab='Nomogram-Predicted Probability of 1-5 year Survival Probability',#便签
     ylab='Actual 1-5 year Survival Probability',#标签
     col="#E05B5B",#设置一个颜色
     xlim = c(0,1),ylim = c(0,1)) ##x轴和y轴范围
plot(cal_5,
     add = T,
     subtitles = F,
     lwd=2,lty=1, ##设置线条形状和尺寸
     errbar.col=c(rgb(0,118,192,maxColorValue = 255)), ##设置一个颜色
     xlab='Nomogram-Predicted Probability of 1-5 year Survival Probability',#便签
     ylab='Actual 1-5 year Survival Probability',#标签
     col="#f7a27f",#设置一个颜色
     xlim = c(0,1),ylim = c(0,1)) ##x轴和y轴范围

#加上图例
legend("bottomright", legend=c("1-year", "3-year", "5-year"), 
       col=c("#519CBA", "#E05B5B", "#f7a27f"), 
       lwd=2)
#调整对角线
abline(0,1,lty=5,lwd=2,col="grey")
dev.off()

png(file = '04.calibrate.png',w=8,h=7,units = 'in',res = 600)
par(mar=c(7,4,4,3),cex=1.5)

plot(cal_1,
     subtitles = F,
     lwd=2,lty=1, ##设置线条形状和尺寸
     errbar.col=c(rgb(0,118,192,maxColorValue = 255)), ##设置一个颜色
     xlab='Nomogram-Predicted Probability of 1-5 year Survival Probability)',#便签
     ylab='Actual 1-5 year Survival Probability',#标签
     col="#519CBA",#设置一个颜色
     xlim = c(0,1),ylim = c(0,1)) ##x轴和y轴范围
plot(cal_3,
     add = T,
     subtitles = F,
     lwd=2,lty=1,  ##设置线条宽度和线条类型
     errbar.col=c(rgb(0,118,192,maxColorValue = 255)), ##设置一个颜色
     xlab='Nomogram-Predicted Probability of 1-5 year Survival Probability',#便签
     ylab='Actual 1-5 year Survival Probability',#标签
     col="#E05B5B",#设置一个颜色
     xlim = c(0,1),ylim = c(0,1)) ##x轴和y轴范围
plot(cal_5,
     add = T,
     subtitles = F,
     lwd=2,lty=1, ##设置线条形状和尺寸
     errbar.col=c(rgb(0,118,192,maxColorValue = 255)), ##设置一个颜色
     xlab='Nomogram-Predicted Probability of 1-5 year Survival Probability',#便签
     ylab='Actual 1-5 year Survival Probability',#标签
     col="#f7a27f",#设置一个颜色
     xlim = c(0,1),ylim = c(0,1)) ##x轴和y轴范围

#加上图例
legend("bottomright", legend=c("1-year", "3-year", "5-year"), 
       col=c("#519CBA", "#E05B5B", "#f7a27f"), 
       lwd=2)
#调整对角线
abline(0,1,lty=5,lwd=2,col="grey")
dev.off()


##roc
res.mul = coxph(Surv(time = OS.time, event = OS)~ MMP1+MMP3 +Age+T.stage+N.stage+M.stage, data = train_risk_clinical)
train_risk_clinical$riskScore_prog = predict(res.mul, newdata = train_risk_clinical, type = "lp")
train_risk_clinical$risk_prog = ifelse(train_risk_clinical$riskScore_prog > median(train_risk_clinical$riskScore_prog, na.rm = T), "low", "high")
train_risk_clinical$risk_prog = factor(train_risk_clinical$risk_prog, levels = c("high", "low"), labels = c("High risk", "Low risk"))
train_risk_clinical2 <- train_risk_clinical
library(survival)
library(survminer)

multi_ROC_out <- function(time_vector, risk_score_table){
  library(survivalROC)
  single_ROC <- function(single_time){
    for_ROC <- survivalROC(Stime=risk_score_table$OS.time,
                           status=risk_score_table$OS,
                           marker=risk_score_table$riskScore_prog,
                           predict.time=single_time,method = 'KM')
    data.frame('True_positive'=for_ROC$TP, 'False_positive'=for_ROC$FP,
               'Cut_values'=for_ROC$cut.values, 'Time'=rep(single_time, length(for_ROC$TP)),
               'AUC'=rep(for_ROC$AUC, length(for_ROC$TP)))
  }
  multi_ROC_list <- lapply(time_vector, single_ROC)
  do.call(rbind, multi_ROC_list)
}

for_multi_ROC <- multi_ROC_out(time_vector = c(365*seq(1,5,2)), 
                               risk_score_table = train_risk_clinical2)
for_multi_ROC$Time <- factor(for_multi_ROC$Time)

# 画ROC曲线 
library(scales)
library(geomROC)
library(plotROC)
auc_y1 <- round(for_multi_ROC[which(for_multi_ROC$Time==365),5][1],2)
auc_y3 <- round(for_multi_ROC[which(for_multi_ROC$Time==1095),5][1],2)
auc_y5 <- round(for_multi_ROC[which(for_multi_ROC$Time==1825),5][1],2)

ROC <- ggplot(for_multi_ROC, aes(x=False_positive,
                                 y=True_positive, 
                                 label=Cut_values, 
                                 color=Time)) + 
  scale_color_manual(breaks = c('365',"1095","1825"),
                     labels = c("1 years", "3 years","5 years"),
                     values = c("#4682B4", "#E05B5B", "#A28CC2")) +
  geom_roc(labels = F, stat = 'identity') + 
  style_roc() + 
  geom_abline(slope = 1, intercept = 0, color = 'gray', linetype=2) +
  theme_bw() +
  labs(title = "ROC curve") +
  theme(panel.background = element_rect(colour = 'black', size = 1, fill = 'white'),
        panel.grid = element_blank(),
        plot.title = element_text(size = 15, hjust = 0.5, face = "bold")) +
  annotate("text", x=0.75, y=c(0.25, 0.15, 0.05),
           label = c(paste('AUC of 1 years =', format(auc_y1,nsmall=2)),
                     paste('AUC of 3 years =', format(auc_y3,nsmall=2)),
                     paste('AUC of 5 years =', format(auc_y5,nsmall=2))))
ROC
ggsave('05.Prog_ROC.png', ROC,width = 5, height = 4,units = 'in',dpi = 600)
ggsave('05.Prog_ROC.pdf', ROC,width = 5, height = 4)

## DCA-----
library(rmda)
dca.dat <- train_risk_clinical
colnames(dca.dat)
complex<-decision_curve(OS ~MMP1+MMP3 +Age+T.stage+N.stage+M.stage,data =dca.dat,family = binomial(link ='logit'),
                        thresholds = seq(0,1, by = 0.01),
                        confidence.intervals= 0.95,
                        study.design = 'case-control',
                        population.prevalence= 0.3
)
MMP1.dca <- decision_curve(OS ~MMP1,data =dca.dat,family = binomial(link ='logit'),
                                thresholds = seq(0,1, by = 0.01),
                                confidence.intervals= 0.95,
                                study.design = 'case-control',
                                population.prevalence= 0.3
)
MMP3.dca <- decision_curve(OS ~MMP3,data =dca.dat,family = binomial(link ='logit'),
                           thresholds = seq(0,1, by = 0.01),
                           confidence.intervals= 0.95,
                           study.design = 'case-control',
                           population.prevalence= 0.3
)
age.dca <- decision_curve(OS ~Age,data =dca.dat,family = binomial(link ='logit'),
                             thresholds = seq(0,1, by = 0.01),
                             confidence.intervals= 0.95,
                             study.design = 'case-control',
                             population.prevalence= 0.3
)
tstage.dca <- decision_curve(OS ~T.stage,data =dca.dat,family = binomial(link ='logit'),
                             thresholds = seq(0,1, by = 0.01),
                             confidence.intervals= 0.95,
                             study.design = 'case-control',
                             population.prevalence= 0.3
)
nstage.dca <- decision_curve(OS ~N.stage,data =dca.dat,family = binomial(link ='logit'),
                             thresholds = seq(0,1, by = 0.01),
                             confidence.intervals= 0.95,
                             study.design = 'case-control',
                             population.prevalence= 0.3
)
mstage.dca <- decision_curve(OS ~M.stage,data =dca.dat,family = binomial(link ='logit'),
                             thresholds = seq(0,1, by = 0.01),
                             confidence.intervals= 0.95,
                             study.design = 'case-control',
                             population.prevalence= 0.3
)

dca.list <- list(MMP1.dca,MMP3.dca,age.dca,tstage.dca,nstage.dca,mstage.dca,complex)


pdf("06.DCA.pdf", width = 6, height = 6)
plot_decision_curve(dca.list,
                    curve.names=c('MMP1','MMP3','Age','T.stage','N.Stage','M.stage','Nomogram'),
                    confidence.intervals=FALSE,
                    standardize = FALSE)
dev.off()
png("06.DCA.png", width = 6, height = 6,units = 'in',res = 600)
plot_decision_curve(dca.list,
                    curve.names=c('MMP1','MMP3','Age','T.stage','N.Stage','M.stage','Nomogram'),
                    confidence.intervals=FALSE,
                    standardize = FALSE)
dev.off()
