rm(list = ls())
setwd('/home/data/project/30.JNPJ-6001-1/')
file.dir = '17_progmodel_CRC_sup'
if (! dir.exists(file.dir)){
  dir.create(file.dir)
}
setwd(file.dir)
train_phenotype <- read.csv('../11_CRCsubtype/04_clinical/clinical_subtype.csv') 

train_phenotype$OS<-as.numeric(train_phenotype$OS)
train_phenotype$OS.time<-as.numeric(train_phenotype$OS.time)
train_phenotype2<-train_phenotype
colnames(train_phenotype2)
table(train_phenotype2$Stage)
train_phenotype2$Stage<-gsub('Stage ','',train_phenotype2$Stage)

table(train_phenotype2$Age)
train_phenotype2$Age<-ifelse(train_phenotype2$Age=='>65','1','0')

colnames(train_phenotype2)
dat<-read.csv('../00_rawdata/dat.fpkm.tcga.csv',row.names = 1,check.names = F)
dat <- log2(dat +1)
hubgene <- read.csv('../06_113model/bestModel.gene.csv')
sub_dat <- dat[hubgene$symbol,]
sub_dat <- t(sub_dat) %>% as.data.frame() %>% rownames_to_column(var = 'sample')
train_phenotype3 <- merge(sub_dat,
                          train_phenotype2,
                          by = "sample")
train_risk_clinical <- train_phenotype3
train_risk_clinical = subset(train_risk_clinical, select = -c(sample))
dim(train_risk_clinical)
colnames_train <- colnames(train_risk_clinical)
covariates_train <- colnames_train[-which(colnames_train %in% c("OS", "OS.time"))]
train_risk_clinical$T.stage <- factor(train_risk_clinical$T.stage)
train_risk_clinical$N.stage <- factor(train_risk_clinical$N.stage)

cox_more <- coxph(Surv(time = OS.time, event = OS) ~ MMP1+MMP3+Age+Stage, data = train_risk_clinical)
cox_zph <- cox.zph(cox_more)
cox_table <- cox_zph$table[-nrow(cox_zph$table),]     
cox_table
write.csv(cox_table,file = 'PH.tab.csv',row.names = T)
cox_summary <- summary(cox_more)
cindex <- cox_summary$concordance[1]
cindex_se <- cox_summary$concordance[2]
cindex
cindex_se
cindex <- data.frame(cindex = cindex,cindex_se = cindex_se)
write.csv(cindex,file = 'cindex.csv',row.names = F)
res.mul = coxph(Surv(time = OS.time, event = OS) ~ MMP1+MMP3+Age+Stage, data = train_risk_clinical)%>% summary
res.mul = cbind(res.mul$conf.int[,-2], res.mul$coefficients[,5]) %>% as.data.frame()
rownames(res.mul)
res.ref = c(1,1,1,NA)
res.mul = rbind(res.mul[c(1:3),], res.ref, res.mul[c(4:6),])

res.mul$Indicators = c("MMP1","MMP3",'Age','Stage1(Reference)','Stage2','Stage3','Stage4')
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
hz[c(4)] <- ""
tabletext <- cbind(c(NA,rownames(multi_res)),
                   c("P value",ifelse(multi_res$p.value<0.0001,
                                      "< 0.0001",
                                      round(multi_res$p.value,4))),
                   c("Hazard Ratio(95% CI)",hz))
nrow(tabletext) + 1
pdf(file = '01.mulcoxforest.pdf',height = 4, width = 12, onefile = F)
forestplot(labeltext=tabletext, 
           graph.pos=4,  
           is.summary = c(TRUE,FALSE,FALSE,FALSE,TRUE,FALSE,FALSE,FALSE,TRUE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE),
           col=fpColors(box="red", lines="royalblue", zero = "gray50"),
           mean=c(NA,multi_res$HR),
           lower=c(NA,multi_res$HR.95L),
           upper=c(NA,multi_res$HR.95H), 
           boxsize=0.1,lwd.ci=3,   
           ci.vertices.height = 0.08,ci.vertices=TRUE, 
           zero=1,lwd.zero=0.5, 
           colgap=unit(5,"mm"),   
           xticks = c(seq(0,24,2)),
           lwd.xaxis=2,            
           lineheight = unit(0.8,"cm"), 
           graphwidth = unit(.6,"npc"), 
           cex=0.9, fn.ci_norm = fpDrawCircleCI, 
           hrzl_lines = list("2" = gpar(col = "black", lty = 1, lwd = 2)),
           txt_gp=fpTxtGp(label=gpar(cex=1),
                          ticks=gpar(cex=0.8, fontface = "bold"),
                          xlab=gpar(cex = 1, fontface = "bold"),
                          title=gpar(cex = 1.25, fontface = "bold")),
           xlab="Hazard Ratio",
           grid = T,
           title = "Multivariate") 
dev.off()

png(file = '01.mulcoxforest.png',height = 4, width = 12,units = 'in',res = 600)
forestplot(labeltext=tabletext, 
           graph.pos=4, 
           is.summary = c(TRUE,FALSE,FALSE,FALSE,TRUE,FALSE,FALSE,FALSE,TRUE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE),
           col=fpColors(box="red", lines="royalblue", zero = "gray50"),
           mean=c(NA,multi_res$HR),
           lower=c(NA,multi_res$HR.95L),
           upper=c(NA,multi_res$HR.95H),
           boxsize=0.1,lwd.ci=3, 
           ci.vertices.height = 0.08,ci.vertices=TRUE, 
           zero=1,lwd.zero=0.5,   
           colgap=unit(5,"mm"),    
           xticks = c(seq(0,24,2)), 
           lwd.xaxis=2,           
           lineheight = unit(0.8,"cm"), 
           graphwidth = unit(.6,"npc"), 
           cex=0.9, fn.ci_norm = fpDrawCircleCI, 
           hrzl_lines = list("2" = gpar(col = "black", lty = 1, lwd = 2)),
           txt_gp=fpTxtGp(label=gpar(cex=1),
                          ticks=gpar(cex=0.8, fontface = "bold"),
                          xlab=gpar(cex = 1, fontface = "bold"),
                          title=gpar(cex = 1.25, fontface = "bold")),
           xlab="Hazard Ratio",
           grid = T,
           title = "Multivariate") 
dev.off()
multi_cov<-c('MMP1','MMP3','Age','Stage')
cox_data_prog <- as.formula(paste0('Surv(OS.time, OS)~',
                                   paste(multi_cov,
                                         sep = '',
                                         collapse = '+')))


cox_more_prog <- coxph(cox_data_prog,
                       data = as.data.frame(train_risk_clinical))

library(rms)
ddist <- datadist(train_risk_clinical)
options(datadist='ddist')
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
png(filename = "02.nomogram_line_points.png", height = 8, width = 17,units = 'in',res = 600)
plot(nom.cox, cex.axis  = 1.5, cex.var = 1.7)
dev.off()
pdf(file = "02.nomogram_line_points.pdf", height = 8, width = 17)
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
pdf(file = '03.calibrate.pdf',w=8,h=7)
par(mar=c(7,4,4,3),cex=1.5)

plot(cal_1,
     subtitles = F,
     lwd=2,lty=1, 
     errbar.col=c(rgb(0,118,192,maxColorValue = 255)), 
     xlab='Nomogram-Predicted Probability of 1-5 year Survival Probability)',
     ylab='Actual 1-5 year Survival Probability',
     col="#519CBA",
     xlim = c(0,1),ylim = c(0,1))
plot(cal_3,
     add = T,
     subtitles = F,
     lwd=2,lty=1, 
     errbar.col=c(rgb(0,118,192,maxColorValue = 255)),
     xlab='Nomogram-Predicted Probability of 1-5 year Survival Probability',
     ylab='Actual 1-5 year Survival Probability',
     col="#E05B5B",
     xlim = c(0,1),ylim = c(0,1)) 
plot(cal_5,
     add = T,
     subtitles = F,
     lwd=2,lty=1, 
     errbar.col=c(rgb(0,118,192,maxColorValue = 255)), 
     xlab='Nomogram-Predicted Probability of 1-5 year Survival Probability',
     ylab='Actual 1-5 year Survival Probability',
     col="#f7a27f",
     xlim = c(0,1),ylim = c(0,1)) 

legend("bottomright", legend=c("1-year", "3-year", "5-year"), 
       col=c("#519CBA", "#E05B5B", "#f7a27f"), 
       lwd=2)
abline(0,1,lty=5,lwd=2,col="grey")
dev.off()

png(file = '03.calibrate.png',w=8,h=7,units = 'in',res = 600)
par(mar=c(7,4,4,3),cex=1.5)

plot(cal_1,
     subtitles = F,
     lwd=2,lty=1,
     errbar.col=c(rgb(0,118,192,maxColorValue = 255)), 
     xlab='Nomogram-Predicted Probability of 1-5 year Survival Probability)',
     ylab='Actual 1-5 year Survival Probability',
     col="#519CBA",
     xlim = c(0,1),ylim = c(0,1))
plot(cal_3,
     add = T,
     subtitles = F,
     lwd=2,lty=1, 
     errbar.col=c(rgb(0,118,192,maxColorValue = 255)),
     xlab='Nomogram-Predicted Probability of 1-5 year Survival Probability',
     ylab='Actual 1-5 year Survival Probability',
     col="#E05B5B",
     xlim = c(0,1),ylim = c(0,1)) 
plot(cal_5,
     add = T,
     subtitles = F,
     lwd=2,lty=1, 
     errbar.col=c(rgb(0,118,192,maxColorValue = 255)),
     xlab='Nomogram-Predicted Probability of 1-5 year Survival Probability',
     ylab='Actual 1-5 year Survival Probability',
     col="#f7a27f",
     xlim = c(0,1),ylim = c(0,1)) 


legend("bottomright", legend=c("1-year", "3-year", "5-year"), 
       col=c("#519CBA", "#E05B5B", "#f7a27f"), 
       lwd=2)

abline(0,1,lty=5,lwd=2,col="grey")
dev.off()

res.mul = coxph(Surv(time = OS.time, event = OS)~ MMP1+MMP3 +Age+Stage, data = train_risk_clinical)
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
ggsave('04.Prog_ROC.png', ROC,width = 5, height = 4,units = 'in',dpi = 600)
ggsave('04.Prog_ROC.pdf', ROC,width = 5, height = 4)

library(rmda)
dca.dat <- train_risk_clinical
colnames(dca.dat)
complex<-decision_curve(OS ~MMP1+MMP3 +Age+Stage,data =dca.dat,family = binomial(link ='logit'),
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
stage.dca <- decision_curve(OS ~Stage,data =dca.dat,family = binomial(link ='logit'),
                             thresholds = seq(0,1, by = 0.01),
                             confidence.intervals= 0.95,
                             study.design = 'case-control',
                             population.prevalence= 0.3
)

dca.list <- list(MMP1.dca,MMP3.dca,age.dca,stage.dca,complex)


pdf("05.DCA.pdf", width = 6, height = 6)
plot_decision_curve(dca.list,
                    curve.names=c('MMP1','MMP3','Age','Stage','Nomogram'),
                    confidence.intervals=FALSE,
                    standardize = FALSE)
dev.off()
png("05.DCA.png", width = 6, height = 6,units = 'in',res = 600)
plot_decision_curve(dca.list,
                    curve.names=c('MMP1','MMP3','Age','Stage','Nomogram'),
                    confidence.intervals=FALSE,
                    standardize = FALSE)
dev.off()
