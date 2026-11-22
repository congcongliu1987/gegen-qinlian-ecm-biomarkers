rm(list = ls()) 
setwd('/home/data/project/30.JNPJ-6001-1/')
file.dir = '06_113model'
if (! dir.exists(file.dir)){
  dir.create(file.dir)
}
setwd(file.dir)

# 代码简介：12种具有变量筛选或预测模型构建的算法组合
# 
# 算法包括：Lasso、Ridge、Enet、Stepglm、SVM、glmBoost、LDA、plsRglm、RandomForest、GBM、XGBoost、NaiveBayes
# 
# 具体内容：
# 在交叉验证的框架下使用一种算法进行变量选择而使用另一种算法构建分类预测模型，
# 并计算所使用的模型组合(共113种)在外部数据集(或包括训练集)上的ROC曲线下面积(AUC)，最终通过热图将模型的评估结果可视化
# 
# 算法输入：
# 训练集表达谱：行为特征(如SYMBOL/ENSEMBL等，但需与测试集有一定交集)，列为样本的表达矩阵，格式见InputData文件中Train_expr.txt
# 训练集分类信息：行为样本，列包含至少一个需要预测的二分类变量(仅支持[0，1]格式，代码中需分别指定相对应的变量名)，格式见InputData文件中Train_class.txt
# 测试集表达谱：行为特征(如SYMBOL/ENSEMBL等，但需与训练集有一定交集)，列为样本的表达矩阵，格式见InputData文件中Test_expr.txt
# 测试集生存信息：行为样本，列包含至少一个需要预测的二分类变量(仅支持[0，1]格式，代码中需分别指定相对应的变量名)，以及一列用于指定队列信息的变量，格式见InputData文件中Test_class.txt
# 注意：建议对表达谱进行预筛选，保留几十至几千量级的感兴趣特征用于最优模型筛选，否则运行时间将大大延长
# 
# 本代码适用于：
# 存在大量对分类预测意义不明的变量时，需同时进行变量筛选，模型筛选和分类模型构建的场景
# 
# 代码输出结果包括：
# 113种算法在训练集上获得的模型(model.rds)，在测试集(或包括训练集)上的评估结果(AUC_mat.txt)，以及评估结果的热图展示(AUC.pdf)


# 12种具有变量筛选或预测模型构建的算法组合
library(openxlsx)
library(seqinr)
library(plyr)
library(randomForestSRC)
library(glmnet)
library(plsRglm)
library(gbm)
library(caret)
library(mboost)
library(e1071)
library(BART)
library(MASS)
library(snowfall)
library(xgboost)
library(ComplexHeatmap)
library(RColorBrewer)
library(pROC)

# 设置工作路径
work.path <- "../06_113model/"; 
# setwd(work.path) 

# 设置其他路径
code.path <- file.path(work.path, "Codes")
data.path <- file.path(work.path, "InputData")
res.path <- file.path(work.path, "Results")
fig.path <- file.path(work.path, "Figures")
if (!dir.exists(data.path)) dir.create(data.path)
if (!dir.exists(res.path)) dir.create(res.path)
if (!dir.exists(fig.path)) dir.create(fig.path)
if (!dir.exists(code.path)) dir.create(code.path)



# 加载模型训练以及模型评估的脚本
source(file.path(code.path, "ML.R"))
# 选择最后生成的模型类型：panML代表生成由不同算法构建的模型； multiLogistic表示抽取其他模型所用到的变量并建立多变量logistic模型
FinalModel <- c("panML", "multiLogistic")[2]

genes = read.csv("../04_candigene/candigene.csv")$symbol


## Training Cohort ---------------------------------------------------------
# 训练集表达谱是行为基因（感兴趣的基因集），列为样本的表达矩阵（基因名与测试集保持相同类型，如同为SYMBOL或ENSEMBL等）
# 行为样本，列包含至少一个需要预测的二分类变量(仅支持[0，1]格式)
# 提取训练集的共有样本
dat <- read.csv('../00_rawdata/dat.GSE87466.csv',row.names = 1)
group <- read.csv('../00_rawdata/group.GSE87466.csv')
table(group$group)
Train_expr <- as.data.frame(dat[genes,])
Train_class <- group
Train_class$group <- ifelse(Train_class$group == 'Normal',0,1)
colnames(Train_class) <- c('Sample','outcome')
Train_class$outcome <- as.numeric(Train_class$outcome)
comsam <- intersect(Train_class$Sample, colnames(Train_expr))
Train_expr <- Train_expr[,comsam]
Train_class <- Train_class[Train_class$Sample == comsam,,drop = F]

write.csv(Train_class,file = 'InputData/Train_class.csv',row.names = F,quote = F)
write.csv(Train_expr,file = 'InputData/Train_expr.csv',row.names = T,quote = F)


## Validation Cohort -------------------------------------------------------
# 测试集表达谱是行为基因（感兴趣的基因集），列为样本的表达矩阵（基因名与训练集保持相同类型，如同为SYMBOL或ENSEMBL等）
# 行为样本，列包含至少一个需要预测的二分类变量(仅支持[0，1]格式)，以及一列用于指定队列信息的变量
# 提取测试集的共有样本
dat.va <- read.csv('../00_rawdata/dat.GSE75214.csv',row.names = 1)
group.va <- read.csv('../00_rawdata/group.GSE75214.csv')
table(group.va$group)
Test_expr <- as.data.frame(dat.va[genes,])

Test_class <- data.frame(
  Sample = group.va$sample,
  outcome = group.va$group,
  Cohort = 'GSE75214'
)
Test_class$outcome <- ifelse(Test_class$outcome == 'Normal',0,1)
Test_class$outcome <- as.numeric(Test_class$outcome)
comsam <- intersect(Test_class$Sample, colnames(Test_expr))
Test_expr <- Test_expr[,comsam]
Test_class <- Test_class[Test_class$Sample == comsam,,drop = F]
write.csv(Test_class,file = 'InputData/Test_class.csv',row.names = F,quote = F)
write.csv(Test_expr,file = 'InputData/Test_expr.csv',row.names = T,quote = F)



# 提取相同基因
comgene <- intersect(rownames(Train_expr),rownames(Test_expr))
Train_expr <- t(Train_expr[comgene,]) # 输入模型的表达谱行为样本，列为基因
Test_expr <- t(Test_expr[comgene,]) # 输入模型的表达谱行为样本，列为基因


# 按队列对数据分别进行标准化（根据情况调整centerFlags和scaleFlags）
## data: 需要表达谱数据（行为样本，列为基因） 
## cohort：样本所属队列，为向量，不输入值时默认全表达矩阵来自同一队列
## centerFlag/scaleFlags：是否将基因均值/标准差标准化为1；
##        默认参数为NULL，表示不进行标准化；
##        为T/F时，表示对所有队列都进行/不进行标准化
##        输入由T/F组成的向量时，按顺序对队列进行处理，向量长度应与队列数一样
##        如centerFlags = c(F, F, F, T, T)，表示对第4、5个队列进行标准化，此时flag顺序应当与队列顺序一致
##        如centerFlags = c("A" = F, "C" = T, "B" = F)，表示对队列C进行标准化，此时不要求flag顺序与data一致

Train_set = scaleData(data = Train_expr, centerFlags = F, scaleFlags = F) # 数据的标准化，如果后续计算预测概率时遇到问题可以考虑替换centerFlags = F, scaleFlags = F
names(x = split(as.data.frame(Test_expr), f = Test_class$Cohort)) # 注意测试集标准化顺序与此一致
Test_set = scaleData(data = Test_expr, cohort = Test_class$Cohort, centerFlags = F, scaleFlags = F) # 数据的标准化，如果后续计算预测概率时遇到问题可以考虑替换centerFlags = F, scaleFlags = F

# Model training and validation -------------------------------------------

## method list --------------------------------------------------------
# 此处记录需要运行的模型，格式为：算法1名称[算法参数]+算法2名称[算法参数]
# 目前仅有Stepglm和Enet支持输入算法参数
methods <- read.xlsx(file.path(code.path, "methods.xlsx"), startRow = 2)
methods <- methods$Model
methods <- gsub("-| ", "", methods)

## Train the model --------------------------------------------------------

classVar = "outcome" # 设置所要预测的变量名（仅支持[0,1]二元变量格式）

classVar = "outcome" # 设置所要预测的变量名（仅支持[0,1]二元变量格式）
min.selected.var = 3 # 设置模型最少纳入的变量数
max.selected.var = 10 # 设置模型最多纳入的变量数

## Pre-training 将各方法所用到的变量筛选过程汇总，以减少计算量
Variable = colnames(Train_set)
preTrain.method =  strsplit(methods, "\\+") # 检视所有方法，分析各方法是否需要进行变量预筛选(pre-training)
preTrain.method = lapply(preTrain.method, function(x) rev(x)[-1]) # 删除各方法用于构建分类模型的算法，保留用于变量筛选的算法
preTrain.method = unique(unlist(preTrain.method)) # 汇总所有变量筛选算法，去除重复计算

preTrain.var <- list() # 用于保存各算法筛选的变量
set.seed(seed = 888) # 设置建模种子，使得结果可重复
for (method in preTrain.method){
  preTrain.var[[method]] = RunML(method = method, # 变量筛选所需要的机器学习方法
                                 Train_set = Train_set, # 训练集有潜在预测价值的变量
                                 Train_label = Train_class, # 训练集分类标签
                                 mode = "Variable",       # 运行模式，Variable(筛选变量)和Model(获取模型)
                                 classVar = classVar) # 用于训练的分类变量，必须出现在Train_class中
}
preTrain.var[["simple"]] <- colnames(Train_set)# 记录未经筛选的变量集（以便后续代码撰写），可视为使用simple方法（无筛选功能）的变量筛选结果

## Model training
model <- list() # 用于保存各模型的所有信息
set.seed(seed = 888) # 设置建模种子，使得结果可重复
Train_set_bk = Train_set # RunML有一个函数(plsRglm)无法正常传参，需要对训练集数据进行存档备份
for (method in methods){
  cat(match(method, methods), ":", method, "\n")
  method_name = method # 本轮算法名称
  method <- strsplit(method, "\\+")[[1]] # 各步骤算法名称
  
  if (length(method) == 1) method <- c("simple", method) # 如果本方法没有预筛选变量，则认为本方法使用simple方法进行了变量筛选
  
  Variable = preTrain.var[[method[1]]] # 根据方法名称的第一个值，调用先前变量筛选的结果
  Train_set = Train_set_bk[, Variable]   # 对训练集取子集，因为有一个算法原作者写的有点问题，无法正常传参
  Train_label = Train_class            # 所以此处需要修改变量名称，以免函数错误调用对象
  model[[method_name]] <- RunML(method = method[2],        # 根据方法名称第二个值，调用构建的函数分类模型
                                Train_set = Train_set,     # 训练集有潜在预测价值的变量
                                Train_label = Train_label, # 训练集分类标签
                                mode = "Model",            # 运行模式，Variable(筛选变量)和Model(获取模型)
                                classVar = classVar)       # 用于训练的分类变量，必须出现在Train_class中
  
  # 如果最终模型纳入的变量数小于预先设定的下限，则认为该算法输出的结果是无意义的
  if(length(ExtractVar(model[[method_name]])) <= min.selected.var) {
    model[[method_name]] <- NULL
  }
  # 如果最终模型纳入的变量数高于预先设定的上限，则认为该算法输出的结果是无意义的
  if(length(ExtractVar(model[[method_name]])) > max.selected.var) {
    model[[method_name]] <- NULL
  }
}
length(model) # 查看保留模型数量
##98
Train_set = Train_set_bk; rm(Train_set_bk) # 将数据还原，并移除备份
saveRDS(model, file.path(res.path, "model.rds")) # 报错各模型的所有中间过程

## Evaluate the model -----------------------------------------------------
# 读取已保存的模型列表
model <- readRDS(file.path(res.path, "model.rds"))
# model <- readRDS(file.path(res.path, "logisticmodel.rds")) # 若希望使用多变量保存最终以多变量Logistic模型计算得分，请运行此行
methodsValid <- names(model)

# 根据给定表达量计算样本风险评分
# 预测概率
RS_list <- list()
for (method in methodsValid){
  RS_list[[method]] <- CalPredictScore(fit = model[[method]], 
                                       new_data = rbind.data.frame(Train_set,Test_set)) # 2.0更新
}
RS_mat <- as.data.frame(t(do.call(rbind, RS_list)))
write.table(RS_mat, file.path(res.path, "RS_mat.txt"),sep = "\t", row.names = T, col.names = NA, quote = F) # 输出风险评分文件


# 根据给定表达量预测分类
Class_list <- list()
for (method in methodsValid){
  Class_list[[method]] <- PredictClass(fit = model[[method]], 
                                       new_data = rbind.data.frame(Train_set,Test_set)) # 2.0更新
}
Class_mat <- as.data.frame(t(do.call(rbind, Class_list)))
#Class_mat <- cbind.data.frame(Test_class, Class_mat[rownames(Class_mat),]) # 若要合并测试集本身的样本信息文件可运行此行
write.table(Class_mat, file.path(res.path, "Class_mat.txt"), # 测试集经过算法预测出的二分类结果
            sep = "\t", row.names = T, col.names = NA, quote = F)


# 提取所筛选的变量（列表格式）
fea_list <- list()
for (method in methodsValid) {
  fea_list[[method]] <- ExtractVar(model[[method]])
}

# 提取所筛选的变量（数据框格式）
fea_df <- lapply(model, function(fit){
  data.frame(ExtractVar(fit))
})
fea_df <- do.call(rbind, fea_df)
fea_df$algorithm <- gsub("(.+)\\.(.+$)", "\\1", rownames(fea_df))
colnames(fea_df)[1] <- "features"
write.table(fea_df, file.path(res.path, "fea_df.txt"), # 两列，包含算法以及算法所筛选出的变量
            sep = "\t", row.names = F, col.names = T, quote = F)



# 对各模型计算C-index

# 重新定义ROC曲线计算函数
RunEval = function(fit, 
                   Test_set = NULL, 
                   Test_label = NULL, 
                   Train_set = NULL, 
                   Train_label = NULL, 
                   Train_name = NULL,
                   cohortVar = "Cohort",
                   classVar){
  if(!is.element(cohortVar, colnames(Test_label))) {
    stop(paste0("There is no [", cohortVar, "] indicator, please fill in one more column!"))
  } 
  if((!is.null(Train_set)) & (!is.null(Train_label))) {
    new_data <- rbind.data.frame(Train_set[, fit$subFeature],
                                 Test_set[, fit$subFeature])
    if(!is.null(Train_name)) {
      Train_label$Cohort <- Train_name
    } else {
      Train_label$Cohort <- "Training"
    }
    colnames(Train_label)[ncol(Train_label)] <- cohortVar
    Test_label <- rbind.data.frame(Train_label[,c(cohortVar, classVar)],
                                   Test_label[,c(cohortVar, classVar)])
    Test_label[,1] <- factor(Test_label[,1], 
                             levels = c(unique(Train_label[,cohortVar]), setdiff(unique(Test_label[,cohortVar]),unique(Train_label[,cohortVar]))))
  } else {
    new_data <- Test_set[, fit$subFeature]
  }
  RS <- suppressWarnings(CalPredictScore(fit = fit, new_data = new_data))
  Predict.out <- Test_label
  Predict.out$RS <- as.vector(RS)
  # 计算Train的AUC
  roc_train = as.numeric(auc(suppressMessages(roc(Predict.out[1:dim(Train_label)[1],][[classVar]], Predict.out$RS[1:dim(Train_label)[1]]))))
  # 计算Test的AUC
  roc_test =  as.numeric(auc(suppressMessages(roc(Predict.out[(dim(Train_label)[1]+1):dim(Test_label)[1],][[classVar]], Predict.out$RS[(dim(Train_label)[1]+1):dim(Test_label)[1]]))))
  # 计算合并数据的AUC
  roc = as.numeric(auc(suppressMessages(roc(Predict.out[[classVar]], Predict.out$RS))))
  # 以表格形式输出
  x = data.frame(
    Train = roc_train,
    Test = roc_test
  )
  return(x)
}
AUC_list = list()
for (method in methodsValid){
  AUC_list[[method]] <- RunEval(fit = model[[method]],     # 分类预测模型
                                Test_set = Test_set,      # 测试集预测变量，应当包含训练集中所有的变量，否则会报错
                                Test_label = Test_class,   # 训练集分类数据，应当包含训练集中所有的变量，否则会报错
                                Train_set = Train_set,     # 若需要同时评估训练集，则给出训练集表达谱，否则置NULL
                                Train_label = Train_class, # 若需要同时评估训练集，则给出训练集分类数据，否则置NULL
                                Train_name = "Train",        # 若需要同时评估训练集，可给出训练集的标签，否则按“Training”处理
                                # cohortVar = cohortVar,      # 重要：用于指定队列的变量，该列必须存在且指定[默认为“Cohort”]，否则会报错
                                classVar = classVar)       # 用于评估的二元分类变量，必须出现在Test_class中
}
AUC_mat <- do.call(rbind, AUC_list)
write.table(AUC_mat, file.path(res.path, "AUC_mat.txt"),
            sep = "\t", row.names = T, col.names = T, quote = F)
# dim(AUC_mat)

# Plot --------------------------------------------------------------------

AUC_mat <- read.table(file.path(res.path, "AUC_mat.txt"),sep = "\t", row.names = 1, header = T,check.names = F,stringsAsFactors = F)
avg_AUC <- apply(AUC_mat, 1, mean)           # 计算每种算法在所有队列中平均AUC
avg_AUC <- sort(avg_AUC, decreasing = T)     # 对各算法AUC由高到低排序
AUC_mat <- AUC_mat[names(avg_AUC), ]      # 对AUC矩阵排序
fea_sel <- fea_list[[rownames(AUC_mat)[1]]]  # 最优模型（测试集AUC均值最大）所筛选的特征
avg_AUC <- as.numeric(format(avg_AUC, digits = 4, nsmall = 4)) # 保留三位小数


if(ncol(AUC_mat) < 3) { # 如果用于绘图的队列小于3个
  CohortCol <- c("red","blue") # 则给出两个颜色即可（可自行替换颜色）
} else { # 否则通过brewer.pal赋予超过3个队列的颜色
  CohortCol <- brewer.pal(n = ncol(AUC_mat), name = "Paired") # 设置队列颜色
}
names(CohortCol) <- colnames(AUC_mat)


cellwidth = 1; cellheight = 0.5
hm <- SimpleHeatmap(AUC_mat, # 主矩阵
                    avg_AUC, # 侧边柱状图
                    CohortCol, "steelblue", # 列标签颜色，右侧柱状图颜色
                    cellwidth = cellwidth, cellheight = cellheight, # 热图每个色块的尺寸
                    cluster_columns = F, cluster_rows = F) # 是否对行列进行聚
pdf(file.path(fig.path, "01.AUC.pdf"), width = cellwidth * ncol(AUC_mat) + 5, height = cellheight * nrow(AUC_mat) * 0.45)
draw(hm)
invisible(dev.off())
png(file.path(fig.path, "01.AUC.png"), width = cellwidth * ncol(AUC_mat) + 5, height = cellheight * nrow(AUC_mat) * 0.45,units = 'in',res = 600,bg = 'white')
draw(hm)
invisible(dev.off())

# 最优模型筛选
bestModel =rownames(AUC_mat)[1]
bestmodel = model[[bestModel]]
saveRDS(bestmodel,file = 'bestModel.rds')
##提取模型基因
modgene <- data.frame(symbol = bestmodel$subFeature)
write.csv(modgene,file = 'bestModel.gene.csv',row.names = F,quote = F)
# 
# # 定义ROC曲线的构建
RunROC = function(fit,
                  Test_set = NULL,
                  Test_label = NULL,
                  Train_set = NULL,
                  Train_label = NULL,
                  Train_name = NULL,
                  cohortVar = "Cohort",
                  classVar){
  if(!is.element(cohortVar, colnames(Test_label))) {
    stop(paste0("There is no [", cohortVar, "] indicator, please fill in one more column!"))
  }
  if((!is.null(Train_set)) & (!is.null(Train_label))) {
    new_data <- rbind.data.frame(Train_set[, fit$subFeature],
                                 Test_set[, fit$subFeature])
    if(!is.null(Train_name)) {
      Train_label$Cohort <- Train_name
    } else {
      Train_label$Cohort <- "Training"
    }
    colnames(Train_label)[ncol(Train_label)] <- cohortVar
    Test_label <- rbind.data.frame(Train_label[,c(cohortVar, classVar)],
                                   Test_label[,c(cohortVar, classVar)])
    Test_label[,1] <- factor(Test_label[,1],
                             levels = c(unique(Train_label[,cohortVar]), setdiff(unique(Test_label[,cohortVar]),unique(Train_label[,cohortVar]))))
  } else {
    new_data <- Test_set[, fit$subFeature]
  }
  RS <- suppressWarnings(CalPredictScore(fit = fit, new_data = new_data))
  Predict.out <- Test_label
  Predict.out$RS <- as.vector(RS)
  roc_train = as.numeric(auc(suppressMessages(roc(Predict.out[1:dim(Train_label)[1],][[classVar]], Predict.out$RS[1:dim(Train_label)[1]]))))
  roc_test =  as.numeric(auc(suppressMessages(roc(Predict.out[(dim(Train_label)[1]+1):dim(Test_label)[1],][[classVar]], Predict.out$RS[(dim(Train_label)[1]+1):dim(Test_label)[1]]))))
  roc = as.numeric(auc(suppressMessages(roc(Predict.out[[classVar]], Predict.out$RS))))
  Curve_Train = roc(Predict.out[1:dim(Train_label)[1],][[classVar]], Predict.out$RS[1:dim(Train_label)[1]])
  Curve_Test = roc(Predict.out[(dim(Train_label)[1]+1):dim(Test_label)[1],][[classVar]], Predict.out$RS[(dim(Train_label)[1]+1):dim(Test_label)[1]])
  ROC = list()
  ROC[["Train"]] = Curve_Train
  ROC[["Test"]] = Curve_Test
  return(ROC)
}

# 绘制ROC曲线
ROC = RunROC(fit = bestmodel,
             Test_set = Test_set,
             Test_label = Test_class,
             Train_set = Train_set,
             Train_label = Train_class,
             Train_name = "Train",
             # cohortVar = cohortVar,
             classVar = classVar)

# Plot ROC
Train_ROC = ROC[["Train"]]
pdf(file.path(fig.path, "02.BestModel_ROC_Train.pdf"),width = 4.5,height = 5)
plot(1-Train_ROC$specificities,Train_ROC$sensitivities,type ="l",
     col = "red",
     lty = 1,
     xlab = "1-Specifictity",
     ylab = "Sensitivity",
     lwd =2 )
abline(0,1)
legend("bottomright", legend=c(paste0("Corhort：GSE87466"),paste0("Model: ",bestModel),paste0("AUC=",round(as.numeric(Train_ROC$auc),4))))
dev.off()
png(file.path(fig.path, "02.BestModel_ROC_Train.png"),width = 4.5,height = 5,units = 'in',res = 600)
plot(1-Train_ROC$specificities,Train_ROC$sensitivities,type ="l",
     col = "red",
     lty = 1,
     xlab = "1-Specifictity",
     ylab = "Sensitivity",
     lwd =2 )
abline(0,1)
legend("bottomright", legend=c(paste0("Corhort：GSE87466"),paste0("Model: ",bestModel),paste0("AUC=",round(as.numeric(Train_ROC$auc),4))))
dev.off()


Test_ROC = ROC[["Test"]]

pdf(file.path(fig.path, "03.BestModel_ROC_Test.pdf"),width = 4.5,height = 5)
plot(1-Test_ROC$specificities,Test_ROC$sensitivities,type ="l",
     col = "red",
     lty = 1,
     xlab = "1-Specifictity",
     ylab = "Sensitivity",
     lwd =2 )
abline(0,1)
legend("bottomright", legend=c(paste0("Corhort：GSE75214"),paste0("Model: ",bestModel),paste0("AUC=",round(as.numeric(Test_ROC$auc),4))))
dev.off()

png(file.path(fig.path, "03.BestModel_ROC_Test.png"),width = 4.5,height = 5,units = 'in',res = 600)
plot(1-Test_ROC$specificities,Test_ROC$sensitivities,type ="l",
     col = "red",
     lty = 1,
     xlab = "1-Specifictity",
     ylab = "Sensitivity",
     lwd =2 )
abline(0,1)
legend("bottomright", legend=c(paste0("Corhort：GSE75214"),paste0("Model: ",bestModel),paste0("AUC=",round(as.numeric(Test_ROC$auc),4))))
dev.off()

modgene <- read.csv('bestModel.gene.csv')

##TCM-活性成分-关键基因⽹络图
all.net <- read.csv('../01_Target/allnet.csv',row.names = 1)
map <- all.net[which(all.net$symbol %in% modgene$symbol),]
map <- map[,c(2,3,5)]
map$link <- paste0(map$herbname,map$Molecule.Name.x,map$symbol)
map <- map[!duplicated(map$link),]
length(unique(map$herbname))
##4
length(unique(map$Molecule.Name.x))
##74
length(unique(map$symbol))
colnames(map)[2] <- 'Molecule.Name'
map <- map[,-4]

write.csv(map,file = 'maingenenet.csv',row.names = F,quote = T)


##ecm的功能网络
ecm.all <- read.csv('../01_Target/ecm.all.csv')
ecm.map <- ecm.all[which(ecm.all$Symbol %in% modgene$symbol),]
write.csv(ecm.map,file = 'mainECM.map.csv',row.names = F,quote = F)
table(ecm.map$Term)
##做个富集kegg-----

library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(enrichplot)
library(GOplot)
gene_transform <- bitr(unique(modgene$symbol),
                       fromType = "SYMBOL",
                       toType = c("ENTREZID",'ENSEMBL'),
                       OrgDb = "org.Hs.eg.db")
ego <- enrichGO(gene = gene_transform$ENTREZID,
                OrgDb = org.Hs.eg.db,
                keyType = "ENTREZID",
                ont = "ALL",
                pAdjustMethod = "BH",
                pvalueCutoff = 0.05,
                qvalueCutoff = 1,
                readable = TRUE)
write.csv(ego,file = "GO.csv",row.names = F)
# 展示富集最显著的 GO term
go_bar<-barplot(ego, showCategory=5, split="ONTOLOGY",label_format = 70) +
  facet_grid(ONTOLOGY ~ ., scales = "free",space = 'free')
go_bar
# ggsave('01.GO_bar.png',go_bar,width =8,height = 6,units = 'in',dpi = 600)
# ggsave('01.GO_bar.pdf',go_bar,width =8,height = 6)

## KEGG富集分析（气泡图）
kk <- enrichKEGG(gene = gene_transform$ENTREZID,
                 keyType = "kegg",
                 organism = "hsa",
                 pAdjustMethod = "BH",
                 pvalueCutoff = 0.05,
                 qvalueCutoff = 1)
kk <- setReadable(kk, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
# write.csv(kk@result,file = "KEGG.csv",row.names = F)
kk_dot <- dotplot(kk, showCategory=15,label_format = 50)
kk_dot
# ggsave('02.KEGG_dot.png',kk_dot,width =8,height = 6,units = 'in',dpi = 600)
# ggsave('02.KEGG_dot.pdf',kk_dot,width =8,height = 6)




