# 清空环境变量
rm(list = ls())
# 自动设置工作目录（仅在 RStudio 中运行）
if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
}
# 日志文件路径
log_file <- "machine_log.txt"
# 开始记录日志
sink(log_file)
# 记录脚本开始时间
cat("Script started at:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")



library(randomForest)
library(ggplot2)
library(pheatmap)
library(pROC)
library(caret)
library(e1071)
library(party)
library(rpart)
# 训练模型
# 读取分组
design = read.table("model.txt", header = T, row.names = 1)
design$Type = as.factor(design$Type)
otu_table = read.table("matrix.txt", header = T, row.names = 1)
design_sub1 = subset(design, Group %in% c("group1"))
summary(design_sub1)


#判断"design_sub"的行名（细菌编号）是否在结构域矩阵"out_table"中
idx = rownames(design_sub1) %in% colnames(otu_table)
# 取出Group1对应的数据集
design_sub = design_sub1[idx,]
otu_sub = otu_table[, rownames(design_sub)]
summary(design_sub)
# 找到在 design_sub1 中存在但在 design_sub 中不存在的名称
missing_names <- rownames(design_sub1)[!rownames(design_sub1) %in% rownames(design_sub)]

# 打印结果
cat("在 design_sub1 中存在但在 design_sub 中不存在的名称：\n")
print(missing_names)

# 训练构建模型
set.seed(1001)
rf = randomForest(t(otu_sub), design_sub$Type, importance=TRUE, proximity=T, ntree = 2000)
print(rf)
# 保存模型
save(rf,file = 'model1_group3_1.RData')
# 直接加载模型
#load('Gram_model.RData')
# 交叉验证选择Features
set.seed(827) # 随机数据保证结果可重复，必须
# rfcv是随机森林交叉验证函数：Random Forest Cross Validation
result = rfcv(t(otu_sub), design_sub$Type, cv.fold=5)
save(result,file = 'result_rfcv1_group3_1.RData')

# 查看错误率表，68时错误率最低，为最佳模型
result$error.cv
error_data <- as.data.frame(result$error.cv)
write.table(error_data,file = 'error1_group3_1.txt',sep = '\t',row.names = T,
            quote = F,col.names = NA)
# 绘制验证结果 
with(result,plot(n.var, error.cv, log="x", type="o", lwd=2))# 交叉验证的结果建议多做5-6次，将结果统一在一张图上

# 导出训练集观测值与预测结果
train.p = predict(rf, type = "response")
df = data.frame(observed = design_sub$Type, predict = train.p)  #提取数据的子集作为另一部分数据再预测一下
# 保存预测结果与真实结果比较
write.table(df,file = "train_predict1_group3_1.txt",quote = F,sep = '\t', row.names = T, col.names = T)



# 导出feature重要性
imp= as.data.frame(rf$importance)
imp = imp[order(imp[,1],decreasing = T),]
head(imp,n=10)
write.table(imp,file = "importance_class.txt",quote = F,sep = '\t', row.names = T, col.names = T)
# 简单可视化
varImpPlot(rf, main = "Top 10 - Feature importance",n.var = 10, bg = par("bg"), color = par("fg"), gcolor = par("fg"), lcolor = "gray" )



# ggplot2美化feature贡献度
# 读取所有feature贡献度
imp = read.table("importance_class.txt", header=T, row.names= 1, sep="\t") 
# 分析选择top23分组效果最好
imp = head(imp, n=23)
# 反向排序X轴，让柱状图从上往下画
imp = imp[order(1:23,decreasing = T),]
#将imp的按第三列从小到大排序
imp = imp[order(imp[,3]),]
# 取出列名
imp$Domain = gsub("","",rownames(imp),perl=TRUE) 

imp$Domain=factor(imp$Domain,levels = imp$Domain)

# 图1. feature重要性柱状图
library(ggplot2)
p=ggplot(data = imp, mapping = aes(x=Domain,y=MeanDecreaseAccuracy,fill=Domain)) + 
  geom_bar(stat="identity")+coord_flip()+theme_bw()
ggsave(p,filename = "imp_shape1.pdf",width = 16,height = 9)


# group2验证
# design = read.table("group_Gr_all.txt",header = T, row.names = 1)
design_test = subset(design, Group %in% c("group2")) 
summary(design_test)
idx = rownames(design_test) %in% colnames(otu_table)
design_test = design_test[idx,]
otu_sub = otu_table[,rownames(design_test)]
summary(design_test)
# 转置，并添加分组信息
otutab_t = as.data.frame(t(otu_sub))
# 将Group2的分组信息添加到domain矩阵中
# 表示按照行名将矩阵"design"中"Group"列的内容添加到矩阵"otutab_t"中并将列名设置为"otutab_t$Group中的"Group",
otutab_t$Type = design[rownames(otutab_t),]$Type

set.seed(13)
otutab.pred = predict(rf, t(otu_sub) )  
pre_tab = table(observed=otutab_t[,"Type"],predicted=otutab.pred) 
pre_tab
# 整理样本原始分组和预测分类
predict = data.frame(Type = otutab_t[,"Type"], predicted=otutab.pred)
# 保存预测结果表
write.table("SampleID\t", file=paste("RF_prediction_binary.txt",sep=""),append = F, quote = F, eol = "", row.names = F, col.names = F)
write.table(predict, file = "RF_prediction_binary.txt",append = T, quote = F, row.names = T, col.names = T, sep = "\t")

# "========================================================================================================================"
# Load required libraries for new metrics
library(pROC)
library(caret)

# For training data evaluation
# Calculate confusion matrix and derived metrics
rf_train_conf_matrix <- confusionMatrix(data = as.factor(train.p), 
                                     reference = as.factor(design_sub$Type), 
                                     positive = "Spore_forming")
print(rf_train_conf_matrix)

# Save confusion matrix plot for training data
pdf("rf_train_confusion_matrix.pdf", width = 8, height = 6)
fourfoldplot(rf_train_conf_matrix$table, color = c("#e0f3f7", "#1a3b7a"),
             conf.level = 0, margin = 1, main = "Confusion Matrix - Training Data")
dev.off()

# Calculate ROC curve and AUC for training data
rf_train_roc <- roc(as.numeric(design_sub$Type), as.numeric(as.factor(train.p)))
rf_train_auc <- auc(rf_train_roc)

pdf("rf_train_roc_curve.pdf", width = 8, height = 8)
plot(rf_train_roc, main = paste("ROC Curve - Training Data\nAUC =", round(rf_train_auc, 3)), 
     col = "orange", lwd = 2,
     xlab = "False Positive Rate", ylab = "True Positive Rate")
abline(a = 0, b = 1, lty = 2, col = "gray")
dev.off()

# Save training metrics
rf_train_metrics <- data.frame(
  Dataset = "Training",
  Accuracy = rf_train_conf_matrix$overall["Accuracy"],
  Precision = rf_train_conf_matrix$byClass["Pos Pred Value"],
  Recall = rf_train_conf_matrix$byClass["Sensitivity"],
  F1_Score = rf_train_conf_matrix$byClass["F1"],
  AUC = as.numeric(rf_train_auc)
)

# For test data evaluation (group2)
# Calculate confusion matrix and derived metrics
rf_test_conf_matrix <- confusionMatrix(data = as.factor(otutab.pred), 
                                    reference = as.factor(otutab_t[,"Type"]), 
                                    positive = "Spore_forming")
print(rf_test_conf_matrix)


# Save confusion matrix plot for test data
pdf("rf_test_confusion_matrix.pdf", width = 8, height = 6)
fourfoldplot(rf_test_conf_matrix$table, color = c("#e0f3f7", "#1a3b7a"),
             conf.level = 0, margin = 1, main = "Confusion Matrix - Test Data")
dev.off()

# Calculate ROC curve and AUC for test data
rf_test_roc <- roc(as.numeric(as.factor(otutab_t[,"Type"])), 
                as.numeric(as.factor(otutab.pred)))
rf_test_auc <- auc(rf_test_roc)

pdf("rf_test_roc_curve.pdf", width = 8, height = 8)
plot(rf_test_roc, main = paste("ROC Curve - Test Data\nAUC =", round(rf_test_auc, 3)), 
     col = "orange", lwd = 2,
     xlab = "False Positive Rate", ylab = "True Positive Rate")
abline(a = 0, b = 1, lty = 2, col = "gray")
dev.off()

# Save test metrics
rf_test_metrics <- data.frame(
  Dataset = "Test",
  Accuracy = rf_test_conf_matrix$overall["Accuracy"],
  Precision = rf_test_conf_matrix$byClass["Pos Pred Value"],
  Recall = rf_test_conf_matrix$byClass["Sensitivity"],
  F1_Score = rf_test_conf_matrix$byClass["F1"],
  AUC = as.numeric(rf_test_auc)
)

# Combine and save all metrics
rf_all_metrics <- rbind(rf_test_metrics)   #rf_all_metrics <- rbind(rf_train_metrics, rf_test_metrics)
write.table(rf_all_metrics, file = "metrics_rf.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = TRUE)




library(e1071)
library(caret)
library(pROC)
# 读取矩阵
raw_data <- read.csv('matrix.txt', header = T, sep = '\t', row.names = 1, check.names = F)
raw_data <- t(raw_data)
dim(raw_data)
# 读取实验设计表
design <- read.csv('model.txt', header = T, sep = '\t', row.names = 1)
head(design)
# 以 group1 为训练集数据，group2 为测试集数据
design_group1 <- subset(design, Group %in% c('group1')) 
index_group1 <- rownames(design_group1) %in% rownames(raw_data)
table(index_group1)
# 确认分组后行数
cat("Group1 样本数：", nrow(design_group1), "\n")
train_data <- raw_data[rownames(design_group1), , drop = FALSE]
train_data <- as.data.frame(apply(train_data, 2, as.numeric))
# 构建测试集数据
design_group2 <- subset(design, Group %in% c('group2'))

# 检查 design_group2 中的行名是否在 raw_data 中存在
missing_rows <- rownames(design_group2)[!rownames(design_group2) %in% rownames(raw_data)]
print(missing_rows)



test_data <- raw_data[rownames(design_group2), , drop = FALSE]
test_data <- as.data.frame(apply(test_data, 2, as.numeric))

# 确保 Type 列为因子
design_group1$Type <- factor(design_group1$Type)
design_group2$Type <- factor(design_group2$Type)



##########################################################################################################
##############################################       SVM建模       #######################################
##########################################################################################################
set.seed(1213)
#####训练 SVM 模型
svm_model <- svm(design_group1$Type ~ ., train_data, probability = TRUE)

########## 训练集预测
svm_train_result <- predict(svm_model, train_data)
########## 测试集预测
svm_test_result <- predict(svm_model, test_data)




#####训练集预测混淆矩阵
svm_train_conf_matrix <- confusionMatrix(data = as.factor(svm_train_result), reference = as.factor(design_group1$Type), positive = "Spore_forming")
print(svm_train_conf_matrix)

#####测试集预测混淆矩阵
svm_test_conf_matrix <- confusionMatrix(data = as.factor(svm_test_result),  reference = as.factor(design_group2$Type), positive = "Spore_forming")
print(svm_test_conf_matrix)

########## 测试集混淆矩阵
# test_conf_matrix <- confusionMatrix(svm_test_result, design_group2$Type ,positive = "Spore_forming" )




######训练集混淆矩阵
pdf("svm_train_confusion_matrix.pdf", width = 8, height = 6)
fourfoldplot(svm_train_conf_matrix$table, color = c("#e0f3f7", "#1a3b7a"),
             conf.level = 0, margin = 1, main = "Confusion Matrix - Training Data")
dev.off()

##########训练集roc
svm_train_roc <- roc(as.numeric(design_group1$Type), as.numeric(as.factor(svm_train_result)))
svm_train_auc <- auc(svm_train_roc)

pdf("svm_train_roc_curve.pdf", width = 8, height = 8)
plot(svm_train_roc, main = paste("ROC Curve - Training Data\nAUC =", round(svm_train_auc, 3)), 
     col = "orange", lwd = 2,
     xlab = "False Positive Rate", ylab = "True Positive Rate")
abline(a = 0, b = 1, lty = 2, col = "gray")
dev.off()

#######训练集混淆矩阵
svm_train_metrics <- data.frame(
  Dataset = "Training",
  Accuracy = svm_train_conf_matrix$overall["Accuracy"],
  Precision = svm_train_conf_matrix$byClass["Pos Pred Value"],
  Recall = svm_train_conf_matrix$byClass["Sensitivity"],
  F1_Score = svm_train_conf_matrix$byClass["F1"],
  AUC = as.numeric(svm_train_auc)
)

pdf("svm_test_confusion_matrix.pdf", width = 8, height = 6)
fourfoldplot(svm_test_conf_matrix$table, color = c("#e0f3f7", "#1a3b7a"),
             conf.level = 0, margin = 1, main = "Confusion Matrix - Test Data")
dev.off()

########## 测试集roc
svm_test_roc <- roc(as.numeric(design_group2$Type), as.numeric(as.factor(svm_test_result)))
svm_test_auc <- auc(svm_test_roc)

pdf("svm_test_roc_curve.pdf", width = 8, height = 8)
plot(svm_test_roc, main = paste("ROC Curve - Test Data\nAUC =", round(svm_test_auc, 3)), 
     col = "orange", lwd = 2,
     xlab = "False Positive Rate", ylab = "True Positive Rate")
abline(a = 0, b = 1, lty = 2, col = "gray")
dev.off()

########## 测试集指标矩阵
svm_test_metrics <- data.frame(
  Dataset = "Test",
  Accuracy = svm_test_conf_matrix$overall["Accuracy"],
  Precision = svm_test_conf_matrix$byClass["Pos Pred Value"],
  Recall = svm_test_conf_matrix$byClass["Sensitivity"],
  F1_Score = svm_test_conf_matrix$byClass["F1"],
  AUC = as.numeric(svm_test_auc)
)

svm_all_metrics <- rbind(svm_test_metrics)   #rf_all_metrics <- rbind(rf_train_metrics, rf_test_metrics)
write.table(svm_all_metrics, file = "metrics_svm.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = TRUE)



##########################################################################################################
##############################################       NB建模       ########################################
##########################################################################################################
set.seed(1213)
# Naive Bayes（朴素贝叶斯算法分类）
nb_model <- naiveBayes(design_group1$Type ~ ., train_data)
nb_train_result <- predict(nb_model, train_data, type='class')
nb_test_result <- predict(nb_model, test_data, type = 'class')


#####训练集预测混淆矩阵
nb_train_conf_matrix <- confusionMatrix(data = as.factor(nb_train_result), reference = as.factor(design_group1$Type), positive = "Spore_forming")
print(nb_train_conf_matrix)

#####测试集预测混淆矩阵
nb_test_conf_matrix <- confusionMatrix(data = as.factor(nb_test_result),  reference = as.factor(design_group2$Type), positive = "Spore_forming")
print(nb_test_conf_matrix)

########## 测试集混淆矩阵
# test_conf_matrix <- confusionMatrix(nb_test_result, design_group2$Type ,positive = "Spore_forming" )


######训练集混淆矩阵
pdf("nb_train_confusion_matrix.pdf", width = 8, height = 6)
fourfoldplot(nb_train_conf_matrix$table, color = c("#e0f3f7", "#1a3b7a"),
             conf.level = 0, margin = 1, main = "Confusion Matrix - Training Data")
dev.off()

##########训练集roc
nb_train_roc <- roc(as.numeric(design_group1$Type), as.numeric(as.factor(nb_train_result)))
nb_train_auc <- auc(nb_train_roc)

pdf("nb_train_roc_curve.pdf", width = 8, height = 8)
plot(nb_train_roc, main = paste("ROC Curve - Training Data\nAUC =", round(nb_train_auc, 3)), 
     col = "orange", lwd = 2,
     xlab = "False Positive Rate", ylab = "True Positive Rate")
abline(a = 0, b = 1, lty = 2, col = "gray")
dev.off()

#######训练集混淆矩阵
nb_train_metrics <- data.frame(
  Dataset = "Training",
  Accuracy = nb_train_conf_matrix$overall["Accuracy"],
  Precision = nb_train_conf_matrix$byClass["Pos Pred Value"],
  Recall = nb_train_conf_matrix$byClass["Sensitivity"],
  F1_Score = nb_train_conf_matrix$byClass["F1"],
  AUC = as.numeric(nb_train_auc)
)

pdf("nb_test_confusion_matrix.pdf", width = 8, height = 6)
fourfoldplot(nb_test_conf_matrix$table, color = c("#e0f3f7", "#1a3b7a"),
             conf.level = 0, margin = 1, main = "Confusion Matrix - Test Data")
dev.off()

########## 测试集roc
nb_test_roc <- roc(as.numeric(design_group2$Type), as.numeric(as.factor(nb_test_result)))
nb_test_auc <- auc(nb_test_roc)

pdf("nb_test_roc_curve.pdf", width = 8, height = 8)
plot(nb_test_roc, main = paste("ROC Curve - Test Data\nAUC =", round(nb_test_auc, 3)), 
     col = "orange", lwd = 2,
     xlab = "False Positive Rate", ylab = "True Positive Rate")
abline(a = 0, b = 1, lty = 2, col = "gray")
dev.off()

########## 测试集指标矩阵
nb_test_metrics <- data.frame(
  Dataset = "Test",
  Accuracy = nb_test_conf_matrix$overall["Accuracy"],
  Precision = nb_test_conf_matrix$byClass["Pos Pred Value"],
  Recall = nb_test_conf_matrix$byClass["Sensitivity"],
  F1_Score = nb_test_conf_matrix$byClass["F1"],
  AUC = as.numeric(nb_test_auc)
)

nb_all_metrics <- rbind(nb_test_metrics)   #rf_all_metrics <- rbind(rf_train_metrics, rf_test_metrics)
write.table(nb_all_metrics, file = "metrics_nb.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = TRUE)



##########################################################################################################
##############################################       CTree建模       #####################################
##########################################################################################################
set.seed(1213)
ctree_model <- ctree(design_group1$Type ~ ., train_data)
ct_train_result <- predict(ctree_model, train_data)
ct_test_result <- predict(ctree_model, test_data)


#####训练集预测混淆矩阵
ct_train_conf_matrix <- confusionMatrix(data = as.factor(ct_train_result), reference = as.factor(design_group1$Type), positive = "Spore_forming")
print(ct_train_conf_matrix)

#####测试集预测混淆矩阵
ct_test_conf_matrix <- confusionMatrix(data = as.factor(ct_test_result),  reference = as.factor(design_group2$Type), positive = "Spore_forming")
print(ct_test_conf_matrix)

########## 测试集混淆矩阵
# test_conf_matrix <- confusionMatrix(ct_test_result, design_group2$Type ,positive = "Spore_forming" )




######训练集混淆矩阵
pdf("ct_train_confusion_matrix.pdf", width = 8, height = 6)
fourfoldplot(ct_train_conf_matrix$table, color = c("#e0f3f7", "#1a3b7a"),
             conf.level = 0, margin = 1, main = "Confusion Matrix - Training Data")
dev.off()

##########训练集roc
ct_train_roc <- roc(as.numeric(design_group1$Type), as.numeric(as.factor(ct_train_result)))
ct_train_auc <- auc(ct_train_roc)

pdf("ct_train_roc_curve.pdf", width = 8, height = 8)
plot(ct_train_roc, main = paste("ROC Curve - Training Data\nAUC =", round(ct_train_auc, 3)), 
     col = "orange", lwd = 2,
     xlab = "False Positive Rate", ylab = "True Positive Rate")
abline(a = 0, b = 1, lty = 2, col = "gray")
dev.off()

#######训练集混淆矩阵
ct_train_metrics <- data.frame(
  Dataset = "Training",
  Accuracy = ct_train_conf_matrix$overall["Accuracy"],
  Precision = ct_train_conf_matrix$byClass["Pos Pred Value"],
  Recall = ct_train_conf_matrix$byClass["Sensitivity"],
  F1_Score = ct_train_conf_matrix$byClass["F1"],
  AUC = as.numeric(ct_train_auc)
)

pdf("ct_test_confusion_matrix.pdf", width = 8, height = 6)
fourfoldplot(ct_test_conf_matrix$table, color = c("#e0f3f7", "#1a3b7a"),
             conf.level = 0, margin = 1, main = "Confusion Matrix - Test Data")
dev.off()

########## 测试集roc
ct_test_roc <- roc(as.numeric(design_group2$Type), as.numeric(as.factor(ct_test_result)))
ct_test_auc <- auc(ct_test_roc)

pdf("ct_test_roc_curve.pdf", width = 8, height = 8)
plot(ct_test_roc, main = paste("ROC Curve - Test Data\nAUC =", round(ct_test_auc, 3)), 
     col = "orange", lwd = 2,
     xlab = "False Positive Rate", ylab = "True Positive Rate")
abline(a = 0, b = 1, lty = 2, col = "gray")
dev.off()

########## 测试集指标矩阵
ct_test_metrics <- data.frame(
  Dataset = "Test",
  Accuracy = ct_test_conf_matrix$overall["Accuracy"],
  Precision = ct_test_conf_matrix$byClass["Pos Pred Value"],
  Recall = ct_test_conf_matrix$byClass["Sensitivity"],
  F1_Score = ct_test_conf_matrix$byClass["F1"],
  AUC = as.numeric(ct_test_auc)
)

ct_all_metrics <- rbind(ct_test_metrics)   #rf_all_metrics <- rbind(rf_train_metrics, rf_test_metrics)
write.table(ct_all_metrics, file = "metrics_ct.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = TRUE)




##########################################################################################################
##############################################      Decision_tree建模       ##############################
##########################################################################################################
set.seed(1213)
# 决策树训练模型
dtree_model <- rpart(design_group1$Type ~ ., train_data, method = 'class')
print(dtree_model$cptable)  # 打印复杂度参数表
# 决策树模型剪枝
dtree_model <- prune(dtree_model, cp = 0.01)
# 训练集预测
dt_train_result <- predict(dtree_model, train_data, type = 'class')
# 测试集预测
dt_test_result <- predict(dtree_model, test_data, type = 'class')




#####训练集预测混淆矩阵
dt_train_conf_matrix <- confusionMatrix(data = as.factor(dt_train_result), reference = as.factor(design_group1$Type), positive = "Spore_forming")
print(dt_train_conf_matrix)

#####测试集预测混淆矩阵
dt_test_conf_matrix <- confusionMatrix(data = as.factor(dt_test_result),  reference = as.factor(design_group2$Type), positive = "Spore_forming")
print(dt_test_conf_matrix)

########## 测试集混淆矩阵
# test_conf_matrix <- confusionMatrix(dt_test_result, design_group2$Type ,positive = "Spore_forming" )




######训练集混淆矩阵
pdf("dt_train_confusion_matrix.pdf", width = 8, height = 6)
fourfoldplot(dt_train_conf_matrix$table, color = c("#e0f3f7", "#1a3b7a"),
             conf.level = 0, margin = 1, main = "Confusion Matrix - Training Data")
dev.off()

##########训练集roc
dt_train_roc <- roc(as.numeric(design_group1$Type), as.numeric(as.factor(dt_train_result)))
dt_train_auc <- auc(dt_train_roc)

pdf("dt_train_roc_curve.pdf", width = 8, height = 8)
plot(dt_train_roc, main = paste("ROC Curve - Training Data\nAUC =", round(dt_train_auc, 3)), 
     col = "orange", lwd = 2,
     xlab = "False Positive Rate", ylab = "True Positive Rate")
abline(a = 0, b = 1, lty = 2, col = "gray")
dev.off()

#######训练集混淆矩阵
dt_train_metrics <- data.frame(
  Dataset = "Training",
  Accuracy = dt_train_conf_matrix$overall["Accuracy"],
  Precision = dt_train_conf_matrix$byClass["Pos Pred Value"],
  Recall = dt_train_conf_matrix$byClass["Sensitivity"],
  F1_Score = dt_train_conf_matrix$byClass["F1"],
  AUC = as.numeric(dt_train_auc)
)

pdf("dt_test_confusion_matrix.pdf", width = 8, height = 6)
fourfoldplot(dt_test_conf_matrix$table, color = c("#e0f3f7", "#1a3b7a"),
             conf.level = 0, margin = 1, main = "Confusion Matrix - Test Data")
dev.off()

########## 测试集roc
dt_test_roc <- roc(as.numeric(design_group2$Type), as.numeric(as.factor(dt_test_result)))
dt_test_auc <- auc(dt_test_roc)

pdf("dt_test_roc_curve.pdf", width = 8, height = 8)
plot(dt_test_roc, main = paste("ROC Curve - Test Data\nAUC =", round(dt_test_auc, 3)), 
     col = "orange", lwd = 2,
     xlab = "False Positive Rate", ylab = "True Positive Rate")
abline(a = 0, b = 1, lty = 2, col = "gray")
dev.off()

########## 测试集指标矩阵
dt_test_metrics <- data.frame(
  Dataset = "Test",
  Accuracy = dt_test_conf_matrix$overall["Accuracy"],
  Precision = dt_test_conf_matrix$byClass["Pos Pred Value"],
  Recall = dt_test_conf_matrix$byClass["Sensitivity"],
  F1_Score = dt_test_conf_matrix$byClass["F1"],
  AUC = as.numeric(dt_test_auc)
)

dt_all_metrics <- rbind(dt_test_metrics)   #rf_all_metrics <- rbind(rf_train_metrics, rf_test_metrics)
write.table(dt_all_metrics, file = "metrics_dt.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = TRUE)




###################################################################
#############################           合并测试集的ROC曲线
###################################################################
# 设置PDF输出
pdf("combined_test_roc_curves.pdf", width = 8, height = 8)

# 初始化绘图区域
plot(rf_test_roc, col = "red", 
     main = "Comparison of the Test_Set ROC Curves",
     lwd = 3,
     xlab = "False Positive Rate", 
     ylab = "True Positive Rate")

# 添加其他模型的ROC曲线
plot(svm_test_roc, add = TRUE, col = "royalblue", lwd = 3)
plot(nb_test_roc, add = TRUE, col = "forestgreen", lwd = 3)
plot(ct_test_roc, add = TRUE, col = "mediumorchid", lwd = 3)
plot(dt_test_roc, add = TRUE, col = "darkorange", lwd = 3)

# 添加对角线（参考线）
abline(a = 0, b = 1, lty = 2, col = "gray")

# 添加图例
legend("bottomright",
       legend = c(paste("Random_Forest (AUC =", sprintf("%.2f", rf_test_auc), ")"),
                  paste("SVM (AUC =", sprintf("%.2f", svm_test_auc), ")"),
                  paste("Naive_Bayes (AUC =", sprintf("%.2f", nb_test_auc), ")"),
                  paste("CTree (AUC =", sprintf("%.2f", ct_test_auc), ")"),
                  paste("Decision_tree (AUC =", sprintf("%.2f", dt_test_auc), ")")),
       col = c("red", "royalblue", "forestgreen", "mediumorchid", "darkorange"),
       lwd = 3, cex = 1.0)

# 关闭图形设备
dev.off()




library(ggplot2)
library(dplyr)
library(tidyr)

# 定义文件和模型名称
files <- c("metrics_rf.txt", "metrics_svm.txt", 
           "metrics_nb.txt", "metrics_ct.txt", 
           "metrics_dt.txt")
models <- c("Random_Forest", "SVM", "Naive_Bayes", "CTree", "Decision_tree")

# 读取并合并数据
combined_data <- data.frame()

for(i in seq_along(files)) {
  data <- read.table(files[i], header = TRUE, sep = "\t", quote = "\"", 
                     stringsAsFactors = FALSE)
  data$Model <- models[i]  # 添加模型名称列
  combined_data <- bind_rows(combined_data, data)
}

# 整理数据格式
combined_data <- combined_data %>%
  select(Dataset, Accuracy, Precision, Recall, F1_Score, AUC, Model) %>%
  pivot_longer(
    cols = c(Accuracy, Precision, Recall, F1_Score, AUC),
    names_to = "Metric",
    values_to = "Value"
  )

# 生成格式化标签（处理 NA 值）
combined_data_label <- combined_data %>%
  mutate(
    Label = case_when(
      is.na(Value) ~ "N/A",  # 将 NA 替换为 "N/A"
      Metric == "AUC" ~ sprintf("%.2f", Value),  # AUC 保留 2 位小数
      TRUE ~ sprintf("%.2f%%", Value * 100)     # 其他指标显示百分比
    )
  )

# 调整输出格式：将每个模型的指标放在一行
combined_data_wide <- combined_data_label %>%
  select(-Value) %>%  # 移除原始数值列
  pivot_wider(names_from = Metric, values_from = Label) %>%
  arrange(Model)  # 按模型名称排序

# 保存格式化文件
write.table(
  combined_data_wide,
  file = "combined_metrics.txt",
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

# 绘制热图（处理 NA 值）
p <- ggplot(combined_data_label, aes(x = Metric, y = Model, fill = Value)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = Label), color = "black", size = 5) +  # 使用格式化标签
  scale_fill_gradient2(
    low = "#2166AC",  # 浅蓝色
    mid = "white",  # 橙色
    high = "#B2182B",  # 深红色）
    midpoint = 0.5,   # 中间值
    limits = c(0, 1),  # 设置渐变范围为 0 到 1
    na.value = "gray"  # 将 NA 值的颜色设置为灰色
  ) +
  theme_minimal(base_size = 14) +
  labs(title = "Model Performance Metrics", x = NULL, y = NULL) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, color = "black", face = "bold"),
    axis.text.y = element_text(color = "black", face = "bold"),
    panel.grid = element_blank(),
    legend.position = "right"
  )

# 保存热图
ggsave("combined_metrics_heatmap.pdf", plot = p, width = 8, height = 6)


# 记录脚本结束时间
cat("\nScript finished at:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
# 关闭日志文件
sink()

