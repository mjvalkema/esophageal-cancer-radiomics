# Model extension / updating
# Authors: M.J. Valkema, L.R. de Ruiter, January 2022, based on scripts created by A. Chatterjee.

################# Initialization ###################
setwd(dir=getwd())
source('DataPrepModUpdating.R')# allData - this is the standardized data you use in this scrpt, available from DataPrepModUpdating.R
library(caret)
library(rms)
library(ROCR)
library(pROC)
library(ggplot2)
library(dplyr)
library(CalibrationCurves)
source('FeatureSelection.R')
source('ROCPlot.R')
library(ggplot2)
source('SigmoidFunction.R')
colors <- c("#4E79A7","#F28E2B", "#E15759")

# Organise data and split into different types of data
# For modeling
outcomes <- allData[, "outcome"] # TRG 1-2 is 0, TRG 3-4 is 1
features <- allData[, 2:102] # all the normalized features without PatientID
features <- features[, order(names(features))] # order alphabetically, so feature selection in correlation matrix is not dependent on order of features
# For modeling
clinical <- allData[, c('cT', 'cN', 'Gender', 'Age', 'Histology')]

# N.B. The dataset  is split in a Training set (for developing) and Validation set (for internal validation). 
# N.B. The nomenclature in the functions is Train:Test. For the current dataset, this means Train:Validate since an independent Test set was not available.

# Feature selection step
nTrainTestSplits <- 100 
selectedFeatures <- RadiomicsFeatureSelection(features, outcomes)

# Make train-test splits for model building
trainTestSplits <- CreateTrainTestSplits(outcomes, nTrainTestSplits = 100, trainTestFraction = 0.67)

# Rank selected features based on AUC, to see which features get highest AUCs
arrayAUC <- array(NA, dim=c(length(selectedFeatures), nTrainTestSplits))

for (n.split in 1:nTrainTestSplits) {
  split <- trainTestSplits[, n.split]
  for (n.feature in 1:ncol(selectedFeatures)) { #loopen over aantal features
    feature <- selectedFeatures[, n.feature] # select column number to get a feature
    AUC <- auc(response = outcomes[split], predictor = feature[split], quiet=T)
    corr <- cor(x=feature[split], y=outcomes[split], use = "everything", method = "spearman")
    if (AUC < 0.5 & corr < 0) {
      AUC <- (1 - AUC)
    }
    arrayAUC[n.feature, n.split] <- AUC
  }
}

medians <- apply(arrayAUC, MARGIN=1, FUN=median)
ranking <- as.data.frame(cbind(names(selectedFeatures), medians))
ranking <- ranking[order(-medians),]
namesRankedFeatures <- ranking$V1[1:7]
ranking$medians <- round(as.numeric(ranking$medians), digits =2)
write.csv(ranking, "output_mu/names_selectedFeatures.csv")

# Create the model data based on selected features
modelData <- cbind(clinical, selectedFeatures, outcomes)
modelFeatures <- cbind(clinical, selectedFeatures)
modelOutcomes <- as.numeric(as.character(outcomes))
randOutcome <- round(runif(length(modelOutcomes))) 

##### Feature selection with LASSO ###########
############# Multivariate variable selection with LASSO ############
set.seed(1234)
library(glmnet)
rm(coefs)
for (split in trainTestSplits) {
  cv <- cv.glmnet(data.matrix(modelFeatures[split,]), modelOutcomes[split], alpha = 1)
  fit <- glmnet(data.matrix(modelFeatures)[split,], modelOutcomes[split], alpha = 1, lambda = cv$lambda.min)
  print(cv$lambda.min)
  temp <- as.data.frame(as.matrix(coef(fit)))
  temp$feature <- rownames(temp)
  
  if (!exists('coefs')) {
    coefs <- temp
  } else {
    coefs <- merge(coefs, temp, by="feature", all.x = TRUE) # Don't worry about the warnings!
  }
}
coefs[coefs == 0] <- NA
rownames(coefs) <- coefs$feature
coefs$feature <- NULL
colnames(coefs) <- 1:ncol(trainTestSplits)

results <- data.frame(count=rowSums(!is.na(coefs)))
rownames(results) <- rownames(coefs)

results$mean <- apply(coefs, MARGIN=1, FUN=mean, na.rm = T)
results$min <- apply(coefs, MARGIN=1, FUN=min, na.rm = T)
results$max <- apply(coefs, MARGIN=1, FUN=max, na.rm = T)
results$median <- apply(coefs, MARGIN=1, FUN=median, na.rm = T)
results$lower <- apply(coefs, MARGIN=1, FUN=quantile, probs = 0.025, na.rm = T)
results$upper <- apply(coefs, MARGIN=1, FUN=quantile, probs = 0.975, na.rm = T)

# Show results of the LASSO coefficients
head(results[order(-results$count),])
LASSOresults <- data.frame(results[order(-results$count),])
write.csv(LASSOresults, "output_mu/LASSO_results.csv")
finalFeatureNames <- rownames(head(results[order(-results$count),][-1,])) 

#Boxplots to take a glance at the normalized features for the TRG outcomes
outcomeBoxplot <- car::recode(modelData$outcomes, "c(0)='TRG 1-2'; c(1)='TRG 3-4'")
outcomeBoxplot <- as.factor(outcomeBoxplot)
columnsOfInterestLabels <- c("flatness", 
                             "minimum intensity", 
                             "Gearys C measure")

BoxPlotColumn <- function (column, data, y_label) {
  ggplot(data, aes_string(x=outcomeBoxplot, y=column))+
    geom_boxplot() +
    labs(y=y_label) +
    theme(legend.position="none") +
    xlab(NULL) +
    theme(axis.text.x = element_text(angle = 0, vjust=0, size = 12)) +
    theme(axis.title.y = element_text(size = 12)) +
    geom_signif(comparisons = list(c("TRG 1-2", "TRG 3-4")), 
                test = "wilcox.test", map_signif_level=TRUE, vjust = 0.4, tip_length = 0)
}

png(filename = "figures_mu/TRGversusOutcome_LR_LASSO.png", units = "cm", width=17, height=10, res=300)
figures <- list()
for (i in 1:length(finalFeatureNames[c(1,3,4)])) {
  figures[[i]] <- BoxPlotColumn(finalFeatureNames[c(1, 3, 4)][i], data = modelData, columnsOfInterestLabels[i])
}
figure <- ggarrange(plotlist = figures,
                    labels = c("A", "B", "C"),
                    ncol =3, nrow =1,
                    font.label = list(size = 12),
                    common.legend= TRUE, 
                    legend = "bottom")
figure
dev.off()

################# exploreSimpleModels ###############
# Make 3 simple models: LASSO in combination with logistic regression, SVM, Naive Bayes
metrics <- c("threshold", "youden", "sensitivity", "specificity", "ppv", "npv", "accuracy")

# Logistic regression model - with 4 variables as chosen from LASSO
featureNames <- finalFeatureNames[1:4] # number of variables in final model, not too many to keep model simple (Avishek paper An Empirical Approach)
metricsArrayTrain <- array(NA, 
                           dim = c(100, length(metrics)), 
                           dimnames = list('split' = 1:100, 'metric' = metrics))
metricsArrayTest <- array(NA, 
                          dim = c(100, length(metrics)), 
                          dimnames = list('split' = 1:100, 'metric' = metrics))
AUCTrain <- array(NA,
                  dim = 100,
                  dimnames = list('split' = 1:100))

AUCTest <- array(NA,
                 dim = 100,
                 dimnames = list('split' = 1:100))
                 
coefficients_array <- array(NA, dim = c(100, (length(featureNames)+1)), dimnames = list('split' = 1:100, 'feature' = c('intercept', featureNames)))


for (i in 1:ncol(trainTestSplits)) {
  # Select data
  split <- trainTestSplits[, i]
  modelTrainData <- modelFeatures[split,]
  modelTestData <- modelFeatures[!split,]
  modelOutcomesTrain <- as.factor(modelOutcomes[split])
  modelOutcomesTest <- as.factor(modelOutcomes[!split])
  
  # Apply model to training data, threshold is set with Youden
  formula <- as.formula(paste('modelOutcomesTrain ~', paste(featureNames, collapse = ' + ')))
  model <- glm(formula, data = modelTrainData, family = "binomial")
  predTrain <- predict.glm(model, type = "response", newdata = NULL)
  roc_object<- roc(modelOutcomesTrain, as.numeric(as.character(predTrain)))
  metricsArrayTrain[i, ] <- coords(roc_object, x="best", best.method = "youden", transpose=FALSE, as.matrix = TRUE, ret = metrics)[1, ]
  AUCTrain[i] <- auc(roc_object)
  
  # output Test for every split
  predTest <- predict.glm(model, newdata = modelTestData, type = "response")
  roc_object <- roc(modelOutcomesTest, as.numeric(as.character(predTest)))
  metricsArrayTest[i, ] <- coords(roc_object, x="best", best.method = "youden", transpose=FALSE, as.matrix = TRUE, ret = metrics)[1, ]
  AUCTest[i] <- auc(roc_object)
  coefficients_array[i, ] <- coef(model)

}

MetricsSummary <- function(obj) {
  # Writes summary in matrix form for summary of metricsArrayTrainobjects with numbers instead of strings
  out = array(as.numeric(sapply(strsplit(summary(obj), ":"), "[[", 2)), 
          dim=dim(summary(obj)),
          dimnames=list(row = c('min', '25%', 'median', 'mean', '75%', 'max'),
                        metric = dimnames(obj)$metric))
  return(out)
}

write.csv(MetricsSummary(metricsArrayTrain), "output_mu/metrics_train_LR.csv")
write.csv(MetricsSummary(metricsArrayTest), "output_mu/metrics_test_LR.csv")
summary(AUCTrain)
summary(AUCTest)
coefficients_array
boxplot(coefficients_array)

dfSummary <- as.data.frame(MetricsSummary(metricsArrayTrain))
dfSummary$AUC <- c(summary(AUCTrain))
dfSummary <- as.data.frame(t(dfSummary))
dfSummary
dfSummary$model <- c("Logistic regression")
dfSummary$metric <- rownames(dfSummary)
dfSummary$type <- c("Training")
dfSummary

dfSummary2 <- as.data.frame(MetricsSummary(metricsArrayTest))
dfSummary2$AUC <- c(summary(AUCTest))
dfSummary2 <- as.data.frame(t(dfSummary2))
dfSummary2
dfSummary2$model <- c("Logistic regression")
dfSummary2$metric <- rownames(dfSummary2)
dfSummary2$type <- c("Validation")
dfSummary2

dfSummary3 <- rbind(dfSummary, dfSummary2)
dfSummary3

png(filename = "figures_mu/LR_medianvalues_100.png", units = "cm", width=10, height=10, res=300)
ggplot(data=dfSummary3[dfSummary3$metric %in% c("sensitivity", "specificity", "accuracy", "AUC"),], 
       aes(x = metric, y=median, fill=type)) + 
  geom_bar(stat="identity", position=position_dodge(0.6), width = 0.5) +
  ylim(0, 1) +
  scale_fill_manual("type", values = colors) +
  theme(legend.title = element_blank(), axis.title.x = element_blank(), plot.title = element_text(hjust = 0.5), axis.text = element_text(size=12)) +
  ggtitle("Logistic regression") +
  labs(y="median value") +
  theme(legend.position = c(0.98, 0.98), legend.justification = c("right", "top"), legend.box.just = "right", legend.margin = margin(2,4,4,4))
dev.off()

# Logistic regression model with outcome random - with 4 variables as chosen from LASSO
# - see if you are able to get a model whose AUC is better than random at 95% CI. If you cannot, it means the workflow is robust against false discovery. 
featureNames <- finalFeatureNames[1:4] # number of variables in final model, not too many to keep model simple (Avishek paper An Empirical Approach)
metricsArrayTrain <- array(NA, 
                           dim = c(100, length(metrics)), 
                           dimnames = list('split' = 1:100, 'metric' = metrics))
metricsArrayTest <- array(NA, 
                          dim = c(100, length(metrics)), 
                          dimnames = list('split' = 1:100, 'metric' = metrics))
AUCTrain <- array(NA,
                  dim = 100,
                  dimnames = list('split' = 1:100))

AUCTest <- array(NA,
                 dim = 100,
                 dimnames = list('split' = 1:100))

coefficients_array <- array(NA, dim = c(100, (length(featureNames)+1)), dimnames = list('split' = 1:100, 'feature' = c('intercept', featureNames)))

for (i in 1:ncol(trainTestSplits)) {
  # Select data
  split <- trainTestSplits[, i]
  modelTrainData <- modelFeatures[split,]
  modelTestData <- modelFeatures[!split,]
  modelOutcomesTrain <- as.factor(randOutcome[split])
  modelOutcomesTest <- as.factor(randOutcome[!split])
  
  # Apply model to training data, threshold is set at 0.5
  formula <- as.formula(paste('modelOutcomesTrain ~', paste(featureNames, collapse = ' + ')))
  model <- glm(formula, data = modelTrainData, family = "binomial")
  predTrain <- predict.glm(model, type = "response", newdata = NULL)
  roc_object<- roc(modelOutcomesTrain, as.numeric(as.character(predTrain)))
  metricsArrayTrain[i, ] <- coords(roc_object, x="best", best.method = "youden", transpose=FALSE, as.matrix = TRUE, ret = metrics)[1, ]
  AUCTrain[i] <- auc(roc_object)
  
  # output Test for every split
  predTest <- predict.glm(model, newdata = modelTestData, type = "response")
  roc_object <- roc(modelOutcomesTest, as.numeric(as.character(predTest)))
  metricsArrayTest[i, ] <- coords(roc_object, x="best", best.method = "youden", transpose=FALSE, as.matrix = TRUE, ret = metrics)[1, ]
  AUCTest[i] <- auc(roc_object)
  coefficients_array[i, ] <- coef(model)
  
}

write.csv(MetricsSummary(metricsArrayTrain), "output_mu/metrics_train_LR_randomOutcome.csv")
write.csv(MetricsSummary(metricsArrayTest), "output_mu/metrics_test_LR_randomOutcome.csv")

summary(AUCTrain)
summary(AUCTest)
coefficients_array
boxplot(coefficients_array)

######### svmModel #######
library(e1071)
featureNames <- finalFeatureNames[1:4] # number of variables in final model, not too many to keep model simple (Avishek paper An Empirical Approach)
metrics <- c("threshold", "youden", "sensitivity", "specificity", "ppv", "npv", "accuracy")
metricsArrayTrain <- array(NA, 
                           dim = c(100, length(metrics)), 
                           dimnames = list('split' = 1:100, 'metric' = metrics))
metricsArrayTest <- array(NA, 
                          dim = c(100, length(metrics)), 
                          dimnames = list('split' = 1:100, 'metric' = metrics))
AUCTrain <- array(NA,
                  dim = 100,
                  dimnames = list('split' = 1:100))

AUCTest <- array(NA,
                 dim = 100,
                 dimnames = list('split' = 1:100))


for (i in 1:ncol(trainTestSplits)) {
  # Select data
  split <- trainTestSplits[, i]
  modelTrainData <- modelFeatures[split,]
  modelTestData <- modelFeatures[!split,]
  modelOutcomesTrain <- as.factor(modelOutcomes[split])
  modelOutcomesTest <- as.factor(modelOutcomes[!split])
  
  # Apply model to training data, threshold is set at 0.5
  formula <- as.formula(paste('modelOutcomesTrain ~', paste(featureNames, collapse = ' + ')))
  model <- svm(formula, data = modelTrainData, kernel = "linear", cost = 1, scale = FALSE)
  predTrain <- predict(model)
  roc_object<- roc(modelOutcomesTrain, as.numeric(as.character(predTrain)))
  metricsArrayTrain[i, ] <- coords(roc_object, x=0.5, transpose=FALSE, as.matrix = TRUE, ret = metrics)[1, ]
  AUCTrain[i] <- auc(roc_object)
  
  # output Test for every split
  predTest <- predict(model, newdata = modelTestData)
  roc_object <- roc(modelOutcomesTest, as.numeric(as.character(predTest)))
  metricsArrayTest[i, ] <- coords(roc_object, x=0.5, transpose=FALSE, as.matrix = TRUE, ret = metrics)[1, ]
  AUCTest[i] <- auc(roc_object)
}

write.csv(MetricsSummary(metricsArrayTrain), "output_mu/metrics_train_SVM.csv")
write.csv(MetricsSummary(metricsArrayTest), "output_mu/metrics_test_SVM.csv")
summary(AUCTrain)
summary(AUCTest)

dfSummary4 <- as.data.frame(MetricsSummary(metricsArrayTrain))
dfSummary4$AUC <- c(summary(AUCTrain))
dfSummary4 <- as.data.frame(t(dfSummary4))
dfSummary4
dfSummary4$model <- c("SVM")
dfSummary4$metric <- rownames(dfSummary4)
dfSummary4$type <- c("Training")
dfSummary4

dfSummary5 <- as.data.frame(MetricsSummary(metricsArrayTest))
dfSummary5$AUC <- c(summary(AUCTest))
dfSummary5 <- as.data.frame(t(dfSummary5))
dfSummary5
dfSummary5$model <- c("SVM")
dfSummary5$metric <- rownames(dfSummary5)
dfSummary5$type <- c("Validation")
dfSummary5

dfSummary6 <- rbind(dfSummary4, dfSummary5)
dfSummary6

png(filename = "figures_mu/SVM_medianvalues_100.png", units = "cm", width=10, height=10, res=300)
ggplot(data=dfSummary6[dfSummary6$metric %in% c("sensitivity", "specificity", "accuracy", "AUC"),], 
       aes(x = metric, y=median, fill=type)) + 
  geom_bar(stat="identity", position=position_dodge(0.6), width = 0.5) +
  ylim(0, 1) +
  scale_fill_manual("type", values = colors) +
  theme(legend.title = element_blank(), axis.title.x = element_blank(), plot.title = element_text(hjust = 0.5), axis.text = element_text(size=12)) +
  ggtitle("SVM") +
  labs(y="median value") +
  theme(legend.position = c(0.98, 0.98), legend.justification = c("right", "top"), legend.box.just = "right", legend.margin = margin(2,4,4,4))
dev.off()

######## Naive Bayes Model #########
library(klaR)
featureNames <- finalFeatureNames[1:4] # number of variables in final model, not too many to keep model simple (Avishek paper An Empirical Approach)
metrics <- c("threshold", "youden", "sensitivity", "specificity", "ppv", "npv", "accuracy")
metricsArrayTrain <- array(NA, 
                           dim = c(100, length(metrics)), 
                           dimnames = list('split' = 1:100, 'metric' = metrics))
metricsArrayTest <- array(NA, 
                          dim = c(100, length(metrics)), 
                          dimnames = list('split' = 1:100, 'metric' = metrics))
AUCTrain <- array(NA,
                  dim = 100,
                  dimnames = list('split' = 1:100))

AUCTest <- array(NA,
                 dim = 100,
                 dimnames = list('split' = 1:100))


for (i in 1:ncol(trainTestSplits)) {
  # Select data
  split <- trainTestSplits[, i]
  modelTrainData <- modelFeatures[split,]
  modelTestData <- modelFeatures[!split,]
  modelOutcomesTrain <- as.factor(modelOutcomes[split])
  modelOutcomesTest <- as.factor(modelOutcomes[!split])
  
  # Apply model to training data, threshold is set with Youden method
  formula <- as.formula(paste('modelOutcomesTrain ~', paste(featureNames, collapse = ' + ')))
  model <- NaiveBayes(formula, data = modelTrainData)
  predTrain <- predict(model)
  roc_object<- roc(modelOutcomesTrain, predTrain$posterior[,2])
  metricsArrayTrain[i, ] <- coords(roc_object, x="best", transpose=FALSE, best.method = "youden", as.matrix = TRUE, ret = metrics)[1, ]
  AUCTrain[i] <- auc(roc_object)
  
  # output Test for every split
  predTest <- predict(model, newdata = modelTestData)
  roc_object <- roc(modelOutcomesTest, predTest$posterior[,2])
  metricsArrayTest[i, ] <- coords(roc_object, x="best", transpose=FALSE, best.method = "youden", as.matrix = TRUE, ret = metrics)[1, ]
  AUCTest[i] <- auc(roc_object)
}

tab <- table(predTest$class, modelOutcomesTest)
caret::confusionMatrix(tab) 

write.csv(MetricsSummary(metricsArrayTrain), "output_mu/metrics_train_NBayes.csv")
write.csv(MetricsSummary(metricsArrayTest), "output_mu/metrics_test_NBayes.csv")
summary(AUCTrain)
summary(AUCTest)

dfSummary7 <- as.data.frame(MetricsSummary(metricsArrayTrain))
dfSummary7$AUC <- c(summary(AUCTrain))
dfSummary7 <- as.data.frame(t(dfSummary7))
dfSummary7
dfSummary7$model <- c("Naïve Bayes")
dfSummary7$metric <- rownames(dfSummary7)
dfSummary7$type <- c("Training")
dfSummary7

dfSummary8 <- as.data.frame(MetricsSummary(metricsArrayTest))
dfSummary8$AUC <- c(summary(AUCTest))
dfSummary8 <- as.data.frame(t(dfSummary8))
dfSummary8
dfSummary8$model <- c("Naïve Bayes")
dfSummary8$metric <- rownames(dfSummary8)
dfSummary8$type <- c("Validation")
dfSummary8

dfSummary9 <- rbind(dfSummary7, dfSummary8)
dfSummary9

png(filename = "figures_mu/NB_medianvalues_100.png", units = "cm", width=10, height=10, res=300)
ggplot(data=dfSummary9[dfSummary9$metric %in% c("sensitivity", "specificity", "accuracy", "AUC"),], 
       aes(x = metric, y=median, fill=type)) + 
  geom_bar(stat="identity", position=position_dodge(0.6), width = 0.5) +
  ylim(0, 1) +
  scale_fill_manual("type", values = colors) +
  theme(legend.title = element_blank(), axis.title.x = element_blank(), plot.title = element_text(hjust = 0.5), axis.text = element_text(size=12)) +
  ggtitle("Naïve Bayes") +
  labs(y="median value") +
  theme(legend.position = c(0.98, 0.98), legend.justification = c("right", "top"), legend.box.just = "right", legend.margin = margin(2,4,4,4))
dev.off()

# The best performing simple model was logistic regression"
####### Develop 1 logistic regression model using the entire set #######
featureNames <- finalFeatureNames[1:4]
formula <- as.formula(paste('modelOutcomes ~', paste(featureNames, collapse = ' + ')))
model2 <- lrm(formula, data = modelData, x=TRUE, y=TRUE) # same but with other function, to get Nagelkerke R2
validate(model2, B = 2000)
# take the index. corrected values
# (Dxy/2) +0.5 = C-index
model <- glm(formula, data = modelData, family = "binomial") # same but with other function, to use with ROC object
predFinal <- predict.glm(model, type = "response")
roc_objectFinal <- roc(modelOutcomes, predFinal)
coordsFinal <- coords(roc_objectFinal, x="best", transpose=FALSE, best.method = "youden", ret = metrics)[1,]
cutoffFinal <- coordsFinal$threshold
aucFinal <- auc(modelOutcomes, predFinal)
coef(model)
summary(model)

png(filename = "figures_mu/Hist_LR_entire_dataset.png", units = "cm", width=14, height=14, res=300)
hist(predFinal, breaks = 30, xlab = "Predicted probability",
     ylim = c(0,30),
     xlim = c(0, 1),
     col = colors[1],
     main = NULL)
abline(v = cutoffFinal, col = colors[1])
dev.off()

ROCPlot(predFinal, modelOutcomes, col="blue")

# Get ROC metrics with pROC package
AUC <- as.numeric(roc_objectFinal$auc)
AUC
ci.auc(roc_objectFinal, method = "bootstrap")

#calculate odds ratio
OddsRatio <- function(logodds) {
  OR <- exp(logodds)
  return(OR)
}
OddsRatio(0.3982657)
# Calculate OR with 95% CI, output lower limit - OR - upper limit
exp(summary(model)$coefficients["(Intercept)",1] + qnorm(c(0.025,0.5,0.975)) * summary(model)$coefficients["(Intercept)",2])
exp(summary(model)$coefficients["flatness",1] + qnorm(c(0.025,0.5,0.975)) * summary(model)$coefficients["flatness",2])
exp(summary(model)$coefficients["cTcT3-4a",1] + qnorm(c(0.025,0.5,0.975)) * summary(model)$coefficients["cTcT3-4a",2])
exp(summary(model)$coefficients["minimum",1] + qnorm(c(0.025,0.5,0.975)) * summary(model)$coefficients["minimum",2])
exp(summary(model)$coefficients["Gearys_C_measure",1] + qnorm(c(0.025,0.5,0.975)) * summary(model)$coefficients["Gearys_C_measure",2])

#calibration
png(filename = "figures_mu/LR_entire_dataset.png", units = "cm", width=14, height=14, res=300)
val.prob.ci.2(y = modelOutcomes, p = predFinal, CL.smooth=TRUE, CL.BT = T)
dev.off()

####### Random forest #######
# see script Python