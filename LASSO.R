# Model extension using LASSO
# code written by L.R. de Ruiter, D. van Klaveren and M.J. Valkema, June 2022

################# Initialization ###################
setwd(dir=getwd())
source('DataPrepModUpdating.R') # to get allData - this is the standardized data to use with this script, available from DataPrepModUpdating.R
library(caret)
library(glmnet) # LASSO
library(rms)
library(ROCR)
library(pROC) # to compute auc
library(ggplot2)
library(dplyr)
library(car)
source('ROCPlot.R')
library(ggplot2)
colors <- c("#4E79A7","#F28E2B", "#E15759")
metrics <- c("threshold", "youden", "sensitivity", "specificity", "ppv", "npv", "accuracy")

# Define data
#outcomes <- allData[, "outcome"] # TRG 1-2 is 0, TRG 3-4 is 1
outcomes <- allData[, "TRG1_TRG234"] # TRG 1 is 0, TRG 2-3-4 is 1
features <- allData[, 2:102] # all the normalized features without PatientID
features <- features[, order(names(features))] # order alphabetically
clinical <- allData[, c('cT', 'cN_grouped', 'Gender', 'Age', 'Histology')]

modelFeatures <- cbind(clinical, features)
modelOutcomes <- as.numeric(as.character(outcomes))

# Explore correlations in the dataset
corrs <- cor(data.matrix(modelFeatures))
View(rowSums(abs(corrs) > 0.7))
write.csv(corrs,"output_mu/corr.csv")

# Remove highly correlated features
corrFeatures <- findCorrelation(corrs, cutoff = 0.9, exact = TRUE, names = TRUE) # returns names of features to remove
sort(corrFeatures)
modelFeatures <- modelFeatures[, !names(modelFeatures) %in% corrFeatures]
names(modelFeatures)

################# Use LASSO on entire dataset using radiomic and clinical features ###################
set.seed(42)
cv <- cv.glmnet(data.matrix(modelFeatures), modelOutcomes, alpha = 1, family='binomial') # For binary outcomes, we need family='binomial'
png(filename = "figures_mu/lasso_allfeatures.png", units = "cm", width=22, height=17, res=600)
plot(cv)
dev.off()
log(cv$lambda.min)
coef(cv, s = "lambda.min")
preds <- plogis(predict(cv,newx=data.matrix(modelFeatures), s = "lambda.min"))
score <- auc(modelOutcomes, preds[,1])

# LASSO coefficients
print(cv$lambda.min) # min lambda value
coefs <- coef(cv, s = "lambda.min")

# Export model coefficients
coefs_df <- as.data.frame(summary(coefs))
coef_names <- names(modelFeatures)[summary(coefs)$i - 1]
coefs_df$feature <- c('intercept', coef_names)
write.csv(coefs_df[c('feature', 'x')], "output_mu/LASSO_results_radiomics_clinical.csv")

# Threshold and performance metrics
roc_objectLASSO <- roc(modelOutcomes, preds[,1])
coordsLASSO <- coords(roc_objectLASSO, x=.90, input="sensitivity", transpose=FALSE, ret = metrics)[1,] #get specificity at 90% sensitivity
coordsLASSO <- coords(roc_objectLASSO, x= "best", best.method = "youden", as.matrix = TRUE, transpose=FALSE, ret = metrics)[1,] #get metrics at maximized Youden's index
coordsLASSO <- coords(roc_objectLASSO, "all", transpose=FALSE, ret = metrics) # get metrics for all thresholds
coordsLASSO

################# Examine optimism using bootstrapping, model radiomics + clinical variables ###################
getBootstrapSample <- function(features, outcomes) {
  data_size <- length(outcomes)  
  indices <- sample(data_size, data_size, replace=TRUE)
  sample_features <- features[indices,]
  sample_outcomes <- outcomes[indices]
  
  return(list('features'=sample_features,
              'outcomes'=sample_outcomes,
              'indices'=indices))
}

# Bootstrapped LASSO
BOOTSTRAP_REPETITIONS <- 200

sample_aucs <- array(NA, dim = BOOTSTRAP_REPETITIONS)
original_aucs <- array(NA, dim = BOOTSTRAP_REPETITIONS)

for (i in 1:BOOTSTRAP_REPETITIONS){
  # Fit model on sample
  sample <- getBootstrapSample(modelFeatures, modelOutcomes)
  cv <- cv.glmnet(data.matrix(sample$features), sample$outcomes, alpha = 1, family='binomial')
  
  # Compute AUC on sample and complete dataset
  sample$preds <- plogis(predict(cv,newx=data.matrix(sample$features), s = "lambda.min"))
  original_preds <- plogis(predict(cv,newx=data.matrix(modelFeatures), s = "lambda.min"))
  
  sample_aucs[i] <- auc(sample$outcomes, sample$preds)
  original_aucs[i] <- auc(modelOutcomes, original_preds)
  
}

optimism <- mean(sample_aucs - original_aucs)
mean(sample_aucs) # mean sample AUC
mean(original_aucs) # mean original (test performance) AUC
internally_validated_LASSO <- score - optimism
internally_validated_LASSO

################# Reference LASSO on entire dataset using clinical features ###################
cv <- cv.glmnet(data.matrix(clinical), modelOutcomes, alpha = 1, family='binomial')
png(filename = "figures_mu/lasso_clinical.png", units = "cm", width=22, height=17, res=600)
plot(cv)
dev.off()
preds_clinical <- plogis(predict(cv,newx=data.matrix(clinical), s = "lambda.min"))
score_clinical <- auc(modelOutcomes, preds_clinical[,1])

# LASSO coefficients
print(cv$lambda.min) # min lambda value
coefs <- coef(cv, s = "lambda.min")

# Export model coefficients
coefs_df <- as.data.frame(summary(coefs))
coef_names <- names(clinical)[summary(coefs)$i - 1]
coefs_df$feature <- c('intercept', coef_names)
write.csv(coefs_df[c('feature', 'x')], "output_mu/LASSO_results_clinical_only.csv")

# Threshold and performance metrics
roc_objectLASSO_clinical <- roc(modelOutcomes, preds_clinical[,1])
coordsLASSO_clinical <- coords(roc_objectLASSO_clinical, x=.90, input="sensitivity", transpose=FALSE, ret = metrics)[1,] #get specificity at 90% sensitivity
coordsLASSO_clinical <- coords(roc_objectLASSO_clinical, x= "best", best.method = "youden", as.matrix = TRUE, transpose=FALSE, ret = metrics)[1,] #get metrics at maximized Youden's index
coordsLASSO_clinical <- coords(roc_objectLASSO_clinical, "all", transpose=FALSE, ret = metrics) # get metrics for all thresholds
coordsLASSO_clinical

################# Examine optimism using bootstrapping for clinical reference model ###################
BOOTSTRAP_REPETITIONS <- 200

set.seed(42)
sample_aucs_clinical <- array(NA, dim = BOOTSTRAP_REPETITIONS)
original_aucs_clinical <- array(NA, dim = BOOTSTRAP_REPETITIONS)

for (i in 1:BOOTSTRAP_REPETITIONS){
  
  # Fit model on sample
  sample <- getBootstrapSample(clinical, modelOutcomes)
  cv <- cv.glmnet(data.matrix(sample$features), sample$outcomes, alpha = 1)
  
  # Compute AUC on sample and complete dataset
  sample$preds <- plogis(predict(cv,newx=data.matrix(sample$features), s = "lambda.min"))
  original_preds <- plogis(predict(cv,newx=data.matrix(clinical), s = "lambda.min"))
  
  sample_aucs_clinical[i] <- auc(sample$outcomes, sample$preds)
  original_aucs_clinical[i] <- auc(modelOutcomes, original_preds)
  
}

optimism_clinical <- mean(sample_aucs_clinical - original_aucs_clinical)
mean(sample_aucs_clinical) # mean sample AUC
mean(original_aucs_clinical) # mean original AUC

internally_validated_LASSO_clinical <- score_clinical - optimism_clinical
internally_validated_LASSO_clinical


#################  ROC for LASSO with radiomic + clinical features and clinical features only ###################
png(filename = "figures_mu/ROCs_extendedmodel.png", units = "cm", width=22, height=17, res=600)
ROCPlot(preds[,1], modelOutcomes, col=colors[1])
ROCPlot(preds_clinical[,1], modelOutcomes, col=colors[2], add = TRUE)
legend("bottomright", 
       legend = c(paste0("Radiomics + clinical LASSO model,", " ", "AUC:", " ", round(internally_validated_LASSO, digits = 2), " ", "(", round(score, digits = 2), ")"), 
                  paste0("Clinical LASSO model,", " ", "AUC:", " ", round(internally_validated_LASSO_clinical, digits = 2), " ", "(", round(score_clinical, digits = 2), ")")), 
       col = colors, lty = 1)
dev.off()
