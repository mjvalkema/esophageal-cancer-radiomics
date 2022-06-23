# External validation script
# Authors: M.J. Valkema, L.R. de Ruiter, June 2022

################# Initialization ###################
setwd(dir=getwd())
source('DataPrepExtVal.R') # use script DataPrepExtVal to get the  dataset dataAll: 
library(rms) # to create calibration
library(devtools)
require(devtools)
library(sva)
library(dplyr)
library(ggplot2)
library(ggpubr) # for function ggarrange
library(car) # for recoding
colors <- c("#4E79A7","#F28E2B", "#E15759", "#76b7b2") # for colors in figures

# Baseline characteristics using dataAll of the script DataPrepExtVal
library(tableone)
baselineVars <- c("Gender", "Age", "Histology", "cT_original", "cT", "cN", "cN_grouped", "TRG", "outcome", "Institution", "CRE")
catVars <- c("Gender", "Histology", "cT_original", "cT", "cN", "cN_grouped", "TRG", "outcome","Institution", "CRE")
tableBl <- CreateTableOne(vars = baselineVars, factorVars = catVars, data = dataAll, strata = c("id"),
                          testApprox = chisq.test)
tableBl <- print(tableBl, nonnormal = "Age", quote = FALSE, noSpaces = TRUE)
write.csv(tableBl, file = "output_extval/tableOne.csv")
# Dubbel check statistical tests
chisq.test(table(dataAll$Gender, dataAll$id))
wilcox.test(Age ~ id, data = dataAll)
chisq.test(table(dataAll$Histology, dataAll$id))
chisq.test(table(dataAll$cT, dataAll$id))
chisq.test(table(dataAll$cN_grouped, dataAll$id))
chisq.test(table(dataAll$outcome, dataAll$id))
# Comparison of proportion patients with TRG 1 out of patients with cT3-4a for development and ext. validation cohort respectively
prop.test(x = c(9, 28), n = c(64, 141))
# Comparison of proportion patients with TRG 1 out of patients with cT1-2 for development and ext. validation cohort respectively
prop.test(x = c(7, 12), n = c(9, 41))

# Describe scan parameters for the external validation cohort
tableScannerTypes <- table(dataAll$ManufacturerModelName, dataAll$Manufacturer)
write.csv(tableScannerTypes, file = "output_extval/tableScannerTypes.csv")
scanParmVars <- c("Interval_CRE_PET", "AcquisitionTime", "SliceThickness",
                  "Manufacturer", "ManufacturerModelName","Institution", "Reconstruction_CRE_PET", 
                  "Glucose_CRE_PET", "MBq", "PixelSize")
scanParmVarsCat <- c("Manufacturer", "ManufacturerModelName","Institution", "Reconstruction_CRE_PET", "SliceThickness")
tableScanParam <- CreateTableOne(vars = scanParmVars, factorVars = scanParmVarsCat, data = dataAll[dataAll$id == "External validation cohort",])
tableScanParam <- print(tableScanParam, nonnormal = c("Interval_CRE_PET", "AcquisitionTime", "Glucose_CRE_PET", "MBq","PixelSize"), quote = FALSE, noSpaces = TRUE)
write.csv(tableScanParam, file = "output_extval/tableScanParams.csv")

# Describe feature values for normalized data and non-normalized data
columnsOfInterest = c("joint_maximum_3D_comb_norm", # Pyradiomics original_glcm_MaximumProbability"
                      "median_absolute_deviation_norm", #Pyradiomics original_firstorder_MeanAbsoluteDeviation
                      "joint_entropy_3D_comb_norm", # Pyradiomics original_glcm_JointEntropy"
                      "sum_entropy_3D_comb_norm", #Pyradiomics original_glcm_SumEntropy"
                      "angular_second_moment_3D_comb_norm", #Pyradiomics original_glcm_JointEnergy"
                      "inverse_variance_3D_comb_norm") # Pyrdiomics original_glcm_InverseVariance"

columnsOfInterestNoNorm = c("joint_maximum_3D_comb", # Pyradiomics original_glcm_MaximumProbability"
                            "median_absolute_deviation", #Pyradiomics original_firstorder_MeanAbsoluteDeviation
                            "joint_entropy_3D_comb", # Pyradiomics original_glcm_JointEntropy"
                            "sum_entropy_3D_comb", #Pyradiomics original_glcm_SumEntropy"
                            "angular_second_moment_3D_comb", #Pyradiomics original_glcm_JointEnergy"
                            "inverse_variance_3D_comb") # Pyrdiomics original_glcm_InverseVariance"

columnsOfInterestLabels <- c("joint maximum", 
                             "median absolute deviation", 
                             "joint entropy",
                             "sum entropy",
                             "angular second moment",
                             "inverse variance") 

tableFeatures <- CreateTableOne(vars = columnsOfInterest, data = dataAll,  strata = c("id"))
tableFeatures <- print(tableFeatures, nonnormal = columnsOfInterest, quote = FALSE, noSpaces = TRUE) # Normalized data
write.csv(tableFeatures, file = "output_extval/featureDescriptivesNormalized.csv")

tableFeatures <- CreateTableOne(vars = columnsOfInterestNoNorm, data = dataAll,  strata = c("id"))
tableFeatures <- print(tableFeatures, nonnormal = columnsOfInterestNoNorm, quote = FALSE, noSpaces = TRUE) # Non-normalized data
write.csv(tableFeatures, file = "output_extval/featureDescriptivesWithoutNormalization.csv")

tableFeatures <- CreateTableOne(vars = columnsOfInterestNoNorm, data = dataAll[(dataAll$Manufacturer == "Siemens"),],  strata = c("id"))
tableFeatures <- print(tableFeatures, nonnormal = columnsOfInterestNoNorm, quote = FALSE, noSpaces = TRUE) # Non-normalized data, only SIEMENS features
write.csv(tableFeatures, file = "output_extval/featureDescriptivesWithoutNormalizationSIEMENS.csv")

# Comparison of feature values between development cohort and external validation cohort ("id")
library(ggsignif)
BoxPlotColumnID <- function (column, data, y_label) {
  ggplot(data, aes_string(x='id', y=column))+
    geom_boxplot() +
    labs(y=y_label) +
    theme(legend.position="none") +
    xlab(NULL) +
    theme(axis.text.x = element_text(angle = 0, vjust=0, size = 12)) +
    theme(axis.title.y = element_text(size = 12)) +
    scale_x_discrete(labels =c("devel.", "ext. val.")) +
    geom_signif(comparisons = list(c("Development cohort", "External validation cohort")), 
                test = "wilcox.test", map_signif_level=TRUE, vjust = 0.4, tip_length = 0)
}

png(filename = "figures_extval/Boxplots_with_normalisation.png", units = "cm", width=22, height=17, res=300)
figures <- list()
for (i in 1:length(columnsOfInterest)) {
  figures[[i]] <- BoxPlotColumnID(columnsOfInterest[i], dataAll, columnsOfInterestLabels[i])
}
figure <- ggarrange(plotlist = figures,
                    labels = c("A", "B", "C", "D", "E", "F"),
                    ncol =3, nrow =2,
                    common.legend= TRUE, 
                    legend = "bottom")
figure
dev.off()

###### External validation ########
library(CalibrationCurves)
source('ROCPlot.R')
source('SigmoidFunction.R')
library(caret) # for confusion matrix
library(pROC) # for ROC plot and coordinates
library(ROCR) # for prediction function

outcome <- "TRG1_TRG234"
features <- c("id", columnsOfInterest, "cT", "TRG")

# Choose the data
model_data <- dataAll[dataAll$id == "External validation cohort", c(features, outcome)] # only external validation cohort

### Visualize the data ###
# Boxplots to take a glance at the normalized features for the TRG outcomes
model_data$outcome <- as.factor(model_data$TRG1_TRG234)
model_data$outcome <- car::recode(model_data$outcome, "c(0)='TRG 1'; c(1)='TRG 2-3-4'")

BoxPlotColumn <- function (column, data, y_label) {
  ggplot(data, aes_string(x='outcome', y=column))+
    geom_boxplot() +
    labs(y=y_label) +
    theme(legend.position="none") +
    xlab(NULL) +
    theme(axis.text.x = element_text(angle = 0, vjust=0, size = 12)) +
    theme(axis.title.y = element_text(size = 12)) +
    geom_signif(comparisons = list(c("TRG 1", "TRG 2-3-4")), 
                test = "wilcox.test", map_signif_level=TRUE, vjust = 0.4, tip_length = 0)
}

png(filename = "figures_extval/Boxplots_TRGvsOutcome.png", units = "cm", width=17.4, height=15, res=300)
figures <- list()
for (i in 1:length(columnsOfInterest)) {
  figures[[i]] <- BoxPlotColumn(columnsOfInterest[i], model_data, columnsOfInterestLabels[i])
}
figure <- ggarrange(plotlist = figures,
                    #labels = c("A", "B", "C", "D", "E", "F"),
                    ncol =3, nrow =2,
                    common.legend= TRUE, 
                    legend = "bottom")
figure
dev.off()

##### Externally validate models #####
##### Model validation function ######
ValidateFeature <- function(intercept, coefficientFeature, feature, coefficientT) {
  logodds <- intercept + (feature * coefficientFeature) + coefficientT
  return(logodds)
}

# Source coefficients: Beukinga et al. Radiology 2018.
coefficients = list(A = list(cT12 = 1.0, cT34 = -2.81, intercept = -1.28, radiomics_feature = 0.83),
                    B = list(cT12 = 1.0, cT34 = -2.78, intercept = -1.27, radiomics_feature = -0.74),
                    C = list(cT12 = 1.0, cT34 = -2.76, intercept = -1.27, radiomics_feature = -0.68),
                    D = list(cT12 = 1.0, cT34 = -2.70, intercept = 0.85, radiomics_feature = -0.68),
                    E = list(cT12 = 1.0, cT34 = -2.71, intercept = 0.88, radiomics_feature = 0.74),
                    F = list(cT12 = 1.0, cT34 = -2.92, intercept = 0.90, radiomics_feature = -0.71))
radiomics_features = list(A = "joint_maximum_3D_comb_norm",
                          B = "median_absolute_deviation_norm",
                          C = "joint_entropy_3D_comb_norm",
                          D = "sum_entropy_3D_comb_norm",
                          E = "angular_second_moment_3D_comb_norm",
                          F = "inverse_variance_3D_comb_norm")

# Apply the model to predict TRG 2-3-4 (N.B. opposite of the paper by Beukinga et al., which predicts TRG 1)
ApplyModel <- function(data, coefficients, radiomics_feature_name) {
  data$cT_coefficient <- NA
  data[data$cT == "cT1-2",]$cT_coefficient <- coefficients$cT12
  data[data$cT == "cT3-4a",]$cT_coefficient <- coefficients$cT34
  logodds <- ValidateFeature(intercept = coefficients$intercept, 
                             coefficientFeature = coefficients$radiomics_feature, 
                             feature = data[,radiomics_feature_name], 
                             coefficientT = data$cT_coefficient)
  return(1 - Sigmoid(logodds))
}

metrics <- c("threshold", "youden", "sensitivity", "specificity", "ppv", "npv", "accuracy")

# Choose how metrics are collected: at maximized Youden's index, at sensitivity 90% (preSANO trial), or show at all thresholds
# Maximum Youden's index
CollectMetrics <- function(true_labels, probabilities, model_name) {
  roc_object <- roc(true_labels, probabilities)
  ci.auc(roc_object, method = "bootstrap")
  coords <- coords(roc_object, x="best", transpose=FALSE, best.method = "youden", ret = metrics)
  coords$AUC <-roc_object$auc[1]
  coords$model <- model_name
  return(coords)
}

# At 90% sensitivity
CollectMetrics <- function(true_labels, probabilities, model_name) {
  roc_object <- roc(true_labels, probabilities)
  ci.auc(roc_object, method = "bootstrap")
  coords <- coords(roc_object, x=0.90, input = "sensitivity", transpose=FALSE, ret = metrics)
  coords$AUC <-roc_object$auc[1]
  coords$model <- model_name
  return(coords)
}

# See what the output is when all thresholds are shown
CollectMetrics <- function(true_labels, probabilities, model_name) {
  roc_object <- roc(true_labels, probabilities)
  ci.auc(roc_object, method = "bootstrap")
  coords <- coords(roc_object, "all", transpose=FALSE, ret = metrics)
  coords$AUC <-roc_object$auc[1]
  coords$model <- model_name
  return(coords)
}

# Further define the data
model_data1 <- dataAll[dataAll$id == "Development cohort", c(features, outcome)] # only Groningen cohort
model_data2 <- dataAll[dataAll$id == "External validation cohort", c(features, outcome)] # only external validation cohort
model_data3 <- dataAll[(dataAll$id == "External validation cohort" & dataAll$Manufacturer == "Siemens"), c(features, outcome)] # to test for one vendor
model_data4 <- dataAll[(dataAll$id == "External validation cohort" & dataAll$Histology == "AC"), c(features, outcome)] # to test for adenocarcinoma

# Apply models: calculate probabilities for the development Groningen cohort
model_data1$probabilityA <- ApplyModel(model_data1, coefficients$A, columnsOfInterest[1])
model_data1$probabilityB <- ApplyModel(model_data1, coefficients$B, columnsOfInterest[2])
model_data1$probabilityC <- ApplyModel(model_data1, coefficients$C, columnsOfInterest[3])
model_data1$probabilityD <- ApplyModel(model_data1, coefficients$D, columnsOfInterest[4])
model_data1$probabilityE <- ApplyModel(model_data1, coefficients$E, columnsOfInterest[5])
model_data1$probabilityF <- ApplyModel(model_data1, coefficients$F, columnsOfInterest[6])

# Apply models: calculate probabilities for the external validation cohort
model_data2$probabilityA <- ApplyModel(model_data2, coefficients$A, columnsOfInterest[1])
model_data2$probabilityB <- ApplyModel(model_data2, coefficients$B, columnsOfInterest[2])
model_data2$probabilityC <- ApplyModel(model_data2, coefficients$C, columnsOfInterest[3])
model_data2$probabilityD <- ApplyModel(model_data2, coefficients$D, columnsOfInterest[4])
model_data2$probabilityE <- ApplyModel(model_data2, coefficients$E, columnsOfInterest[5])
model_data2$probabilityF <- ApplyModel(model_data2, coefficients$F, columnsOfInterest[6])

# Apply models: calculate probabilities for the external validation SIEMENS cohort
model_data3$probabilityA <- ApplyModel(model_data3, coefficients$A, columnsOfInterest[1])
model_data3$probabilityB <- ApplyModel(model_data3, coefficients$B, columnsOfInterest[2])
model_data3$probabilityC <- ApplyModel(model_data3, coefficients$C, columnsOfInterest[3])
model_data3$probabilityD <- ApplyModel(model_data3, coefficients$D, columnsOfInterest[4])
model_data3$probabilityE <- ApplyModel(model_data3, coefficients$E, columnsOfInterest[5])
model_data3$probabilityF <- ApplyModel(model_data3, coefficients$F, columnsOfInterest[6])

# Apply models: calculate probabilities for the external validation cohort, adenocarcinoma only
model_data4$probabilityA <- ApplyModel(model_data4, coefficients$A, columnsOfInterest[1])
model_data4$probabilityB <- ApplyModel(model_data4, coefficients$B, columnsOfInterest[2])
model_data4$probabilityC <- ApplyModel(model_data4, coefficients$C, columnsOfInterest[3])
model_data4$probabilityD <- ApplyModel(model_data4, coefficients$D, columnsOfInterest[4])
model_data4$probabilityE <- ApplyModel(model_data4, coefficients$E, columnsOfInterest[5])
model_data4$probabilityF <- ApplyModel(model_data4, coefficients$F, columnsOfInterest[6])

# Collect metrics: save all the performance metrics in one table
# Do for development cohort
tabelSum1 <- CollectMetrics(model_data1$TRG1_TRG234, model_data1$probabilityA, "A")
coords <- CollectMetrics(model_data1$TRG1_TRG234, model_data1$probabilityB, "B")
tabelSum1 <- dplyr::bind_rows(tabelSum1, coords)
coords <- CollectMetrics(model_data1$TRG1_TRG234, model_data1$probabilityC, "C")
tabelSum1 <- dplyr::bind_rows(tabelSum1, coords)
coords <- CollectMetrics(model_data1$TRG1_TRG234, model_data1$probabilityD, "D")
tabelSum1 <- dplyr::bind_rows(tabelSum1, coords)
coords <- CollectMetrics(model_data1$TRG1_TRG234, model_data1$probabilityE, "E")
tabelSum1 <- dplyr::bind_rows(tabelSum1, coords)
coords <- CollectMetrics(model_data1$TRG1_TRG234, model_data1$probabilityF, "F")
tabelSum1 <- dplyr::bind_rows(tabelSum1, coords)
is.num <- sapply(tabelSum1, is.numeric)
tabelSum1[is.num] <- lapply(tabelSum1[is.num], round, 2)
write.csv(tabelSum1, file = "output_extval/tabelSum1_DevelopmentCohort.csv")

# Do for external validation cohort
tabelSum2 <- CollectMetrics(model_data2$TRG1_TRG234, model_data2$probabilityA, "A")
coords <- CollectMetrics(model_data2$TRG1_TRG234, model_data2$probabilityB, "B")
tabelSum2 <- dplyr::bind_rows(tabelSum2, coords)
coords <- CollectMetrics(model_data2$TRG1_TRG234, model_data2$probabilityC, "C")
tabelSum2 <- dplyr::bind_rows(tabelSum2, coords)
coords <- CollectMetrics(model_data2$TRG1_TRG234, model_data2$probabilityD, "D")
tabelSum2 <- dplyr::bind_rows(tabelSum2, coords)
coords <- CollectMetrics(model_data2$TRG1_TRG234, model_data2$probabilityE, "E")
tabelSum2 <- dplyr::bind_rows(tabelSum2, coords)
coords <- CollectMetrics(model_data2$TRG1_TRG234, model_data2$probabilityF, "F")
tabelSum2 <- dplyr::bind_rows(tabelSum2, coords)
is.num <- sapply(tabelSum2, is.numeric)
tabelSum2[is.num] <- lapply(tabelSum2[is.num], round, 2)
write.csv(tabelSum2, file = "output_extval/tabelSum2_extvalcohort.csv")

# Do for external validation cohort with one vendor (SIEMENS)
tabelSum3 <- CollectMetrics(model_data3$TRG1_TRG234, model_data3$probabilityA, "A")
coords <- CollectMetrics(model_data3$TRG1_TRG234, model_data3$probabilityB, "B")
tabelSum3 <- dplyr::bind_rows(tabelSum3, coords)
coords <- CollectMetrics(model_data3$TRG1_TRG234, model_data3$probabilityC, "C")
tabelSum3 <- dplyr::bind_rows(tabelSum3, coords)
coords <- CollectMetrics(model_data3$TRG1_TRG234, model_data3$probabilityD, "D")
tabelSum3 <- dplyr::bind_rows(tabelSum3, coords)
coords <- CollectMetrics(model_data3$TRG1_TRG234, model_data3$probabilityE, "E")
tabelSum3 <- dplyr::bind_rows(tabelSum3, coords)
coords <- CollectMetrics(model_data3$TRG1_TRG234, model_data3$probabilityF, "F")
tabelSum3 <- dplyr::bind_rows(tabelSum3, coords)
is.num <- sapply(tabelSum3, is.numeric)
tabelSum3[is.num] <- lapply(tabelSum3[is.num], round, 2)
write.csv(tabelSum3, file = "output_extval/tabelSum3_extvalcohort_onlySIEMENS.csv")

# Do for external validation cohort with only adenocarcinoma patients
tabelSum4 <- CollectMetrics(model_data4$TRG1_TRG234, model_data4$probabilityA, "A")
coords <- CollectMetrics(model_data4$TRG1_TRG234, model_data4$probabilityB, "B")
tabelSum4 <- dplyr::bind_rows(tabelSum4, coords)
coords <- CollectMetrics(model_data4$TRG1_TRG234, model_data4$probabilityC, "C")
tabelSum4 <- dplyr::bind_rows(tabelSum4, coords)
coords <- CollectMetrics(model_data4$TRG1_TRG234, model_data4$probabilityD, "D")
tabelSum4 <- dplyr::bind_rows(tabelSum4, coords)
coords <- CollectMetrics(model_data4$TRG1_TRG234, model_data4$probabilityE, "E")
tabelSum4 <- dplyr::bind_rows(tabelSum4, coords)
coords <- CollectMetrics(model_data4$TRG1_TRG234, model_data4$probabilityF, "F")
tabelSum4 <- dplyr::bind_rows(tabelSum4, coords)
is.num <- sapply(tabelSum4, is.numeric)
tabelSum4[is.num] <- lapply(tabelSum4[is.num], round, 2)
write.csv(tabelSum4, file = "output_extval/tabelSum4_extvalcohort_onlyAC.csv")

# Make ROCs in one figure for development cohort en external validation cohort (+ Siemens and AC at bottom script, Supplementary Material)
png(filename = "figures_extval/ROCs_develop_extval.png", units = "cm", width=22, height=17, res=600)
par(mfrow = c(2, 3))
ROCPlot(model_data1$probabilityA, model_data1$TRG1_TRG234, col=colors[1], main = "A")
ROCPlot(model_data2$probabilityA, model_data2$TRG1_TRG234, col=colors[2], add = TRUE)
legend("bottomright", 
       legend = c(paste("Development,", "AUC:", tabelSum1$AUC[1]), 
                  paste("Ext.val.,", "AUC:", tabelSum2$AUC[1])), 
       col = colors, lty = 1)
ROCPlot(model_data1$probabilityB, model_data1$TRG1_TRG234, col=colors[1], main = "B")
ROCPlot(model_data2$probabilityB, model_data2$TRG1_TRG234, col=colors[2], add = TRUE)
legend("bottomright", 
       legend = c(paste("Development,", "AUC:", tabelSum1$AUC[2]), 
                  paste("Ext.val.,", "AUC:", tabelSum2$AUC[2])), 
       col = colors, lty = 1)
ROCPlot(model_data1$probabilityC, model_data1$TRG1_TRG234, col=colors[1], main = "C")
ROCPlot(model_data2$probabilityC, model_data2$TRG1_TRG234, col=colors[2], add = TRUE)
legend("bottomright", 
       legend = c(paste("Development,", "AUC:", tabelSum1$AUC[3]), 
                  paste("Ext.val.,", "AUC:", tabelSum2$AUC[3])), 
       col = colors, lty = 1)
ROCPlot(model_data1$probabilityD, model_data1$TRG1_TRG234, col=colors[1], main = "D")
ROCPlot(model_data2$probabilityD, model_data2$TRG1_TRG234, col=colors[2], add = TRUE)
legend("bottomright", 
       legend = c(paste("Development,", "AUC:", tabelSum1$AUC[4]), 
                  paste("Ext.val.,", "AUC:", tabelSum2$AUC[4])), 
       col = colors, lty = 1)
ROCPlot(model_data1$probabilityE, model_data1$TRG1_TRG234, col=colors[1], main = "E")
ROCPlot(model_data2$probabilityE, model_data2$TRG1_TRG234, col=colors[2], add = TRUE)
legend("bottomright", 
       legend = c(paste("Development,", "AUC:", tabelSum1$AUC[5]), 
                  paste("Ext.val.,", "AUC:", tabelSum2$AUC[5])), 
       col = colors, lty = 1)
ROCPlot(model_data1$probabilityF, model_data1$TRG1_TRG234, col=colors[1], main = "F")
ROCPlot(model_data2$probabilityF, model_data2$TRG1_TRG234, col=colors[2], add = TRUE)
legend("bottomright", 
       legend = c(paste("Development,", "AUC:", tabelSum1$AUC[6]), 
                  paste("Ext.val.,", "AUC:", tabelSum2$AUC[6])), 
       col = colors, lty = 1)

dev.off()


##### Calibration #####
# For the development cohort
png(filename = "figures_extval/calplot_dev.png", units = "cm", width=22, height=17, res=600)
par(mfrow = c(2, 3))
val.prob.ci.2(y = model_data1$TRG1_TRG234, p = model_data1$probabilityA, smooth = "loess", CL.smooth=TRUE, main = "A")
val.prob.ci.2(y = model_data1$TRG1_TRG234, p = model_data1$probabilityB, smooth = "loess", CL.smooth=TRUE, main = "B")
val.prob.ci.2(y = model_data1$TRG1_TRG234, p = model_data1$probabilityC, smooth = "loess", CL.smooth=TRUE, main = "C")
val.prob.ci.2(y = model_data1$TRG1_TRG234, p = model_data1$probabilityD, smooth = "loess", CL.smooth=TRUE, main = "D")
val.prob.ci.2(y = model_data1$TRG1_TRG234, p = model_data1$probabilityE, smooth = "loess", CL.smooth=TRUE, main = "E")
val.prob.ci.2(y = model_data1$TRG1_TRG234, p = model_data1$probabilityF, smooth = "loess", CL.smooth=TRUE, main = "F")
dev.off()

# For the external validation cohort
png(filename = "figures_extval/calplot_extval.png", units = "cm", width=22, height=17, res=600)
par(mfrow = c(2, 3))
val.prob.ci.2(y = model_data2$TRG1_TRG234, p = model_data2$probabilityA, smooth = "loess", CL.smooth=TRUE, main = "A")
val.prob.ci.2(y = model_data2$TRG1_TRG234, p = model_data2$probabilityB, smooth = "loess", CL.smooth=TRUE, main = "B")
val.prob.ci.2(y = model_data2$TRG1_TRG234, p = model_data2$probabilityC, smooth = "loess", CL.smooth=TRUE, main = "C")
val.prob.ci.2(y = model_data2$TRG1_TRG234, p = model_data2$probabilityD, smooth = "loess", CL.smooth=TRUE, main = "D")
val.prob.ci.2(y = model_data2$TRG1_TRG234, p = model_data2$probabilityE, smooth = "loess", CL.smooth=TRUE, main = "E")
val.prob.ci.2(y = model_data2$TRG1_TRG234, p = model_data2$probabilityF, smooth = "loess", CL.smooth=TRUE, main = "F")
dev.off()

# Make histograms
# For the development cohort
png(filename = "figures_extval/Histall_DevelopmentCohort.png", units = "cm", width=22, height=14, res=600)
par(mfrow = c(2, 3))
hist(model_data1$probabilityA, breaks = 30, xlab = "Predicted probability",
     ylim = c(0,100),
     xlim = c(0, 1),
     main="A", col = colors[1])
abline(v = tabelSum1$threshold[1], col = colors[1]) # threshold is chosen to get 90% sensitivity
hist(model_data1$probabilityB, breaks = 30, xlab = "Predicted probability",
     ylim = c(0,100),
     xlim = c(0, 1),
     main="B", col = colors[1])
abline(v = tabelSum1$threshold[2], col = colors[1])
hist(model_data1$probabilityC, breaks = 30, xlab = "Predicted probability",
     ylim = c(0,100),
     xlim = c(0, 1),
     main="C", col = colors[1])
abline(v = tabelSum1$threshold[3], col = colors[1])
hist(model_data1$probabilityD, breaks = 30, xlab = "Predicted probability",
     ylim = c(0,100),
     xlim = c(0, 1),
     main="D", col = colors[1])
abline(v = tabelSum1$threshold[4], col = colors[1])
hist(model_data1$probabilityE, breaks = 30, xlab = "Predicted probability",
     ylim = c(0,100),
     xlim = c(0, 1),
     main="E", col = colors[1])
abline(v = tabelSum1$threshold[5], col = colors[1])
hist(model_data1$probabilityF, breaks = 30, xlab = "Predicted probability",
     ylim = c(0,100),
     xlim = c(0, 1),
     main="F", col = colors[1])
abline(v = tabelSum1$threshold[6], col = colors[1])
dev.off()

# For the external validation cohort
# thresholds were determined using threshold at 90% sensitivity
png(filename = "figures_extval/Histall_externalvalidation.png", units = "cm", width=22, height=14, res=600)
par(mfrow = c(2, 3))
hist(model_data2$probabilityA, breaks = 30, xlab = "Predicted probability",
     ylim = c(0,150),
     xlim = c(0, 1),
     main="A", col = colors[1])
abline(v = 0.53, col = colors[1])
hist(model_data2$probabilityB, breaks = 30, xlab = "Predicted probability",
     ylim = c(0,150),
     xlim = c(0, 1),
     main="B", col = colors[1])
abline(v = 0.58, col = colors[1])
hist(model_data2$probabilityC, breaks = 30, xlab = "Predicted probability",
     ylim = c(0,150),
     xlim = c(0, 1),
     main="C", col = colors[1])
abline(v = 0.64, col = colors[1])
hist(model_data2$probabilityD, breaks = 30, xlab = "Predicted probability",
     ylim = c(0,150),
     xlim = c(0, 1),
     main="D", col = colors[1])
abline(v = 0.18, col = colors[1])
hist(model_data2$probabilityE, breaks = 30, xlab = "Predicted probability",
     ylim = c(0,150),
     xlim = c(0, 1),
     main="E", col = colors[1])
abline(v = 0.12, col = colors[1])
hist(model_data2$probabilityF, breaks = 30, xlab = "Predicted probability",
     ylim = c(0,150),
     xlim = c(0, 1),
     main="F", col = colors[1])
abline(v = 0.75, col = colors[1])
dev.off()

# Check relation between predicted probabilities and cT stage
table(model_data2$probabilityA > 0.8, model_data2$cT)
table(model_data2$probabilityB > 0.8, model_data2$cT)
table(model_data2$probabilityC > 0.8, model_data2$cT)
table(model_data2$probabilityD > 0.8, model_data2$cT)
table(model_data2$probabilityE > 0.6, model_data2$cT)
table(model_data2$probabilityF > 0.9, model_data2$cT)

# ROCs with the three model_data's, development cohort, external validation cohort, external validation cohort of one vendor,
# and external validation cohort only adenocarcinoma patients
# Make ROCs in one figure
png(filename = "figures_extval/ROCs_subgroups.png", units = "cm", width=26, height=17, res=900)
par(mfrow = c(2, 3))
ROCPlot(model_data1$probabilityA, model_data1$TRG1_TRG234, col=colors[1], main = "A")
ROCPlot(model_data2$probabilityA, model_data2$TRG1_TRG234, col=colors[2], add = TRUE)
ROCPlot(model_data3$probabilityA, model_data3$TRG1_TRG234, col=colors[3], add = TRUE)
ROCPlot(model_data4$probabilityA, model_data4$TRG1_TRG234, col=colors[4], add = TRUE)
legend("bottomright", 
       legend = c(paste("Development,", "AUC:", tabelSum1$AUC[1]), 
                  paste("Ext.val.,", "AUC:", tabelSum2$AUC[1]), 
                  paste("Ext.val., one vendor,", "AUC:", tabelSum3$AUC[1]),
                  paste("Ext.val., adeno.,", "AUC:", tabelSum4$AUC[1])), 
       col = colors, lty = 1)
ROCPlot(model_data1$probabilityB, model_data1$TRG1_TRG234, col=colors[1], main = "B")
ROCPlot(model_data2$probabilityB, model_data2$TRG1_TRG234, col=colors[2], add = TRUE)
ROCPlot(model_data3$probabilityB, model_data3$TRG1_TRG234, col=colors[3], add = TRUE)
ROCPlot(model_data4$probabilityB, model_data4$TRG1_TRG234, col=colors[4], add = TRUE)
legend("bottomright", 
       legend = c(paste("Development,", "AUC:", tabelSum1$AUC[2]), 
                  paste("Ext.val.,", "AUC:", tabelSum2$AUC[2]), 
                  paste("Ext.val., one vendor,", "AUC:", tabelSum3$AUC[2]),
                  paste("Ext.val., adeno.,", "AUC:", tabelSum4$AUC[2])), 
       col = colors, lty = 1)
ROCPlot(model_data1$probabilityC, model_data1$TRG1_TRG234, col=colors[1], main = "C")
ROCPlot(model_data2$probabilityC, model_data2$TRG1_TRG234, col=colors[2], add = TRUE)
ROCPlot(model_data3$probabilityC, model_data3$TRG1_TRG234, col=colors[3], add = TRUE)
ROCPlot(model_data4$probabilityC, model_data4$TRG1_TRG234, col=colors[4], add = TRUE)
legend("bottomright", 
       legend = c(paste("Development,", "AUC:", tabelSum1$AUC[3]), 
                  paste("Ext.val.,", "AUC:", tabelSum2$AUC[3]), 
                  paste("Ext.val., one vendor,", "AUC:", tabelSum3$AUC[3]),
                  paste("Ext.val., adeno.,", "AUC:", tabelSum4$AUC[3])), 
       col = colors, lty = 1)
ROCPlot(model_data1$probabilityD, model_data1$TRG1_TRG234, col=colors[1], main = "D")
ROCPlot(model_data2$probabilityD, model_data2$TRG1_TRG234, col=colors[2], add = TRUE)
ROCPlot(model_data3$probabilityD, model_data3$TRG1_TRG234, col=colors[3], add = TRUE)
ROCPlot(model_data4$probabilityD, model_data4$TRG1_TRG234, col=colors[4], add = TRUE)
legend("bottomright", 
       legend = c(paste("Development,", "AUC:", tabelSum1$AUC[4]), 
                  paste("Ext.val.,", "AUC:", tabelSum2$AUC[4]), 
                  paste("Ext.val., one vendor,", "AUC:", tabelSum3$AUC[4]),
                  paste("Ext.val., adeno.,", "AUC:", tabelSum4$AUC[4])), 
       col = colors, lty = 1)
ROCPlot(model_data1$probabilityE, model_data1$TRG1_TRG234, col=colors[1], main = "E")
ROCPlot(model_data2$probabilityE, model_data2$TRG1_TRG234, col=colors[2], add = TRUE)
ROCPlot(model_data3$probabilityE, model_data3$TRG1_TRG234, col=colors[3], add = TRUE)
ROCPlot(model_data4$probabilityE, model_data4$TRG1_TRG234, col=colors[4], add = TRUE)
legend("bottomright", 
       legend = c(paste("Development,", "AUC:", tabelSum1$AUC[5]), 
                  paste("Ext.val.,", "AUC:", tabelSum2$AUC[5]), 
                  paste("Ext.val., one vendor,", "AUC:", tabelSum3$AUC[5]),
                  paste("Ext.val., adeno.,", "AUC:", tabelSum4$AUC[5])), 
       col = colors, lty = 1)
ROCPlot(model_data1$probabilityF, model_data1$TRG1_TRG234, col=colors[1], main = "F")
ROCPlot(model_data2$probabilityF, model_data2$TRG1_TRG234, col=colors[2], add = TRUE)
ROCPlot(model_data3$probabilityF, model_data3$TRG1_TRG234, col=colors[3], add = TRUE)
ROCPlot(model_data4$probabilityF, model_data4$TRG1_TRG234, col=colors[4], add = TRUE)
legend("bottomright", 
       legend = c(paste("Development,", "AUC:", tabelSum1$AUC[6]), 
                  paste("Ext.val.,", "AUC:", tabelSum2$AUC[6]), 
                  paste("Ext.val., one vendor,", "AUC:", tabelSum3$AUC[6]),
                  paste("Ext.val., adeno.,", "AUC:", tabelSum4$AUC[6])), 
       col = colors, lty = 1)
dev.off()
dev.off()