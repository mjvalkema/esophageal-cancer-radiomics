# Data preparation for model updating
# Authors: M.J. Valkema, L.R. de Ruiter, January 2022

################# Initialization ###################
setwd(dir=getwd())
library(car) # for recoding
library(ggpubr) # for ggarrange

# use dataAll from DataPrepExtVal
source('DataPrepExtVal.R')

# Remove column names with NA's and remove the features that were normalized in the script of external validation.
# In this script normalization will be done using scanner-specific standardization
drop <- colnames(dataAll)[apply(dataAll, 2, anyNA)]
data <- dataAll[,!(names(dataAll) %in% drop)]
drop <- c("joint_maximum_3D_comb_norm", "median_absolute_deviation_norm",
          "joint_entropy_3D_comb_norm", "sum_entropy_3D_comb_norm", 
          "angular_second_moment_3D_comb_norm", "inverse_variance_3D_comb_norm")
data <- data[,!(names(data) %in% drop)] 
# 101 features were calculated by the software as used for the development cohort

# Recode clinical features
data$cT <- as.factor(data$cT)
data$cN <- as.factor(data$cN)
data$cN <- car::recode(data$cN, "c(0)='cN0'; c(1,2)='cN1-2'; c(3)='cN3'; c(9)='cNx'")
data$Gender <- as.factor(data$Gender)
data$Histology <- as.factor(data$Histology)
data$TRG12_TRG34 <- car::recode(data$TRG, "c(1,2)='TRG 1-2'; c(3,4)='TRG 3-4'")
data$TRG12_TRG34 <- as.factor(data$TRG12_TRG34)
data$Manufacturer <- as.factor(data$Manufacturer)
data$ManufacturerModelName <- as.factor(data$ManufacturerModelName)
data$outcome <- car::recode(data$TRG, "c(1,2)='0'; c(3,4)='1'")

# Standardize radiomic features using method published by A. Chatterjee
# Chatterjee, A., et al. Creating robust predictive radiomic models for data from independent institutions using 
# normalization. IEEE Transactions on Radiation and Plasma Medical Sciences 3 (2019): 210-215.
# using ManufacturerModelName as batch
# use only ManufacturerModelNames with >=8 patients
# only normalize once with this script

# Only use data with patients who were scanned on frequently used scanners
smallScannerThreshold <- 8
smallScannerModelNames <- names(table(data$ManufacturerModelName))[table(data$ManufacturerModelName) < smallScannerThreshold]
datStandard <- data[!data$ManufacturerModelName %in% smallScannerModelNames, ] #duplicate data without patients with small scanners

# Split the data for feature standardization
features <- datStandard[, c(1:102, 112, 114)] # Includes "id", all features, "ManufacturerModelName", "PatientID"
clinical <- datStandard[, c(1, 103:116)] # Includes "id", all clinical data, "PatientID"

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

# Load function to make plots
# Visualize features vs. manufacturer
BoxPlotColumnPerManufacturer <- function (column, data, y_label) {
  ggplot(data, aes_string(x='ManufacturerModelName', y=column))+
    geom_boxplot() +
    labs(y=y_label) +
    theme(legend.position="none") +
    xlab(NULL) +
    theme(axis.text.x = element_text(angle = +90, vjust=0, size = 12)) +
    theme(axis.title.y = element_text(size = 12)) +
    scale_x_discrete(labels =c("1080", "Biograph 128 mCT", "Biograph 40 mCT", "Biograph 64 mCT", "Discovery 710", "GEMINI TF TOF 16", "Guardian Body(C)"))
}

png(filename = "figures_mu/Boxplots_Manufacturer_withoutStandardization.png", units = "cm", width=22, height=24, res=300)
figures <- list()
for (i in 1:length(columnsOfInterestNoNorm)) {
  figures[[i]] <- BoxPlotColumnPerManufacturer(columnsOfInterestNoNorm[i], datStandard, columnsOfInterestLabels[i])
}
figure <- ggarrange(plotlist = figures,
                    labels = c("A", "B", "C", "D", "E", "F"),
                    ncol =3, nrow =2,
                    common.legend= TRUE, 
                    legend = "bottom")
figure
dev.off()


RadiomicsFeatureStandardization <- function(features) {
  library(caret)
  process <- preProcess(features, method=c("center", "scale")) # alternative is "range" between 0-1
  normalized <- predict(process, features)
  return(normalized)
}

temp1 <- subset(features, ManufacturerModelName == "Biograph 40_mCT")
temp2 <- subset(features, ManufacturerModelName == "GEMINI TF TOF 16")
temp3 <- subset(features, ManufacturerModelName == "Biograph 128_mCT")
temp4 <- subset(features, ManufacturerModelName == "1080")
temp5 <- subset(features, ManufacturerModelName == "Discovery 710")
temp6 <- subset(features, ManufacturerModelName == "Biograph 64_mCT")
temp7 <- subset(features, ManufacturerModelName == "Guardian Body(C)")

standardizedFeatures1 <- RadiomicsFeatureStandardization(temp1) # using mean 0 and SD 1. Change to range 0-1 if required.
standardizedFeatures2 <- RadiomicsFeatureStandardization(temp2) # using mean 0 and SD 1. Change to range 0-1 if required.
standardizedFeatures3 <- RadiomicsFeatureStandardization(temp3) # using mean 0 and SD 1. Change to range 0-1 if required.
standardizedFeatures4 <- RadiomicsFeatureStandardization(temp4) # using mean 0 and SD 1. Change to range 0-1 if required.
standardizedFeatures5 <- RadiomicsFeatureStandardization(temp5) # using mean 0 and SD 1. Change to range 0-1 if required.
standardizedFeatures6 <- RadiomicsFeatureStandardization(temp6) # using mean 0 and SD 1. Change to range 0-1 if required.
standardizedFeatures7 <- RadiomicsFeatureStandardization(temp7) # using mean 0 and SD 1. Change to range 0-1 if required.

standardizedFeatures <- merge(standardizedFeatures1, standardizedFeatures2, all.y = TRUE, all.x = TRUE)
standardizedFeatures <- merge(standardizedFeatures, standardizedFeatures3, all.y = TRUE, all.x = TRUE)
standardizedFeatures <- merge(standardizedFeatures, standardizedFeatures4, all.y = TRUE, all.x = TRUE)
standardizedFeatures <- merge(standardizedFeatures, standardizedFeatures5, all.y = TRUE, all.x = TRUE)
standardizedFeatures <- merge(standardizedFeatures, standardizedFeatures6, all.y = TRUE, all.x = TRUE)
standardizedFeatures <- merge(standardizedFeatures, standardizedFeatures7, all.y = TRUE, all.x = TRUE)

library(dplyr)
allData <- left_join(standardizedFeatures, clinical, by = "PatientID", suffix = c("", ""))
write.csv(allData, file = "output_mu/allData_used_for_model_updating.csv")
# allData will be used in the script for model updating (createTrainValidate.R)

# Explore the effect of Scanner Specific Standardization
png(filename = "figures_mu/Boxplots_Manufacturer_ScannerSpecificStandardization.png", units = "cm", width=22, height=24, res=300)
figures <- list()
for (i in 1:length(columnsOfInterestNoNorm)) {
  figures[[i]] <- BoxPlotColumnPerManufacturer(columnsOfInterestNoNorm[i], allData, columnsOfInterestLabels[i])
}
figure <- ggarrange(plotlist = figures,
                    labels = c("A", "B", "C", "D", "E", "F"),
                    ncol =3, nrow =2,
                    common.legend= TRUE, 
                    legend = "bottom")
figure
dev.off()

# Normalize using Combat Harmonization, to explore which method of standardization/harmonization works best
library(neuroCombat)
datCombat <- left_join(features, clinical, by = "PatientID", suffix = c("", ""))
dataHarmonized <- neuroCombat(as.data.frame(t(datCombat[,2:102])), batch = as.character(datCombat$ManufacturerModelName), mod = datCombat$TRG12_TRG34, parametric=FALSE)
tempCombat <- as.data.frame(t(dataHarmonized$dat.combat))
tempCombat$PatientID <- datCombat$PatientID
allDataCombat <- merge(tempCombat, clinical, by.x = "PatientID", by.y="PatientID", all.x = TRUE) # merge the data to the clinical

png(filename = "figures_mu/Boxplots_Manufacturer_CombatHarmonization.png", units = "cm", width=22, height=24, res=300)
figures <- list()
for (i in 1:length(columnsOfInterestNoNorm)) {
  figures[[i]] <- BoxPlotColumnPerManufacturer(columnsOfInterestNoNorm[i], allDataCombat, columnsOfInterestLabels[i])
}
figure <- ggarrange(plotlist = figures,
                    labels = c("A", "B", "C", "D", "E", "F"),
                    ncol =3, nrow =2,
                    common.legend= TRUE, 
                    legend = "bottom")
figure
dev.off()