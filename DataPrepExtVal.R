# Data preparation script: conbining development cohort with external validation cohort
# Authors: M.J. Valkema, L.R. de Ruiter, June 2022

################# Initialization ###################
rm(list=ls()) # clear global environment
setwd(dir=getwd())
library(dplyr)
library(car) # for recoding

###### Development cohort  #####
dataGroningen <- read.csv(file="data/GroningenData.csv", header=TRUE)

# Remove baseline and delta features containing T0 and clean up names of features
index <- grep("T0", colnames(dataGroningen))
T0 <- names(dataGroningen[index])
dataGClean <- dataGroningen[,!(names(dataGroningen) %in% T0)]
colnames(dataGClean) <- sapply(strsplit(colnames(dataGClean), split = "_T2"), "[[" , 1)
dataGroningen <- dataGClean

# Define outcome
CreateTarget <- function(x, group0, group1) {
  if (x %in% group0) return(0)
  if (x %in% group1) return(1)
  stop("Error: value found that is not specified")
}

dataGroningen$TRG1_TRG234 <- apply(dataGroningen$response %o% 1, 1, CreateTarget, c(1), c(2,3,4,5))
dataGroningen$outcome <- as.factor(dataGroningen$TRG1_TRG234)
dataGroningen$outcome <- car::recode(dataGroningen$TRG1_TRG234, "c(0)='TRG 1'; c(1)='TRG 2-3-4'")
dataGroningen$Gender <- car::recode(dataGroningen$geslacht, "c(0)=1 ; c(1)=2")
dataGroningen$Age <- dataGroningen$age
dataGroningen$Histology <- car::recode(dataGroningen$histolog, "c(0)='AC' ; c(1)='SCC'")
dataGroningen$cT <- dataGroningen$cT_stage_TNM7
dataGroningen$cN <- dataGroningen$cN_stage_TNM7
dataGroningen$TRG <- car::recode(dataGroningen$response, "c(5)=4")
dataGroningen$Manufacturer <- 'Siemens'
dataGroningen$ManufacturerModelName <- 'Biograph 64_mCT'
dataGroningen$Institution <- 'UMCG'
dataGroningen$PatientID <- c(1:73)
dataGroningen$PatientID <- as.character(dataGroningen$PatientID)

# Clean up features
drop <- c("tumor_length", "lokal", "response", "geslacht", "age", "histolog", "cT_stage_TNM7", "cN_stage_TNM7")
dataGroningen <- dataGroningen[,!(names(dataGroningen) %in% drop)] 

# Normalise features using min-max scaling (X-min(X))/(max(X)-min(X))
## - If only data is supplied, use regular min-max scaling
## - If also min and max values are supplied, use those
NormalizeFeature <- function(X, min=NULL, max=NULL) {
  if (is.null(min) & is.null(max)) {
    return((X - min(X)) / (max(X) - min(X)))
  }
  if (is.null(min) | is.null(max)) {
    return(error("Supply both min and max or neither"))
  }
  return((X - min) / (max - min))
}

dataGroningen$joint_maximum_3D_comb_norm <- NormalizeFeature(dataGroningen$joint_maximum_3D_comb)
dataGroningen$median_absolute_deviation_norm <- NormalizeFeature(dataGroningen$median_absolute_deviation)
dataGroningen$joint_entropy_3D_comb_norm <- NormalizeFeature(dataGroningen$joint_entropy_3D_comb)
dataGroningen$sum_entropy_3D_comb_norm <- NormalizeFeature(dataGroningen$sum_entropy_3D_comb)
dataGroningen$angular_second_moment_3D_comb_norm <- NormalizeFeature(dataGroningen$angular_second_moment_3D_comb)
dataGroningen$inverse_variance_3D_comb_norm <- NormalizeFeature(dataGroningen$inverse_variance_3D_comb)

##### Load data of external validation cohort #####
features <- read.csv(file="data/finalDB_Transposed.csv", header=TRUE)

# Clean up feature list, removing information on descriptives of feature extraction. Remove avg features (should be only the combined features)
index <- grep("initial_image", colnames(features))
inImg <- names(features[index])
dataClean <- features[,!(names(features) %in% inImg)]

index <- grep("ROI", colnames(dataClean))
ROI <- names(dataClean[index])
dataClean <- dataClean[,!(names(dataClean) %in% ROI)]

index <- grep("interpolated_image", colnames(dataClean))
intImg <- names(dataClean[index])
dataClean <- dataClean[,!(names(dataClean) %in% intImg)]

index <- grep("_avg", colnames(dataClean))
avg <- names(dataClean[index])
dataClean <- dataClean[,!(names(dataClean) %in% avg)]

features <- dataClean

clinical <- read.csv(file="data/PET_post_nCRT/clinical_post.csv",header=TRUE)
clinical$Manufacturer <- car::recode(clinical$Manufacturer, "c('GE MEDICAL SYSTEMS')='GE' ; c('Philips Medical Systems')='Philips' ; c('SIEMENS')='Siemens'")
clinical$ManufacturerModelName <- car::recode(clinical$ManufacturerModelName, "c('Biograph40_mCT')='Biograph 40_mCT' ; c('Biograph128_mCT')='Biograph 128_mCT'")
clinical$TRG1_TRG234 <- apply(clinical$TRG %o% 1, 1, CreateTarget, c(1), c(2,3,4))
clinical$outcome <- as.factor(clinical$TRG1_TRG234)
clinical$outcome <- car::recode(clinical$TRG1_TRG234, "c(0)='TRG 1'; c(1)='TRG 2-3-4'")

metadata <- read.csv(file="data/PET_post_nCRT/metadata_scans_PET_post_withMBq.csv",header=TRUE)
metadata <- subset(metadata, select=c("PatientID", "MBq", "PixelSize"))

#Feature normalization of external validation cohort using statistics of Groningen
features$joint_maximum_3D_comb_norm <- NormalizeFeature(features$joint_maximum_3D_comb, min(dataGroningen$joint_maximum_3D_comb), max(dataGroningen$joint_maximum_3D_comb))
features$median_absolute_deviation_norm <- NormalizeFeature(features$median_absolute_deviation, min(dataGroningen$median_absolute_deviation), max(dataGroningen$median_absolute_deviation))
features$joint_entropy_3D_comb_norm <- NormalizeFeature(features$joint_entropy_3D_comb, min(dataGroningen$joint_entropy_3D_comb), max(dataGroningen$joint_entropy_3D_comb))
features$sum_entropy_3D_comb_norm <- NormalizeFeature(features$sum_entropy_3D_comb, min(dataGroningen$sum_entropy_3D_comb), max(dataGroningen$sum_entropy_3D_comb))
features$angular_second_moment_3D_comb_norm <- NormalizeFeature(features$angular_second_moment_3D_comb, min(dataGroningen$angular_second_moment_3D_comb), max(dataGroningen$angular_second_moment_3D_comb))
features$inverse_variance_3D_comb_norm <- NormalizeFeature(features$inverse_variance_3D_comb, min(dataGroningen$inverse_variance_3D_comb), max(dataGroningen$inverse_variance_3D_comb))

clinical <- merge(x = clinical, y = metadata, by.x = "PatientID", all.x = TRUE)
data <- merge(x = features, y = clinical, by.x = "PatientID", all.x = TRUE)

#### Combine the entire two datasets ####
dataAll <- dplyr::bind_rows(dataGroningen, data, .id = "id")
dataAll$id <- car::recode(dataAll$id, "c(1)='Development cohort'; c(2)='External validation cohort'")
dataAll$Gender <- car::recode(dataAll$Gender, "c(1)='Male'; c(2)='Female'")
dataAll$cT_original <- as.factor(dataAll$cT)
dataAll$cT <- car::recode(dataAll$cT, "c(1,2)='cT1-2'; c(3,4)='cT3-4a'")
dataAll$cN_grouped <- as.factor(car::recode(dataAll$cN, "c(0)='cN0'; c(1,2,3,9)='cN+'"))
dataAll$cN <- car::recode(dataAll$cN, "c(9)='cNx'")
dataAll$AcquisitionTime <- (dataAll$AcquisitionTime / 3600) # convert to minutes instead of seconds
write.csv(dataAll, file = "output/dataAll.csv")
