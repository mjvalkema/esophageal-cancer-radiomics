# Source: Chatterjee et al.: Empirical approach for avoiding false discoveries when applying high-dimensional radiomics,
# IEE transactions on radiation and plasma medical sciences, vol. 3, no.2, March 2019
# Original code of above source rewritten by L.R. de Ruiter and M.J. Valkema, November 2021

# N.B. The dataset  is split in a Training set (for developing) and Validation set (for internal validation). 
# N.B. The nomenclature in the functions is Train:Test. For the current dataset, this means Train:Validate since an independent Test set was not available.

RadiomicsFeatureSelection <- function(features, outcomes) {
  ################# createTrainValidate ###############
  # Goal: to get rid of bad features and create 100 train:validation splits
  outcomes <- as.numeric(outcomes)
  
  # Clean up feature set
  thresholdVariance <- 0.05 # default 0.05, minimum relative variation a feature should have 
  thresholdSameValue <- 0.3 # default 0.3, maximum fraction of patients for whom feature can have same value
  
  goodFeatures <- NULL
  badFeatures <- NULL
  for (feature in names(features)) {
    temp <- features[[feature]] #this selects 1 feature for all rows (patients)
    if (sum(is.na(temp)) == 0 & 
        abs(sd(temp)/mean(temp)) > thresholdVariance & 
        sd(temp) > 10^-10 &
        max(table(temp)) / length(temp) < thresholdSameValue
    ) {
      goodFeatures <- c(goodFeatures, feature)
    } else {
      badFeatures <- c(badFeatures, feature)
    }
  }
  features <- features[,goodFeatures]
  print(c("length badFeatures:", length(badFeatures)))
  print(c("names badFeatures:", names(badFeatures)))
  
  # Create train-validate splits for the Clearance step, corrTest. train:validate == 1:1
  trainTestSplits <- CreateTrainTestSplits(outcomes, nTrainTestSplits = 100, trainTestFraction = 0.5)

  ################# corrTest ###############
  # Optional script
  # For each split i in the 100 splits there is S1,i (train) and S2,i (test) for each feature
  # Only consider features with absolute correlation with the outcome >0.2
  ## S1,i - correlation with outcome and S2,i - correlation with outcome
  ## Calculate Rp,i of that
  ## Create a vector with 100 values of Rp,i
  
  Rp_all <- c() #vector with 100 values of Rp,i
  relDiff_all <- c()
  n_features <- c()
  
  for (split in trainTestSplits) {
    
    correlations_train <- c()
    correlations_test <- c()
    
    for (feature in features) {
      corr_train <- cor(x = feature[split], y = outcomes[split], use = "everything", method = "pearson")
      if (abs(corr_train) >= 0.2) {
        corr_test <- cor(x=feature[!split], y=outcomes[!split], use = "everything", method = "pearson")
        correlations_train <- c(correlations_train, corr_train)
        correlations_test <- c(correlations_test, corr_test)
      }
    }
    
    if (length(correlations_train) > 0) {
      Rp <- cor(correlations_train, correlations_test, method = "pearson")
      relDiff <- (correlations_train - correlations_test) / correlations_train #variance between corr train and validate
    } else {
      Rp <- NA
      relDiff <- NA
    }

    Rp_all <- c(Rp_all, Rp) #for every split
    relDiff_all <- c(relDiff_all, relDiff) #for every split
    n_features <- c(n_features, length(correlations_train)) #per split, how many features have corr >0.2
  }
  
  # Check consistency for train:validate sets for features correlating >0.2 with outcome
  mean(Rp_all, na.rm = TRUE) # should be between 0.5 and 1
  sd(Rp_all, na.rm = TRUE)
  mean(relDiff_all, na.rm = TRUE) # mean relative uncertainty, ideal is 0
  sd(relDiff_all)
  print(paste('Mean Rp:', format(round(mean(Rp_all, na.rm = TRUE),2))))
  print(paste('Standard deviation Rp:', format(round(sd(Rp_all, na.rm = TRUE),2))))
  print(paste('Mean relative uncertainty:', format(round(mean(relDiff_all, na.rm = TRUE),2))))
  print(paste('Standard deviation relative uncertainty:', format(round(sd(relDiff_all, na.rm = TRUE),2))))
  
  # use case might be a productive one if mean Rp_all is between 0.5-1 and mean rel_Diff_all is <1
  png(filename = "figures_mu/feature_selection_histograms.png", units = "cm", width=22, height=7, res=600)
  par(mfrow = c(1, 3))
  hist(Rp_all, breaks=30, xlab = "Rp", main = "", xlim = range(0,1), ylim = range(0,50), col = colors[1])
  hist(relDiff_all, breaks=30, xlab = "relative difference of R (delta R (S1, S2)", main = "", xlim = range(-1, 2.5), ylim = range(0,200), col = colors[1])
  hist(n_features, breaks = 30, xlab = "number of promising features per split", main = "", xlim = range(0, 100), ylim = range(0,50), col = colors[1]) # keep track of the number of promising features in train set per split
  dev.off()
  
  ################# featureExploration and aucExploration ###############
  ## Change train:validation splits to 2:1
  trainTestSplits <- CreateTrainTestSplits(outcomes, nTrainTestSplits = 100, trainTestFraction = 0.67)
  
  # Feature Exploration with AUCs for every split in the training set, rows = good features, columns = splits
  # Get a sense of stability of the AUC value
  arrayAUC <- array(NA, dim=c(length(goodFeatures), nTrainTestSplits))
  
  for (n.split in 1:nTrainTestSplits) {
    split <- trainTestSplits[, n.split]
    for (n.feature in 1:ncol(features)) { #loopen over aantal features
      feature <- features[, n.feature] # select column number to get a feature
      AUC <- auc(response = outcomes[split], predictor = feature[split], quiet=T)
      corr <- cor(x=feature[split], y=outcomes[split], use = "everything", method = "spearman")
      if (AUC < 0.5 & corr < 0) {
        AUC <- (1 - AUC)
      }
      arrayAUC[n.feature, n.split] <- AUC
    }
  }
  
  # Apply over rows (dus per feature) to compute the metrics for the stability of AUC in the training set
  means <- apply(arrayAUC, MARGIN=1, FUN=mean)
  mins <- apply(arrayAUC, MARGIN=1, FUN=min)
  maxes <- apply(arrayAUC, MARGIN=1, FUN=max)
  medians <- apply(arrayAUC, MARGIN=1, FUN=median)
  lowers <- apply(arrayAUC, MARGIN=1, FUN=quantile, probs = 0.025)
  uppers <- apply(arrayAUC, MARGIN=1, FUN=quantile, probs = 0.975)
  
  png(filename = "figures_mu/feature_selection_AUCs_over100splits.png", units = "cm", width=22, height=14, res=600)
  par(mfrow = c(2, 2))
  plot(medians, lowers, xlim = c(0.5, 1), ylim = c(0.5, 1), xlab = "median AUC per feature", ylab = "P2.5 AUC per feature")
  plot(medians, mins, xlim = c(0.5, 1), ylim = c(0.5, 1), xlab = "median AUC per feature", ylab = "min AUC per feature")
  plot(medians, uppers, xlim = c(0.5, 1), ylim = c(0.5, 1), xlab = "median AUC per feature", ylab = "P97.5 AUC per feature")
  plot(medians, maxes, xlim = c(0.5, 1), ylim = c(0.5, 1), xlab = "median AUC per feature", ylab = "max AUC per feature")
  dev.off()
  
  ################# featureSelection ###############
  # We use only the training set for selection of features, we don't use splits
  # Keep features based on threshold 1  
  # Drop features based on threshold 2
  threshold1 <- 0.55 # median over 100 training sets should have at least AUC of 0.55
  threshold2 <- 0.7 # features cannot be correlated with each other > 0.7 (Pearson)
  
  correlationMatrix <- abs(cor(features[trainTestSplits$V2,])) #make correlation matrix for training set, e.g. split1
  
  keep <- (medians > threshold1)
  drop <- rep(FALSE, length(keep))
  
  for (m in 1:length(features)) {
    for (n in 1:length(features)) {
      # Only check features in upper half of correlation matrix, excluding diagonal
      if (m <= n) {
        next
      }
      # Only check pairs of features which are both in keep
      if (keep[n] & keep[m]) {
        # Stop if either feature has already been dropped
        if (drop[n] | drop[m]) {
          next
        } else {
          # Compare features based on correlation and median AUC
          if (correlationMatrix[m, n] > threshold2) {
            if (medians[m] > medians[n]) {
              drop[n] <- TRUE
            } else {
              drop[m] <- TRUE
            }
          }
        }
      }
    }
  }
  
  selectedFeatures <- features[, (keep & !drop)]
  return(selectedFeatures)
}

CreateTrainTestSplits <- function(outcomes, nTrainTestSplits, trainTestFraction, seed=3) {
  
  # Convert outcomes to factor
  if (class(outcomes) != 'factor') {
    warning('Converting outcomes to factor...')
    outcomes <- as.factor(outcomes)
  }
  
  # Create train-validation splits to be used in the Clearance Step and Feature selection scripts
  set.seed(seed)
  indicesTrainTestSplits <- createDataPartition(outcomes, p = trainTestFraction, list = F, times = nTrainTestSplits)
  
  # Make an array in which for each split patients in training are TRUE and patients for testing are FALSE 
  trainTestSplits <- array(NA, dim = c(length(outcomes), nTrainTestSplits))
  for (i in 1:nTrainTestSplits) {
    split <- indicesTrainTestSplits[,i] # column of the split
    trainTestSplits[,i] <- (1:length(outcomes) %in% split) # checks if patient is in the split, outputs TRUE or FALSE
  }
  trainTestSplits <- as.data.frame(trainTestSplits) #each split is V1, V2, etc
  print(paste('Dimensions train:', length(outcomes[trainTestSplits$V1]))) # dimensions Training set
  print(paste('Dimensions validation:', length(outcomes[!trainTestSplits$V1]))) # dimensions Validation set
  
  return(trainTestSplits)
}
