# Define a function to make ROCs using performance function from ROCR package
ROCPlot <- function(pred, truth, ...) {
  predob = prediction(pred, truth)
  perf = performance(predob, "tpr", "fpr")
  plot(perf, ...)
  area <- auc(truth, pred)
  area <- format(round(area, 4), nsmall = 4)
  #text(x=0.8, y=0.1, labels = paste("AUC =", area))
  
  # the reference x=y line
  segments(x0=0, y0=0, x1=1, y1=1, col="gray", lty=2)
}
