Call_DeCAF <- function(dataSet) {

    load("../R/DeCAF/decaf_classifier.Rdata")
    classifier = decaf_classifier$classifier2

  source("../R/DeCAF/decaf_functions.R")
  
  rownames(dataSet$ex) <- make.names(dataSet$featInfo$SYMBOL,unique = TRUE)
  predictions <- apply_decaf(data = dataSet$ex, classifier = classifier)
  
  dataSet$sampInfo <- cbind(dataSet$sampInfo, predictions)
  return(dataSet)
}
