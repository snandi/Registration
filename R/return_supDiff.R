#'@export

########################################################################
## Returns sup of D statistic (max pointwise diff) between two matrices
########################################################################
return_supDiff <- function(Mat1, Mat2){
  Mat1 <- as.matrix(Mat1)
  Mat2 <- as.matrix(Mat2)
  
  mean1 <- rowMeans(Mat1)
  mean2 <- rowMeans(Mat2)
  
  Dvals <- abs(mean1 - mean2)
  DSup <- max(Dvals)
  DMedian <- median(Dvals)
  DMean <- mean(Dvals)
  
  return(list(DSup = DSup, DMedian = DMedian, DMean = DMean))
}

