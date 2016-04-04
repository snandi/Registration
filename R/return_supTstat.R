#'@export

########################################################################
## Returns sup of T-type statistic between two matrices
########################################################################
return_supTstat <- function(Mat1, Mat2){
  Mat1 <- as.matrix(Mat1)
  Mat2 <- as.matrix(Mat2)
  
  N1 <- ncol(Mat1)
  N2 <- ncol(Mat2)
  
  mean1 <- rowMeans(Mat1)
  mean2 <- rowMeans(Mat2)
  var1 <- apply(X=Mat1, MARGIN=1, FUN=var)/N1
  var2 <- apply(X=Mat2, MARGIN=1, FUN=var)/N2
  
  Tvals <- abs(mean1 - mean2)/sqrt(var1 + var2)
  TSup <- max(Tvals)
  TMedian <- median(Tvals)
  TMean <- mean(Tvals)
  
  return(list(TSup = TSup, TMedian = TMedian, TMean = TMean))
}
