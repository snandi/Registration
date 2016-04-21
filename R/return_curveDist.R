#'@import fdakma
#'@import RFunctionsSN
#'
#'@export

return_curveDist <- function( Mat1, Mat2, argvals, D1 = TRUE ){
  ## This function returns: Similarity index between the mean curves of
  ## two samples of curves

  if(nrow(Mat1) != nrow(Mat2)){
    stop("Both datasets should have same numbers of observations")
  }
  L <- nrow(Mat1)
  if(L != length(argvals)){
    stop("argvals should be the same as the numbers of rows of matrices")
  }
  
  if(D1){
    Dist_Obs <- fdakma::kma.similarity(x.f = argvals, y1.f = rowMeans( Mat1 ), 
                                       x.g = argvals, y1.g = rowMeans( Mat2 ), 
                                       similarity.method = 'd1.pearson')
  } else{
    Dist_Obs <- fdakma::kma.similarity(x.f = argvals, y0.f = rowMeans( Mat1 ), 
                                       x.g = argvals, y0.g = rowMeans( Mat2 ), 
                                       similarity.method = 'd0.pearson')
  }
  return(Dist_Obs)  
}


