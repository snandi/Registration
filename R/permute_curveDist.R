# x1fd <- Regfd.M$regfd
# x2fd <- Regfd.F$regfd
# 
# Nperm = 200
# rangeobs = x1fd$basis$range
# 
# argvals = seq(rangeobs[1], rangeobs[2], length.out = 101)
# Mat1 <-   eval.fd(evalarg = argvals, x1fd)
# Mat2 <-   eval.fd(evalarg = argvals, x2fd)
# 
# q = 0.05
# returnPlot = TRUE
#'@import fda
#'@import fdakma
#'@import ggplot2
#'@import RFunctionsSN
#'
#'@export
permute_kmaSimilarity <- function(Mat1, Mat2, Nperm = 200, argvals, q = 0.05, D1 = TRUE){

  ## This function returns: 
  ##                      overall similarity measure & p-value
  ## This functions requires two matrices of same number of rows, 
  ## each column being a separate curve; 
  ## GOAL: To establish, nonparametrically, by permutations, if these two
  ## groups of curves are statistically different. It estimates the similarity
  ## measure using kma.similarity, the default d1.pearson distance metric,
  ##
  ## The matrices should be of the first derivatives of the curves, for d1.pearson
  
  if(nrow(Mat1) != nrow(Mat2)){
    stop("Both datasets should have same numbers of observations")
  }
  L <- nrow(Mat1)
  if(L != length(argvals)){
    stop("argvals should be the same as the numbers of rows of matrices")
  }
  
  # q <- 1 - q    ## Not needed
  N1 <- ncol(Mat1)
  N2 <- ncol(Mat2)
  # print(paste(N1, N2))
  MatAll <- cbind(Mat1, Mat2)
  
  Dist_null <- rep(0, Nperm)
  
  i <- 2
  ## For loop for permutations
  for (i in 1:Nperm) {
    MatAll_Perm <- MatAll[, sample(N1 + N2)]
    PermMat1 <- MatAll_Perm[, 1:N1]
    PermMat2 <- MatAll_Perm[, N1+(1:N2)]
    rm(MatAll_Perm)
    if(D1){
      Dist_null[i] <- kma.similarity(x.f = argvals, y1.f = rowMeans(PermMat1), 
                                     x.g = argvals, y1.g = rowMeans(PermMat2), 
                                     similarity.method = 'd1.pearson')
    } else{
      Dist_null[i] <- kma.similarity(x.f = argvals, y0.f = rowMeans(PermMat1), 
                                     x.g = argvals, y0.g = rowMeans(PermMat2), 
                                     similarity.method = 'd0.pearson')
    }
  }

  ## Statistics for observed data
  if(D1){
    Dist_Obs <- fdakma::kma.similarity(x.f = argvals, y1.f = rowMeans(Mat1), 
                                       x.g = argvals, y1.g = rowMeans(Mat2), 
                                       similarity.method = 'd1.pearson')
  } else{
    Dist_Obs <- fdakma::kma.similarity(x.f = argvals, y0.f = rowMeans(Mat1), 
                                       x.g = argvals, y0.g = rowMeans(Mat2), 
                                       similarity.method = 'd0.pearson')
  }
  
  pval <- mean(Dist_Obs > Dist_null)  
  qval <- quantile(Dist_null, q)
  
  ## Plotting the data
  Maintitle <- paste('Permutation test, with L1 norm similarity,', Nperm, 'iterations, p-value', pval)
  Color <- 'royalblue1' # 'turquoise2'
  ## If Dist_Obs is outside the null distribution, adjust the Xlim accordingly
  if(inRange(x = Dist_Obs, Range = range(Dist_null))){
    XMin <- min(Dist_null)
    XMax <- max(Dist_null)
    Xlim <- range(Dist_null)
  } else{
    ## print('Using this one')
    XMin <- min(Dist_null, Dist_Obs)
    XMax <- max(Dist_null, Dist_Obs)
    RangeDiff <- diff(range(Dist_null))/10    
    Xlim <- c(XMin - RangeDiff, XMax + RangeDiff)
    ## print(Xlim)
  }
  
  Xlab <- expression(paste('Null distribution of similarity index,          [Reject when ', T[obs], ' < ', T[crit], ' ]'))
  Binwidth <- min(0.01, diff(Xlim)/30)

  Plot_pval <- qplot() + geom_histogram(aes(Dist_null), fill = Color, binwidth = Binwidth, color = Color) + 
    xlim(Xlim) + 
    geom_vline(xintercept = Dist_Obs, colour = 'orangered', size = 1) +
    geom_vline(xintercept = qval, colour = 'gray10', size = 1) +
    annotate(geom = 'text', x = Dist_Obs, y = Inf, vjust = 2, label = 'Obs', col = 'orangered') + 
    annotate(geom = 'text', x = qval, y = Inf, vjust = 4, label = 'Crit', col = 'gray10') +
    ggtitle(Maintitle) + xlab(label = Xlab)

  return(list(pval = pval, qval = qval, Dist_Obs = Dist_Obs, Dist_null = Dist_null, 
              Plot_pval = Plot_pval))
}
