# x1fd <- Regfd.M$regfd
# x2fd <- Regfd.F$regfd
# 
# Nperm = 200
# rangeobs = x1fd$basis$range
# 
# argvals = seq(rangeobs[1], rangeobs[2], length.out = 101)
# Mat1 <-   eval.fd(evalarg=argvals, x1fd)
# Mat2 <-   eval.fd(evalarg=argvals, x2fd)
# 
# q = 0.05
# returnPlot = TRUE

permute_kmaSimilarity <- function(Mat1, Mat2, Nperm=20000, argvals, q=0.05){
  require(fda)
  require(fdakma)
  ## This function returns: 
  ##                      overall similarity measure & p-value
  ## This functions requires two matrices of same number of rows, 
  ## each column being a separate curve; 
  ## GOAL: To establish, nonparametrically, by permutations, if these two
  ## groups of curves are statistically different. It estimates the similarity
  ## measure using kma.similarity, the default d1.pearson distance metric,
  ##
  ## The matrices should be of the first derivatives of the curves
  
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
  print(paste(N1, N2))
  MatAll <- cbind(Mat1, Mat2)
  
  Dist_null <- rep(0, Nperm)
  
  i <- 2
  ## For loop for permutations
  for (i in 1:Nperm) {
    MatAll_Perm <- MatAll[, sample(N1 + N2)]
    PermMat1 <- MatAll_Perm[, 1:N1]
    PermMat2 <- MatAll_Perm[, N1+(1:N2)]
    rm(MatAll_Perm)
    Dist_null[i] <- kma.similarity(x.f=argvals, y1.f=rowMeans(PermMat1), 
                               x.g=argvals, y1.g=rowMeans(PermMat2), 
                   similarity.method='d1.pearson')
  }

  ## Statistics for observed data
  Dist_Obs <- fdakma::kma.similarity(x.f=argvals, y1.f=rowMeans(Mat1), 
                            x.g=argvals, y1.g=rowMeans(Mat2), 
                            similarity.method='d1.pearson')
  
  pval <- mean(Dist_Obs > Dist_null)  
  qval <- quantile(Dist_null, q)
  
  ## Plotting the data
  Maintitle <- paste('Permutation test, with L1 norm similarity,', Nperm, 'iterations, p-value', pval)
  Color <- 'royalblue1' # 'turquoise2'
  Xlim <- c(min(Dist_null, Dist_Obs)*0.999, 1)
  Plot_pval <- qplot() + geom_histogram(aes(Dist_null), binwidth=0.001, fill=Color, color=Color) + 
    xlim(Xlim) + 
    geom_vline(xintercept=Dist_Obs, colour='orangered', size=1) +
    geom_vline(xintercept=qval, colour='gray10', size=1) +
    annotate(geom='text', x=Dist_Obs, y=Inf, vjust=2, label='Obs', col='orangered') + 
    annotate(geom='text', x=qval, y=Inf, vjust=4, label='Crit', col='gray10') +
    ggtitle(Maintitle)

  return(list(pval = pval, qval = qval, Dist_Obs = Dist_Obs, Dist_null = Dist_null, 
              Plot_pval=Plot_pval))
}
