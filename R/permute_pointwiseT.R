# x1fd <- Regfd.M$regfd
# x2fd <- Regfd.F$regfd
# 
# Nperm = 5000
# rangeobs = x1fd$basis$range
# 
# argvals = seq(rangeobs[1], rangeobs[2], length.out = 101)
# Mat1 <-   eval.fd(evalarg=argvals, x1fd)
# Mat2 <-   eval.fd(evalarg=argvals, x2fd)
# 
# q = 0.05
# returnPlot = TRUE

permute_pointwiseT <- function(Mat1, Mat2, Nperm=200, argvals, q=0.05, returnPlot=TRUE){
  require(fda)
  require(ggplot2)
  ## This function returns: 
  ##                      pointwise T statistic values, p-values
  ##                      overall T value & p-value
  ##                      plot of pointwise T values
  ## This function is adapted from fda:tperm.fd
  ## This functions requires two matrices of same number of rows, 
  ## each column being a separate curve; 
  ## GOAL: To establish, nonparametrically, by permutations if these two
  ## groups of curves are statistically different
  if(nrow(Mat1) != nrow(Mat2)){
    stop("Both datasets should have same numbers of observations")
  }
  L <- nrow(Mat1)
  if(L != length(argvals)){
    stop("argvals should be the same as the numbers of rows of matrices")
  }
  
  q <- 1 - q
  N1 <- ncol(Mat1)
  N2 <- ncol(Mat2)
  MatAll <- cbind(Mat1, Mat2)
  
  Tnull <- rep(0, Nperm)
  Tnullvals <- matrix(data=0, nrow=L, ncol=Nperm)
  
  i <- 1
  ## For loop for permutations
  for (i in 1:Nperm) {
    MatAll_Perm <- MatAll[, sample(N1 + N2)]
    tmean1 <- rowMeans(MatAll_Perm[, 1:N1]) # apply(X=tXmat[], MARGIN=1, FUN=mean)
    tmean2 <- rowMeans(MatAll_Perm[, N1+(1:N2)]) # apply(X=tXmat[, n1 + (1:n2)], MARGIN=1, FUN=mean)
    tvar1 <- apply(X=MatAll_Perm[, 1:N1], MARGIN=1, FUN=var)/N1
    tvar2 <- apply(X=MatAll_Perm[, N1+(1:N2)], MARGIN=1, FUN=var)/N2
    Tnullvals[, i] = abs(tmean1 - tmean2)/sqrt(tvar1 + tvar2)
    Tnull[i] = max(Tnullvals[, i])   ## This is the test statistic for null distribution
  }
  
  ## Statistics for observed data
  mean1 <- rowMeans(Mat1)
  mean2 <- rowMeans(Mat2)
  var1 <- apply(X=Mat1, MARGIN=1, FUN=var)/N1
  var2 <- apply(X=Mat2, MARGIN=1, FUN=var)/N2
  Tvals <- abs(mean1 - mean2)/sqrt(var1 + var2)
  Tobs <- max(Tvals)    ## This is the observed value of the test statistic
  pval <- mean(Tobs < Tnull)  
  qval <- quantile(Tnull, q)
  
  pvals.pts <- apply(X=(Tvals < Tnullvals), MARGIN=1, FUN=mean)  ## Pointwise p-value
  qvals.pts <- apply(X=Tnullvals, MARGIN=1, FUN=quantile, q)    ## Pointwise quantile value
  
  ## Plotting the result of pointwise tests
  DataToPlot <- as.data.frame(rbind(cbind(argvals, TValues=as.vector(Tvals), TValueType='Observed T'), 
                                    cbind(argvals, TValues=as.vector(qval), TValueType='Critical T'), 
                                    cbind(argvals, TValues=as.vector(qvals.pts), TValueType='Pointwise Critical T')))
  DataToPlot <- within(data=DataToPlot,{
    argvals <- as.numeric(as.vector(argvals))
    TValues <- as.numeric(as.vector(TValues))
    TValueType <- factor(TValueType)
  })
  
  Ylim <- c(min(Tvals, qvals.pts) - 0.1, max(Tobs, qval) + 0.1)
  MainTitle <- paste('Pointwise T statistics for', Nperm, 'Permutations')
  if(returnPlot==TRUE){
    TPlot <- qplot(x=argvals, y=TValues, data=DataToPlot) + 
      ylim(Ylim) + 
      geom_point(aes(color=TValueType)) +
      geom_line(aes(color=TValueType)) + 
      ggtitle(label=MainTitle)
  } else{
    TPlot <- 0
  }

  return(list(pval = pval, qval = qval, Tobs = Tobs, Tnull = Tnull, 
              Tvals = Tvals, Tnullvals = Tnullvals, qvals.pts = qvals.pts, 
              pvals.pts = pvals.pts, TPlot=TPlot))
}
