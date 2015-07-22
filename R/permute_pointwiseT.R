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
    # Tnull[i] = max(Tnullvals[, i])   ## This is the test statistic recommended by fda package
    Tnull[i] = median(Tnullvals[, i])   ## This is the test statistic for null distribution
  }
  
  ## Statistics for observed data
  mean1 <- rowMeans(Mat1)
  mean2 <- rowMeans(Mat2)
  var1 <- apply(X=Mat1, MARGIN=1, FUN=var)/N1
  var2 <- apply(X=Mat2, MARGIN=1, FUN=var)/N2
  Tvals <- abs(mean1 - mean2)/sqrt(var1 + var2)
  # Tobs <- max(Tvals)    ## This is the observed value of the test statistic, recommended by fda
  Tobs <- median(Tvals)    ## This is the observed value of the test statistic
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
  
  Ylim <- c(min(Tvals, qvals.pts) - 0.1, max(Tvals, qval) + 0.1)
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

  ## Plotting the pvalue
  Maintitle <- paste('Permutation test, with T-type statistic,', Nperm, 'iterations, p-value', pval)
  Color <- 'royalblue1' # 'turquoise2'
  Binwidth <- round(diff(range(Tnull))/30, 1)
  Xlim <- c(min(Tnull, Tobs)*0.999, max(Tnull, Tobs, 1)*1.001)
  Xlab <- expression(paste('Null distribution of T-type statistic,         [ Reject when ', T[obs], ' > ', T[crit], ' ]'))
  Plot_pval <- qplot() + geom_histogram(aes(Tnull), fill=Color, binwidth=Binwidth, color=Color) + 
    xlim(Xlim) + 
    geom_vline(xintercept=Tobs, colour='orangered', size=1) +
    geom_vline(xintercept=qval, colour='gray10', size=1) +
    annotate(geom='text', x=Tobs, y=Inf, vjust=2, label='Obs', col='orangered') + 
    annotate(geom='text', x=qval, y=Inf, vjust=4, label='Crit', col='gray10') +
    ggtitle(Maintitle) + xlab(label=Xlab)
  
  return(list(pval = pval, qval = qval, Tobs = Tobs, Tnull = Tnull, 
              Tvals = Tvals, Tnullvals = Tnullvals, qvals.pts = qvals.pts, 
              pvals.pts = pvals.pts, TPlot=TPlot, Plot_pval=Plot_pval))
}
