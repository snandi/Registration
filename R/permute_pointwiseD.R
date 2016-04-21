#'@import fda
#'@import fdakma
#'@import ggplot2

#'@export

permute_pointwiseD <- function(
  Mat1, 
  Mat2, 
  Nperm = 200, 
  argvals, 
  q = 0.05,
  returnPlot = TRUE,
  TitleText = ''
){
  ## This function returns: 
  ##                      pointwise D statistic values, p-values
  ##                      overall D value & p-value
  ##                      plot of pointwise D values
  ## This function is similar to permute_pointwiseT, without the denominator 
  ## of the T-type test statistic
  
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
  
  Dnull <- rep(0, Nperm)
  Dnullvals <- matrix(data = 0, nrow = L, ncol = Nperm)
  
  i <- 1
  ## For loop for permutations
  for (i in 1:Nperm) {
    MatAll_Perm <- MatAll[, sample(N1 + N2)]
    tmean1 <- rowMeans(MatAll_Perm[, 1:N1]) # apply(X=tXmat[], MARGIN=1, FUN=mean)
    tmean2 <- rowMeans(MatAll_Perm[, N1+(1:N2)]) # apply(X=tXmat[, n1 + (1:n2)], MARGIN=1, FUN=mean)
    Dnullvals[, i] = abs(tmean1 - tmean2)
    Dnull[i] = max(Dnullvals[, i])  
  }
  
  ## Statistics for observed data
  mean1 <- rowMeans(Mat1)
  mean2 <- rowMeans(Mat2)
  Dvals <- abs(mean1 - mean2)
  Dobs <- max(Dvals)    ## This is the observed value of the test statistic, recommended by fda
  # Tobs <- median(Tvals)    ## This is the observed value of the test statistic
  pval <- mean(Dobs < Dnull)  
  qval <- quantile(Dnull, q)
  
  pvals.pts <- apply(X = (Dvals < Dnullvals), MARGIN=1, FUN=mean)  ## Pointwise p-value
  qvals.pts <- apply(X = Dnullvals, MARGIN = 1, FUN = quantile, q)    ## Pointwise quantile value
  
  ## Plotting the result of pointwise tests
  DataToPlot <- as.data.frame(rbind(cbind(argvals, DValues = as.vector(Dvals), DValueType='Observed D'), 
                                    cbind(argvals, DValues = as.vector(qval), DValueType='Critical D'), 
                                    cbind(argvals, DValues = as.vector(qvals.pts), DValueType='Pointwise Critical D')))
  DataToPlot <- within(data=DataToPlot,{
    argvals <- as.numeric(as.vector(argvals))
    DValues <- as.numeric(as.vector(DValues))
    DValueType <- factor(DValueType)
  })
  
  Ylim <- c(min(Dvals, qvals.pts) - 0.1, max(Dvals, qval) + 0.1)
  MainTitle <- paste('Pointwise D statistics for', Nperm, 'Permutations', '\n', TitleText)
  if(returnPlot == TRUE){
    DPlot <- qplot(x = argvals, y = DValues, data = DataToPlot) + 
      ylim(Ylim) + 
      geom_point(aes(color = DValueType)) +
      geom_line(aes(color = DValueType)) + 
      ggtitle(label = MainTitle) +
        theme(legend.position = 'top')
      
  } else{
    DPlot <- 0
  }

  ## Plotting the pvalue
  Maintitle <- paste('Permutation test, with D (abs diff),', Nperm, 'iterations, p-value', pval,
                     '\n', TitleText)
  Color <- 'royalblue1' # 'turquoise2'
  Binwidth <- round( diff( range( Dnull ) )/30, 1 )
  Xlim <- c( min( Dnull, Dobs )*0.999, max( Dnull, Dobs, 1 )*1.001 )
  Xlab <- expression(paste('Null distribution of D (abs diff),         [ Reject when ', D[obs], ' > ', D[crit], ' ]'))
  Plot_pval <- qplot() + geom_histogram(aes( Dnull ), fill = Color, binwidth = Binwidth, color = Color) + 
    xlim( Xlim ) + 
    geom_vline( xintercept = Dobs, colour = 'orangered', size = 1 ) +
    geom_vline( xintercept = qval, colour = 'gray10', size = 1 ) +
    annotate( geom = 'text', x = Dobs, y = Inf, vjust = 2, label = 'Obs', col = 'orangered' ) + 
    annotate( geom = 'text', x = qval, y = Inf, vjust = 4, label = 'Crit', col = 'gray10') +
    ggtitle( Maintitle ) + xlab( label = Xlab )
  
  return(list(
    pval      = pval, 
    qval      = qval, 
    Dobs      = Dobs, 
    Dnull     = Dnull, 
    Dvals     = Dvals, 
    Dnullvals = Dnullvals, 
    qvals.pts = qvals.pts, 
    pvals.pts = pvals.pts, 
    DPlot     = DPlot, 
    Plot_pval = Plot_pval
  ))
}

################## For debugging ##################
# require(fda)
# require(ggplot2)
# data(growth)
# Mat1 <- growth[['hgtm']]
# Mat2 <- growth[['hgtf']]
# Arguments <- growth[['age']]
# PermTestResults <- permute_pointwiseD(
#   Mat1 = Mat1, 
#   Mat2 = Mat2, 
#   Nperm = 200,
#   argvals = Arguments,
#   q = 0.05,
#   returnPlot = TRUE,
#   TitleText = 'Growths of boys and girls (Berkeley growth data)'
# )
# names(PermTestResults)
##################################################
