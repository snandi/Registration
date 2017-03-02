#'@import fda
#'@import fdakma
#'@import matrixStats
#'@import robustX
#'@export

###############################################################################
## This function returns a weighted consensus for iterated registration
###############################################################################
getWeightedConsensus <- function( Curves, abscissa ){
  medianCurves <- robustX::L1median( t( Curves ) )$estimate
  
  options( warn = -1 )
  Sim_All <- getKmaSimilarityWithTemplate(
    Mat               = Curves, 
    Template          = medianCurves,
    Xaxis             = abscissa,
    similarity.method = c("d1.pearson"),
    Deriv             = FALSE
  )
  options( warn = 0 )
  
  Weights <- Sim_All
  Weights[ Weights < 0 ] <- 0
  medianWeighted <- rowWeightedMedians( x = Curves, w = Weights )
  meanWeighted <- rowWeightedMeans(x = Curves, w = Weights)
  
  return( list( Mean_Weighted = meanWeighted, Median_Weighted = medianWeighted ) )
}
