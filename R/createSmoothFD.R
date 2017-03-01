#'@import fda
#'@export

########################################################################## 
## Create a fd object of a curve, after smoothing it with bsplines
## This works even if Curve is a matrix
########################################################################## 
createSmoothFD <- function( 
  curvesToSmooth, 
  abscissa, 
  lambdas        = exp( -5:5 ), 
  basisBreakFreq = 3,
  basisOrder     = 4,
  pbasis         = NULL
 ){
  
  curvesToSmoothMat <- cbind( curvesToSmooth )
  if( ncol( curvesToSmoothMat ) > 1 ){
    stopifnot( length( abscissa ) == nrow( curvesToSmoothMat ) )
  } else{
    stopifnot( length( abscissa ) == length( curvesToSmoothMat ) )
  }
  
  basisLength <- length( abscissa )/basisBreakFreq
  if( is.null( pbasis ) ){
    basisBreaks <- seq( from = min( abscissa ), to = max( abscissa ), length.out = basisLength )
    pbasis <- create.bspline.basis( rangeval = range( abscissa ), norder = basisOrder, breaks = basisBreaks ) 
  }

  gcvs <- rep( 0, length( lambdas ) )

  for( i in 1:length( lambdas ) ){ 
    pPar <- fdPar( pbasis, int2Lfd( 2 ), lambdas[i] )
    gcvs[i] <- mean( smooth.basis( argvals = abscissa, y = curvesToSmooth, fdParobj = pPar )$gcv )
  }
  best <- which.min( gcvs )
  lambdaBest <- lambdas[best]

  pPar <- fdPar( pbasis, int2Lfd( 2 ), lambda = lambdaBest )
  curvesSmooth <- smooth.basis( argvals = abscissa, y = curvesToSmooth, fdParobj = pPar )
  return( list( curvesSmooth = curvesSmooth, lambdaBest = lambdaBest ) )
}
