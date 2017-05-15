#' @title Create a smooth functional data object using bsplines
#' 
#' @description For a curve of a matrix this function creates a smooth functional data object as per 
#' the definition in R package fda, using bspline basis. 
#' 
#' @param curvesToSmooth A matrix or a vector to smooth
#' @param abscissa The x-axis. Should of of the same length as the vector curvesToSmooth or the 
#' number of rows in matrix curvesToSmooth
#' @param lambdas The different roughness penalty to try during smoothing, for the parameter 
#' \code{lambda} in \code{fda:fdPar}
#' @param basisBreakFreq Defaults to 3. More frequent could overfit, less frequent could over-smooth.
#' @param basisOrder Defaults to 4. If more noisy, then this can be higher.
#' @param pbasis User-defined \code{pbasis} if user-defined basisObject is desired. Else defaults to NULL.
#' 
#' @return A list of two items
#' \item{curvesSmooth}{FDA object of smoothed curves}
#' \item{lambdaBest}{lambda associated with minimum GCV}
#' 
#' @references Graves, S., Hooker, G., & Ramsay, J. (2009). Functional data analysis with R and MATLAB.
#' 
#' @author Subhrangshu Nandi, PhD Statistics, UW Madison; snandi@wisc.edu or nands31@gmail.com
#' 
#' @examples 
#' data( growth, package = 'fda' )
#' Mat1 <- growth[['hgtm']]
#' Mat1Smooth <- createSmoothFD( 
#'    curvesToSmooth = Mat1, 
#'    abscissa = 1:nrow( Mat1 ), 
#'    lambdas = exp( -5:5 ), 
#'    basisBreakFreq = 3, 
#'    basisOrder = 4, 
#'    pbasis = NULL )
#' #plot( (Mat1Smooth$curvesSmooth)$fd )
#' 
#' @seealso \code{fda::create.basis}, \code{fda::create.bspline.basis}
#' 
#' @import fda
#' @export

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
