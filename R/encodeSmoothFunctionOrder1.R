#' @title Encodes a smooth function into a sequence of states
#' 
#' @description Encodes a smooth function by its first derivatives. At each point on
#' the abscissa the function is encoded into a binary variable. These points are represented as
#' sequence of states
#' 
#' @param smoothFunc A smooth function, already converted into an \code{fda::fd} object. 
#' @param abscissa A vector representing the x-axis or abscissa
#' 
#' @details To use this function for a vector use the function \code{createSmoothFD} to convert to a 
#' \code{fd} object.
#' For a smooth function \eqn{f(x)}, the two possible states are as follows:
#' \itemize{
#'    \item \eqn{A=(1): (f'>0)}
#'    \item \eqn{B=(0): (f'\leq 0)}
#' }
#' 
#' @author Subhrangshu Nandi, PhD Statistics; snandi@wisc.edu or nands31@gmail.com
#' 
#' @return A \code{string} of states representing the encoded form of the smooth function, example "AADDCCCB"
#' 
#' @examples
#'  data( growth, package = 'fda')
#'  Mat1 <- growth[['hgtm']]
#'  Arguments <- growth[['age']]
#'  Mean <- rowMeans( Mat1 )
#' MeanFD <- createSmoothFD(
#'     curvesToSmooth = Mean,
#'     abscissa       = Arguments,
#'     lambdas        = exp( -5:5 ),
#'     basisBreakFreq = 3,
#'     basisOrder     = 4,
#'     pbasis         = NULL )
#'  MeanStates <- encodeSmoothFunctionOrder1( smoothFunc = (MeanFD$curvesSmooth)$fd, abscissa = Arguments )
#' 
#' @import fda
#' @export

encodeSmoothFunctionOrder1 <- function( smoothFunc, abscissa ){
  V1 <- V2 <- NULL
  
  if( class( smoothFunc ) != 'fd' ){
    stop( 'smoothFunc needs to be of class fd' )
  }
  
  encodedFunc <- as.data.frame( cbind( 
    D1 = eval.fd( evalarg = abscissa, fdobj = smoothFunc, Lfdobj = 1 )
  ) )
  encodedFunc <- within( data = encodedFunc, {
    D1 = ( round( V1, 4 ) > 0 )
  })
  
  AllStates <- getEncodingStatesOrder1()
  encodedFunc$State <- sapply( X = encodedFunc[, 'D1' ], FUN = getStateName, AllStates = AllStates )
  StateVector <- encodedFunc$State
  names( StateVector ) <- abscissa
  return( StateVector )  
}

getEncodingStatesOrder1 <- function(){
  States <- list( c( 'TRUE' ), 
                  c( 'FALSE' )
  )
  
  names( States ) <- c( 'A', 'B' )
  return( States )
}

getStateName <- function( Row, AllStates ){
  StateName <- ''
  Order <- length( AllStates[[1]] )
  for( i in 1:length( AllStates ) ){
    if ( Row == AllStates[[i]] ){
      StateName <- names( AllStates[i] ) 
    }
  }
  return( StateName )
}

