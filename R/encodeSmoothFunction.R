#' @title Encodes a smooth function into a sequence of states
#' 
#' @description Encodes a smooth function by its values, first and second derivatives. At each point on
#' the abscissa the function is encoded into a three dimensional point. These points are represented as
#' sequence of states
#' 
#' @param smoothFunc A smooth function, already converted into an \code{fda::fd} object. 
#' @param abscissa A vector representing the x-axis or abscissa
#' 
#' @details To use this function for a vector use the function \code{createSmoothFD} to convert to a 
#' \code{fd} object.
#' For a smooth function \eqn{f(x)}, the eight possible states are as follows:
#' \itemize{
#'    \item \eqn{A=(1,1,1): (f>0, f'>0, f''>0)}
#'    \item \eqn{B=(1,0,1): (f>0, f'\leq 0, f''>0)}
#'    \item \eqn{C=(1,1,0): (f>0, f'> 0, f''\leq 0)}
#'    \item \eqn{D=(1,0,0): (f>0, f'\leq 0, f''\leq 0)}
#'    \item \eqn{E=(0,1,1): (f\leq 0, f'> 0, f''>0)}
#'    \item \eqn{G=(0,0,1): (f\leq 0, f'\leq 0, f''> 0)}
#'    \item \eqn{H=(0,1,0): (f\leq 0, f'> 0, f''\leq 0)}
#'    \item \eqn{I=(0,0,0): (f\leq 0, f'\leq 0, f''\leq 0)}
#' }
#' 
#' @author Subhrangshu Nandi; snandi@wisc.edu or nands31@gmail.com
#' 
#' @return A \code{string} states representing the encoded form of the smooth function, example "AAHHFFFGGH"
#' 
#' @examples
#'  data( growth, package = 'fda')
#'  Mat1 <- growth[['hgtm']]
#'  Arguments <- growth[['age']]
#'  Mean <- rowMeans( Mat1 )
#'  MeanFD <- createSmoothFD( 
#'      curvesToSmooth = Mean,
#'      abscissa       = Arguments,
#'      lambdas        = exp( -5:5 ),
#'      basisBreakFreq = 3,
#'      basisOrder     = 4,
#'      pbasis         = NULL )
#'  MeanStates <- encodeSmoothFunction( smoothFunc = (MeanFD$curvesSmooth)$fd, abscissa = Arguments )
#' 
#' @import fda
#' @export

encodeSmoothFunction <- function( smoothFunc, abscissa ){
  V1 <- V2 <- V3 <- NULL
  
  if( class( smoothFunc ) != 'fd' ){
    stop( 'smoothFunc needs to be of class fd' )
  }
  
  encodedFunc <- as.data.frame( cbind( 
    D0 = eval.fd( evalarg = abscissa, fdobj = smoothFunc, Lfdobj = 0 ), 
    D1 = eval.fd( evalarg = abscissa, fdobj = smoothFunc, Lfdobj = 1 ), 
    D2 = eval.fd( evalarg = abscissa, fdobj = smoothFunc, Lfdobj = 2 )
  ) )
  encodedFunc <- within( data = encodedFunc, {
    D0 = ( round( V1, 4 ) > 0 )
    D1 = ( round( V2, 4 ) > 0 )
    D2 = ( round( V3, 4 ) > 0 )
  })

  AllStates <- getEncodingStates()
  encodedFunc$State <- apply( X = encodedFunc[, c( 'D0', 'D1', 'D2' ) ], MARGIN = 1, 
                                FUN = getStateName, AllStates = AllStates )
  StateVector <- encodedFunc$State
  names( StateVector ) <- abscissa
  return( StateVector )  
}

getEncodingStates <- function(){
  States <- list( c( 'TRUE', 'TRUE', 'TRUE' ), 
                  c( 'TRUE', 'FALSE', 'TRUE' ), 
                  c( 'TRUE', 'TRUE', 'FALSE' ), 
                  c( 'TRUE', 'FALSE', 'FALSE' ), 
                  c( 'FALSE', 'TRUE', 'TRUE' ), 
                  c( 'FALSE', 'FALSE', 'TRUE' ), 
                  c( 'FALSE', 'TRUE', 'FALSE' ), 
                  c( 'FALSE', 'FALSE', 'FALSE' ) 
  )
  
  names( States ) <- c( 'A', 'B', 'C', 'D', 'E', 'G', 'H', 'I' )
  return( States )
}

getStateName <- function( Row, AllStates ){
  StateName <- ''
  Order <- length( AllStates[[1]] )
  for( i in 1:length( AllStates ) ){
    if ( sum( unlist( Row ) == AllStates[[i]] ) == Order ){
      StateName <- names( AllStates[i] ) 
    }
  }
  return( StateName )
}

