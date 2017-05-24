#' @title RMSE of a sample of curves with a template
#' 
#' @description This function takes in a sample of curves in the form of a \code{matrix} or \code{data.frame}
#' with each column being a separate curve, and returns the average RMSE (root mean squared error) between
#' the curves and a template.
#' 
#' @param curves A \code{matrix} or \code{data.frame} of curves, each column being a separate curve
#' @param template A template curve. If \code{NULL} the template is estimated as \code{rowMeans(curves)}
#' @param Round An integer value to round off the average rmse.
#' 
#' @author Subhrangshu Nandi; PhD, Statistics; snandi@wisc.edu or nands31@gmail.com
#' 
#' @return Returns a list
#' \item{rmseVector}{A vector of rmse between curves and template}
#' \item{rmseMean}{Average of \code{rmseVector}}
#' 
#' @examples
#' data( growth, package = 'fda' )
#' Mat1 <- growth[['hgtm']]
#' rmseMat1 <- rmseCurves( curves = Mat1 )
#' 
#' @keywords rmse
#'  
#' @importFrom Metrics rmse
#' @export
#' 

rmseCurves <- function( curves, template = NULL, Round = 4 ){
  
  stopifnot( class( curves ) == 'data.frame' | class( curves ) == 'matrix' )
  
  if( is.null( template ) ){
    template <- rowMeans( curves )
  }
  
  rmseVector <- apply( X = curves, MARGIN = 2, FUN = Metrics::rmse, actual = template )
  rmseMean <- round( mean( rmseVector ), Round )
  
  rmseList <- list( rmseVector = rmseVector, rmseMean = rmseMean )
  return( rmseList )
}
