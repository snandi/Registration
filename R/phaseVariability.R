#' @title Phase variability of a sample of curves
#' 
#' @description This function estimates the phase variability of a sample of curves that are considered
#' random realizations of a smooth function
#' 
#' @param curves A \code{matrix} or \code{data.frame} of curves, each column representing a separate curve
#' @param template A template curve. If \code{NULL} the template is estimated as \code{rowMeans(curves)}
#' 
#' @author Subhrangshu Nandi; PhD, Statistics; snandi@wisc.edu or nands31@gmail.com
#' 
#' @return phaseVariability \eqn{\omega}, where \eqn{0 \leq \omega \leq 1}
#' 
#' @details This is the first mathematical definition of phase variability. It uses similarity index 
#' (\eqn{\rho}), defined in \code{fdakma::kma.similarity}. For a sample of \eqn{n} curves that are random 
#' realizations of a smooth template curve, phase variability is defined as 
#' \eqn{\omega = \frac{1 - \rho_m}{2}}, where \eqn{\rho_m = \frac{1}{n}\sum\rho_i}, \eqn{\rho_i} being the 
#' similarity between curve \eqn{i} and the template.
#' 
#' @examples
#' data( growth, package = 'fda' )
#' Mat1 <- growth[['hgtm']]
#' phaseVar <- phaseVariability( curves = Mat1 )
#' 
#' @keywords phase-variability similarity 
#'  
#' @import RFunctionsSN 
#' @export
#' 

phaseVariability <- function( curves, template = NULL ){
  
  curvesNormalized <- normalizeMatrixWithRowSD( matrixData = curves )
  curvesNormalized <- normalizeMatrix( matrixData = curvesNormalized, withSD = F )
  
  if( is.null( template ) ){
    template <- rowMeans( curvesNormalized )
  }
  templateNormalized <- normalize_vector( Vector = template, withSD = FALSE )
  
  abscissa <- 1:length( templateNormalized )

  phaseSimilarity <- getKmaSimilarityWithTemplate(
    Mat               = curvesNormalized,
    Template          = templateNormalized,
    Xaxis             = abscissa,
    similarity.method = "d1.pearson",
    Deriv             = FALSE
  )
  phaseVariability <- 0.5*( 1 - mean( phaseSimilarity ) ) ## Multiplied by 0.5 because -1 < phaseSim < 1
                                                          ## To ensure 0 < phaseVar < 1
  return( phaseVariability )
}
