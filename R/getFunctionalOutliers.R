#' @title Gets functional outliers from a sample of curves
#' 
#' @description Gets functional outliers from a sample of curves using methods described in fda.usc
#' 
#' @param Curves A matrix (or dataframe) of curves, with each column being a separate curve
#' @param Xaxis The abscissa
#' @param Names list(main = 'Main', xlab = 'PixelPosition', ylab = 'Intensity')
#' @param DepthType Either of the four \code{c( 'FM', 'Mode', 'RTukey', 'RProj' )}
#' @param N_Bootstrap Number of boostrap samples. Defaults to 500.
#' @param Trim Whether to Trim the samples or not. 
#' @param TrimPct How much to Trim the samples or not
#' 
#' @return A vector of column names that are detected as outliers.
#' 
#' @references 
#' \itemize{
#'    \item Febrero-Bande, M., & Oviedo de la Fuente, M. (2012). Statistical computing in functional data analysis: the R package fda. usc. Journal of Statistical Software, 51(4), 1-28.
#'    \item Cuevas A, Febrero M, Fraiman R. 2006. On the use of bootstrap for estimating functions with functional data. Computational Statistics and Data Analysis 51: 1063-1074.
#'    \item Febrero-Bande, M., Galeano, P., and Gonzalez-Manteiga, W. (2008). Outlier detection in functional data by depth measures with application to identify abnormal NOx levels. Environmetrics 19, 4, 331-345.
#'    \item Febrero-Bande, M., Galeano, P. and Gonzalez-Manteiga, W. (2007). A functional analysis of NOx levels: location and scale estimation and outlier detection. Computational Statistics 22, 3, 411-427.
#'    \item Febrero-Bande, M., Oviedo de la Fuente, M. (2012). Statistical Computing in Functional Data Analysis: The R Package fda.usc. Journal of Statistical Software, 51(4), 1-28. 
#'    \item \url{http://www.jstatsoft.org/v51/i04/}
#' }
#' 
#' @author Subhrangshu Nandi, PhD Statistics, UW Madison; snandi@wisc.edu or nands31@gmail.com
#' 
#' @seealso \code{outliers.depth.trim}, \code{outliers.depth.pond}
#' 
#' @examples 
#' data( growth, package = 'fda' )
#' Mat1 <- growth[['hgtm']]
#' Arguments <- growth[['age']]
#' getFunctionalOutliers( 
#'   Curves = Mat1, 
#'   Xaxis = Arguments, 
#'   Names = list(main = 'Main', xlab = 'PixelPosition', ylab = 'Intensity'),
#'   DepthType = 'FM',
#'   N_Bootstrap = 500,
#'   Trim = 'Yes',
#'   TrimPct = 0.05
#' )
#' 
#' @keywords functional-data outliers depth
#' 
#' @import fda
#' @import fda.usc
#' @export

###############################################################################
## This function is to get functional outliers, using the package fda.usc
## This should replace fn_getFunctionalOutliers in fn_Library_Simulation.R
## and fn_getOutliers in fn_Library_QualityScore.R
###############################################################################
getFunctionalOutliers <- function( 
  Curves, 
  Xaxis, 
  Names = list(main = 'Main', xlab = 'PixelPosition', ylab = 'Intensity'),
  DepthType = c( 'FM', 'Mode', 'RTukey', 'RProj' ),
  N_Bootstrap = 500,
  Trim = c( 'Yes', 'No' ),
  TrimPct = 0.05
 ){

  Curves.fdata <- fdata( 
    mdata = t( Curves ), 
    argvals = Xaxis,
    rangeval = NULL, 
    names = Names,
    fdata2d = FALSE
  )
  set.seed( seed = 100 )
  if( Trim == 'Yes' ){

    if( DepthType == 'FM' ){
      Outlier <- outliers.depth.trim( fdataobj = Curves.fdata, nb = N_Bootstrap, smo = 0, trim = TrimPct, dfunc = depth.FM )
    } else if( DepthType == 'Mode' ){
      Outlier <- outliers.depth.trim( fdataobj = Curves.fdata, nb = N_Bootstrap, smo = 0, trim = TrimPct, dfunc = depth.mode )
    } else if( DepthType == 'RTukey' ){
      Outlier <- outliers.depth.trim( fdataobj = Curves.fdata, nb = N_Bootstrap, smo = 0, trim = TrimPct, dfunc = depth.RT )
    } else if( DepthType == 'RProj' ){
      Outlier <- outliers.depth.trim( fdataobj = Curves.fdata, nb = N_Bootstrap, smo = 0, trim = TrimPct, dfunc = depth.RP )
    } else{
      Outlier <- outliers.depth.trim( fdataobj = Curves.fdata, nb = N_Bootstrap, smo = 0, trim = TrimPct, dfunc = depth.FM )
    } 
  } else{
    
    if( DepthType == 'FM' ){
      Outlier <- outliers.depth.pond( fdataobj = Curves.fdata, nb = N_Bootstrap, smo = 0, dfunc = depth.FM )
    } else if( DepthType == 'Mode' ){
      Outlier <- outliers.depth.pond( fdataobj = Curves.fdata, nb = N_Bootstrap, smo = 0, dfunc = depth.mode )
    } else if( DepthType == 'RTukey' ){
      Outlier <- outliers.depth.pond( fdataobj = Curves.fdata, nb = N_Bootstrap, smo = 0, dfunc = depth.RT )
    } else if( DepthType == 'RProj' ){
      Outlier <- outliers.depth.pond( fdataobj = Curves.fdata, nb = N_Bootstrap, smo = 0, dfunc = depth.RP )
    } else{
      Outlier <- outliers.depth.pond( fdataobj = Curves.fdata, nb = N_Bootstrap, smo = 0, dfunc = depth.FM )
    } 
  }
  
  return( Outlier )
}
###############################################################################
