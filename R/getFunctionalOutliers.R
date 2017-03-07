#'@import fda
#'@import fda.usc

#'@export

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
