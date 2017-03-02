#'@import fda
#'@import ggplot2

############################################################################### 
## This function conducts iteration registration to weighted average template
## estimated at every iteration
############################################################################### 
registerIterated <- function( 
  dataToRegister, 
  Lambdas_ConstrainedWarping, 
  abscissaFrom,
  abscissaTo,
  abscissaIncrement,
  basisOrder = 4,
  basisBreakFreq = 3,
  Lambdas_Roughness = exp( -5:5 )
){
  abscissa <- seq( from = abscissaFrom, to = abscissaTo, by = abscissaIncrement ) 
  
  dataToRegisterSmooth <- createSmoothFD( 
    curvesToSmooth = dataToRegister, 
    abscissa       = abscissa, 
    lambdas        = Lambdas_Roughness, 
    basisBreakFreq = basisBreakFreq,
    basisOrder     = basisOrder,
    pbasis         = NULL
  )$curvesSmooth
  
}