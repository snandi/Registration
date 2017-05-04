#'@import fda
#'@import ggplot2
#'@export

############################################################################### 
## This function conducts single iteration registration, according to the penalty 
## parameter lambda passed as argument
############################################################################### 
registerSingleIter <- function( 
  dataToRegister, 
  Lambda_ConstrainedWarping, 
  abscissaFrom,
  abscissaTo,
  abscissaIncrement,
  Convergence_Threshold = 1e-05,
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

  if( !( is.na( Lambda_ConstrainedWarping ) ) ){
    Lambda <- Lambda_ConstrainedWarping
  } else {
    Lambda <- 0.005
  }
  
  ## Prepare basis for warping function  
  basisBreaks <- seq( 
    from        = abscissaFrom, 
    to          = abscissaTo, 
    length.out  = length( abscissa )/basisBreakFreq )
  
  basisObj <- create.bspline.basis( 
    rangeval  = range( abscissa ), 
    norder    = basisOrder, 
    breaks    = basisBreaks )
  
  Wfd0 <- fd( matrix( data = 0, nrow = basisObj$nbasis, ncol = 1 ), basisObj )
  WfdParobj <- fdPar( fdobj = Wfd0, Lfdobj = 2, lambda = Lambda )
  
  ## Prepare template to register to
  templateToRegister <- rowMeans( dataToRegister )
  
  templateToRegister_fd <- createSmoothFD( 
    curvesToSmooth = templateToRegister, 
    abscissa       = abscissa, 
    lambdas        = Lambdas_Roughness, 
    basisBreakFreq = basisBreakFreq,
    basisOrder     = basisOrder,
    pbasis         = NULL
  )$curvesSmooth
  
  ## REGISTER
  Regfd_All <- fda::register.fd( 
    y0fd      = templateToRegister_fd$fd, 
    yfd       = dataToRegisterSmooth$fd, 
    WfdParobj = WfdParobj, 
    dbglev    = 0,
    conv      = Convergence_Threshold,
    periodic  = FALSE, 
    crit      = 2 ) 
  
  Regfd_Final <- Regfd_All
  registeredCurves <-   eval.fd( evalarg = abscissa, fdobj = Regfd_Final$regfd )
  registeredCurves.D1 <-   eval.fd( evalarg = abscissa, fdobj = Regfd_Final$regfd, Lfdobj = 1 )
  
  return( list( 
    dataToRegister      = dataToRegister, 
    Regfd_Final         = Regfd_Final,
    registeredCurves    = registeredCurves, 
    registeredCurves.D1 = registeredCurves.D1
  ) ) 
}
################################################################################ 
