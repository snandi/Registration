#'@import fda
#'@import ggplot2
#'@export

############################################################################### 
## This function conducts iteration registration to weighted average template
## estimated at every iteration
############################################################################### 
registerIterated <- function( 
  dataToRegister, 
  Lambdas_ConstrainedWarping = c(0.001, 0.0001, 0.00005), 
  LambdaDefault = 0.00001,
  abscissaFrom,
  abscissaTo,
  abscissaIncrement,
  basisOrder = 5,
  basisBreakFreq = 3,
  Lambdas_Roughness = exp( -5:5 ),
  outlierTrimPct = 0.20,
  RE_REGISTER = FALSE,            ## Re-register registered curves from previous iteration?
  SimMeanDiff_Threshold = 0.0001,
  MinSimilarityThreshold = 0.15,
  MAX_ITERATION = 6
){
  abscissa <- seq( from = abscissaFrom, to = abscissaTo, by = abscissaIncrement ) 
  moleculesForConsensus <- colnames( dataToRegister )
  maxOutliersPerIteration <- round( outlierTrimPct*length( moleculesForConsensus ) )
  
  dataToRegisterSmooth <- createSmoothFD( 
    curvesToSmooth = dataToRegister, 
    abscissa       = abscissa, 
    lambdas        = Lambdas_Roughness, 
    basisBreakFreq = basisBreakFreq,
    basisOrder     = basisOrder,
    pbasis         = NULL
  )$curvesSmooth
  
  Index <- 0
  Iterate <- TRUE
  
  Index <- Index + 1
  Sim_toTemplate <- c() ## Contains the similarities of all curves to their median
  
  while ( Iterate == TRUE ){
    print( paste( 'Starting Iteration', Index ) )
    if( !( is.na( Lambdas_ConstrainedWarping[ Index ] ) ) ){
      Lambda <- Lambdas_ConstrainedWarping[ Index ]
    } else {
      Lambda <- LambdaDefault
    }
    
    basisLength <- length( abscissa )/basisBreakFreq
    basisBreaks <- seq( from = abscissaFrom, to = abscissaTo, length.out  = basisLength )
    basisObj <- create.bspline.basis( rangeval  = range( abscissa ), norder = basisOrder, breaks  = basisBreaks )
    Wfd0 <- fd( matrix( data = 0, nrow = basisObj$nbasis, ncol = 1 ), basisObj )
    WfdParobj <- fdPar( fdobj = Wfd0, Lfdobj = 3, lambda = Lambda )
    
    if( Index == 1 ){
      ## Detect Outliers method "RProj"
      Outliers_RProj_Trim <- getFunctionalOutliers (
        Curves      = dataToRegister, 
        Xaxis       = abscissa, 
        Names       = list( main = 'Main', xlab = 'PixelPosition', ylab = 'Intensity' ),
        DepthType   = 'RProj',
        N_Bootstrap = 500,
        Trim        = 'Yes',
        TrimPct     = outlierTrimPct
      )
      print( Outliers_RProj_Trim$outliers )
      
      ## Detect Outliers method "FM"
      Outliers_FM_Trim <- getFunctionalOutliers (
        Curves = dataToRegister, 
        Xaxis = abscissa, 
        Names = list( main = 'Main', xlab = 'PixelPosition', ylab = 'Intensity' ),
        DepthType = 'FM',
        N_Bootstrap = 500,
        Trim = 'Yes',
        TrimPct = outlierTrimPct
      )
      print( Outliers_FM_Trim$outliers )
      Outliers_Union <- union( x = Outliers_RProj_Trim$outliers, y = Outliers_FM_Trim$outliers )    
      
      ## Detect Outliers method "RTukey"
      Outliers_RTukey_Trim <- getFunctionalOutliers (
        Curves = dataToRegister, 
        Xaxis = abscissa, 
        Names = list( main = 'Main', xlab = 'PixelPosition', ylab = 'Intensity' ),
        DepthType = 'RTukey',
        N_Bootstrap = 500,
        Trim = 'Yes',
        TrimPct = outlierTrimPct
      )
      print( Outliers_RTukey_Trim$outliers )
      Outliers_Union <- union( x = Outliers_Union, y = Outliers_RTukey_Trim$outliers )    
      
      ## Detect Outliers method "Mode"
      Outliers_Mode_Trim <- getFunctionalOutliers (
        Curves = dataToRegister, 
        Xaxis = abscissa, 
        Names = list( main = 'Main', xlab = 'PixelPosition', ylab = 'Intensity' ),
        DepthType = 'RTukey',
        N_Bootstrap = 500,
        Trim = 'Yes',
        TrimPct = outlierTrimPct
      )
      print( Outliers_Mode_Trim$outliers )
      Outliers_Union <- union( x = Outliers_Union, y = Outliers_Mode_Trim$outliers )    
      
      Keep <- moleculesForConsensus %w/o% Outliers_Union
      dataToRegisterNoOutliers <- dataToRegister[ , Keep ]
      
      ## Estimate weighted mean, without the outliers (TEMPLATE)
      Consensus <- getWeightedConsensus( Curves = dataToRegisterNoOutliers, abscissa = abscissa )
      templateToRegister <- Consensus[['Mean_Weighted']]
      templateToRegisterFD <- createSmoothFD(
        curvesToSmooth = templateToRegister, 
        abscissa       = abscissa, 
        lambdas        = Lambdas_Roughness, 
        basisBreakFreq = basisBreakFreq, 
        basisOrder     = basisOrder, 
        pbasis         = NULL
      )
      
      ## Estimate Similarity with Template before registration
      options( warn = -1 )
      Sim_Before_Regist <- getKmaSimilarityWithTemplate(
        Mat               = dataToRegister, 
        Template          = templateToRegister,
        Xaxis             = abscissa,
        similarity.method = c( "d1.pearson" ),
        Deriv             = FALSE
      )
      options( warn = 0 )
      Sim_toTemplate <- round( Sim_Before_Regist, 4 )
      names( Sim_toTemplate ) <- colnames( dataToRegister )
      print( Sim_toTemplate )
      
      ## Register first iteration
      Regfd_All <- register.fd(
        y0fd      = templateToRegisterFD$fd, 
        yfd       = dataToRegisterSmooth$fd, 
        WfdParobj = WfdParobj, 
        dbglev    = 0,
        periodic  = FALSE, 
        crit      = 2
      ) 
      
      Regfd1 <- Regfd_All
      Regfd1_eval <- eval.fd( evalarg = abscissa, fdobj = Regfd1$regfd )
      colnames( Regfd1_eval ) <- moleculesForConsensus
      
      ## Estimate Similarity with Template after registration      
      options( warn = -1 )
      Sim_After_Regist <- getKmaSimilarityWithTemplate(
        Mat               = Regfd1_eval, 
        Template          = templateToRegister,
        Xaxis             = abscissa,
        similarity.method = "d1.pearson",
        Deriv             = FALSE
      )
      options( warn = 0 )
      Sim_toTemplate <- rbind( Sim_toTemplate, round( Sim_After_Regist, 4 ) )
      names( Sim_After_Regist ) <- moleculesForConsensus
      print( paste( 'Iteration', Index, 'Complete' ) )
      
      print( Sim_toTemplate )
      ## End of if statement for first iteration, started in line 52 
    } else{ 
      rm( Regfd_All )
      
      ## Detect Outliers
      Outliers_RProj_Trim <- getFunctionalOutliers (
        Curves      = Regfd1_eval, 
        Xaxis       = abscissa, 
        Names       = list( main = 'Main', xlab = 'PixelPosition', ylab = 'Intensity' ),
        DepthType   = 'RProj',
        N_Bootstrap = 500,
        Trim        = 'Yes',
        TrimPct     = outlierTrimPct
      )
      print( Outliers_RProj_Trim$outliers )
      
      Outliers_FM_Trim <- getFunctionalOutliers (
        Curves      = Regfd1_eval, 
        Xaxis       = abscissa, 
        Names       = list( main = 'Main', xlab = 'PixelPosition', ylab = 'Intensity' ),
        DepthType   = 'FM',
        N_Bootstrap = 500,
        Trim        = 'Yes',
        TrimPct     = outlierTrimPct
      )
      print( Outliers_FM_Trim$outliers )
      
      Outliers_Union <- unique( c( Outliers_RProj_Trim$outliers, Outliers_FM_Trim$outliers, 
                                   names( Sim_After_Regist )[ Sim_After_Regist < MinSimilarityThreshold ] ) )    
      print( Outliers_Union )
      if( length( Outliers_Union ) > maxOutliersPerIteration ){
        OrderedSim_After_Regist <- Sim_After_Regist[ Outliers_Union ][order( Sim_After_Regist[ Outliers_Union ] )]
        Outliers_Union <- names( OrderedSim_After_Regist[1:maxOutliersPerIteration] )
      }
      print( Outliers_Union )
      
      Keep <- moleculesForConsensus %w/o% Outliers_Union
      Regfd1_noOutliers <- Regfd1_eval[ , Keep ]
      
      ## Estimate weighted mean, without the outliers (TEMPLATE)
      Consensus <- getWeightedConsensus( Curves = Regfd1_noOutliers, abscissa = abscissa )
      templateToRegister <- Consensus[['Mean_Weighted']]
      templateToRegisterFD <- createSmoothFD(
        curvesToSmooth = templateToRegister, 
        abscissa       = abscissa, 
        lambdas        = Lambdas_Roughness, 
        basisBreakFreq = basisBreakFreq, 
        basisOrder     = basisOrder, 
        pbasis         = NULL
      )
      
      ## Register subsequent iteration
      if( RE_REGISTER ){
        Regfd_All <- register.fd(
          y0fd      = templateToRegisterFD$fd, 
          yfd       = Regfd1$regfd, 
          WfdParobj = WfdParobj, 
          dbglev    = 0,
          periodic  = FALSE, 
          crit      = 2
        ) 
      } else{
        Regfd_All <- register.fd(
          y0fd      = templateToRegisterFD$fd, 
          yfd       = dataToRegisterSmooth$fd, 
          WfdParobj = WfdParobj, 
          dbglev    = 0,
          periodic  = FALSE, 
          crit      = 2
        ) 
      }
      
      Regfd1 <- Regfd_All
      Regfd1_eval <- eval.fd( evalarg = abscissa, fdobj = Regfd1$regfd )
      colnames( Regfd1_eval ) <- moleculesForConsensus
      
      ## Estimate Similarity with Template after registration      
      options( warn = -1 )
      Sim_After_Regist <- getKmaSimilarityWithTemplate(
        Mat               = Regfd1_eval, 
        Template          = templateToRegister,
        Xaxis             = abscissa,
        similarity.method = "d1.pearson",
        Deriv             = FALSE
      )
      options( warn = 0 )
      Sim_toTemplate <- rbind( Sim_toTemplate, round( Sim_After_Regist, 4 ) )
      names( Sim_After_Regist ) <- moleculesForConsensus
      
      Outliers_Union <- union( 
        x = union( Outliers_RProj_Trim$outliers, Outliers_FM_Trim$outliers ), 
        y = names( Sim_After_Regist )[ Sim_After_Regist < MinSimilarityThreshold ] 
      )
      if( length( Outliers_Union ) > maxOutliersPerIteration ){
        OrderedSim_After_Regist <- Sim_After_Regist[ Outliers_Union ][order( Sim_After_Regist[ Outliers_Union ] )]
        Outliers_Union <- names( OrderedSim_After_Regist[1:maxOutliersPerIteration] )
      }
      print( Outliers_Union )
      
      print( paste( 'Iteration', Index, 'Complete' ) )
      print( Sim_toTemplate )
    } ## End of iterated registration if condition
    
    ## Check if iteration should continue
    notOutlier <- moleculesForConsensus %w/o% Outliers_Union
    
    SimPrevious <- Sim_toTemplate[Index, ]
    
    SimMean_New <- mean( Sim_After_Regist[ notOutlier ] )
    SimMean_Old <- mean( SimPrevious[ notOutlier ] )
    print( paste( 'SimMeans', round( SimMean_New, 5 ), round( SimMean_Old, 5 ) ) )
    
    SimMeanDiff <- abs( SimMean_New - SimMean_Old )
    print( SimMeanDiff )
    
    if( SimMeanDiff < SimMeanDiff_Threshold ) Iterate <- FALSE
    Index <- Index + 1
    if( Index == 2) Iterate <- TRUE   ## This line ensures at least 2 iterations.
    if( Index == MAX_ITERATION) Iterate <- FALSE  ## Maximum 6 iterations
    
  } ## End of while loop started in line 38
  
  N_Iter <- Index - 1
  Regfd_Final <- Regfd_All
  registeredCurves <- eval.fd( evalarg = abscissa, fdobj = Regfd_Final$regfd )
  registeredCurves.D1 <- eval.fd( evalarg = abscissa, fdobj = Regfd_Final$regfd, Lfdobj = 1 )
  SimPrevious <- Sim_toTemplate[ ( Index - 1 ), ]
  
  colnames( registeredCurves ) <- moleculesForConsensus
  colnames( registeredCurves.D1 ) <- moleculesForConsensus
  
  registeredCurvesAll <- registeredCurves
  registeredCurvesAll.D1 <- registeredCurves.D1
  
  rm( registeredCurves, registeredCurves.D1 )
  ## Detect Outliers
  Outliers_RProj_Trim <- getFunctionalOutliers (
    Curves      = registeredCurvesAll, 
    Xaxis       = abscissa, 
    Names       = list( main = 'Main', xlab = 'PixelPosition', ylab = 'Intensity' ),
    DepthType   = 'RProj',
    N_Bootstrap = 500,
    Trim        = 'Yes',
    TrimPct     = outlierTrimPct
  )
  print( Outliers_RProj_Trim$outliers )
  
  Outliers_FM_Trim <- getFunctionalOutliers (
    Curves      = registeredCurvesAll, 
    Xaxis       = abscissa, 
    Names       = list( main = 'Main', xlab = 'PixelPosition', ylab = 'Intensity' ),
    DepthType   = 'FM',
    N_Bootstrap = 500,
    Trim        = 'Yes',
    TrimPct     = outlierTrimPct
  )
  print( Outliers_FM_Trim$outliers )
  
  Outliers_Union <- unique( c( Outliers_RProj_Trim$outliers, Outliers_FM_Trim$outliers, 
                               names( Sim_After_Regist )[ Sim_After_Regist < MinSimilarityThreshold ] ) )    
  
  if( length( Outliers_Union ) > maxOutliersPerIteration ){
    OrderedSim_After_Regist <- Sim_After_Regist[ Outliers_Union ][order( Sim_After_Regist[ Outliers_Union ] )]
    Outliers_Union <- names( OrderedSim_After_Regist[1:maxOutliersPerIteration] )
  }
  
  ## This returns only the curves that should be used to estimate the consensus
  notOutlier <- moleculesForConsensus %w/o% Outliers_Union
  print( c("not outliers", notOutlier ) )
  aboveThreshold <- moleculesForConsensus[ SimPrevious > MinSimilarityThreshold ]
  print( c( "above threshold", aboveThreshold ) )
  curvesForConsensus <- intersect( notOutlier, aboveThreshold )
  print( c( "consensus curves", curvesForConsensus ) )
  
  registeredCurves <- registeredCurvesAll[ , curvesForConsensus ]
  registeredCurves.D1 <- registeredCurvesAll.D1[ , curvesForConsensus ]
  
  ## Estimate weighted mean, without the outliers
  Consensus <- getWeightedConsensus( Curves = registeredCurves, abscissa = abscissa )
  FinalConsensus <- Consensus[['Mean_Weighted']]
  
  ## Final list to return
  return( list( dataToRegister         = dataToRegister, 
                Regfd_Final            = Regfd_Final,
                FinalConsensus         = FinalConsensus,
                curvesForConsensus     = curvesForConsensus,
                registeredCurvesAll    = registeredCurvesAll, 
                registeredCurvesAll.D1 = registeredCurvesAll.D1, 
                registeredCurves       = registeredCurves, 
                registeredCurves.D1    = registeredCurves.D1, 
                Sim_toTemplate         = Sim_toTemplate
  ) )
}
