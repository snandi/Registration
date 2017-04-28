#'@export
prepareDataWithSymmetricKernel <- function( 
  SeqCompData, 
  kernelType = c( 'Uniform', 'Gaussian' ), 
  colnamesNonSeq = c( 'consensusMean', 'Interval', 'subInterval', 'subIntervalPixels' )
){
  if( kernelType == 'Gaussian' ){
    K0 <- KernelGaussian( h = 0 )
    K1 <- KernelGaussian( h = 1 )
    K2 <- KernelGaussian( h = 2 )
  } else{
    K0 <- K1 <- K2 <- 1
  }
  
  # colnamesNonSeq <- c( 'consensusMean', 'Interval', 'subInterval', 'subIntervalPixels' )
  colnamesSeqComp <- colnames( SeqCompData ) %w/o% colnamesNonSeq
  
  fromRow <- 3
  toRow <- nrow(SeqCompData) - 2
  
  dataSubInt <- SeqCompData[ fromRow:toRow, colnamesNonSeq ]
  
  dataK0 <- SeqCompData[ fromRow:toRow, colnamesSeqComp ]
  
  dataKp1 <- K1 * SeqCompData[ (fromRow + 1):(toRow + 1), colnamesSeqComp ]
  dataKm1 <- K1 * SeqCompData[ (fromRow - 1):(toRow - 1), colnamesSeqComp ]
  dataK1 <- dataKp1 + dataKm1
  colnames( dataK1 ) <- paste( colnames( dataK0 ), 'p1', sep = '_' )
  
  dataKp2 <- K2 * SeqCompData[ (fromRow + 2):(toRow + 2), colnamesSeqComp ]
  dataKm2 <- K2 * SeqCompData[ (fromRow - 2):(toRow - 2), colnamesSeqComp ]
  dataK2 <- dataKp2 + dataKm2
  colnames( dataK2 ) <- paste( colnames( dataK0 ), 'p2', sep = '_' )
  
  dataWithKernel <- cbind( dataSubInt, dataK0, dataK1, dataK2 )
  # dataWithKernel <- cbind( consensusMean = SeqCompData$consensusMean[fromRow:toRow], dataWithKernel )

  return( dataWithKernel )  
}
