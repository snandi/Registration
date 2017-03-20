#'@import fdakma
#'@import gtools
#'@export

########################################################################## 
## To use kma.similarity for a matrix                           
########################################################################## 
getKmaSimilarityMatrix <- function( 
  Mat, 
  Xaxis,
  similarity.method = c( "d1.pearson", "d1.L2", "d1.L2.centered" )
 ){
  ## Assume that similarity measure is to be estimated between columns
  ## of the matrix. 
  ## Assumed that the x-axis of all the functions will be the same
  ## Return the average of the similarity index, between all the columns
  ## d1 corresponds to first derivatives. Make sure the matrix Mat is that of
  ## first derivatives and not just the function values
  
  Vector <- c( 1:ncol( Mat ) )
  if( !is.null( colnames( Mat ) ) ){
    Vector <- colnames( Mat )[Vector]
  }
  ColCombs <- as.data.frame( gtools::combinations( n = length( Vector ), r = 2, v = Vector ) )
  ColCombs$Similarity <- 0
  ColCombs$V1 <- as.vector( ColCombs$V1 )
  ColCombs$V2 <- as.vector( ColCombs$V2 )
  
  for( Row in 1:nrow( ColCombs ) ){
    ColCombs[Row,'Similarity'] <- fdakma::kma.similarity ( 
      x.f               = Xaxis, 
      y1.f              = Mat[,ColCombs[Row,1]], 
      x.g               = Xaxis, 
      y1.g              = Mat[,ColCombs[Row,2]], 
      similarity.method = similarity.method )
  }
  return( list( ColCombs = ColCombs, AvgSimilarity = mean( ColCombs$Similarity ) ) )
}
########################################################################## 
