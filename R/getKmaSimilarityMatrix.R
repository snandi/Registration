#' @import fdakma
#' @import gtools
#' @export

#' @title Get a pairwise similarity matrix of a group of curves
#' 
#' @description Returns pairwise similarity between columns of a matrix of first derivatives of curves, 
#' and the average of the similarities
#'
#' @param Mat Matrix with each column representing a separate curve
#' 
#' @param Xaxis Vector of numbers, the abscissa
#' 
#' @param similarity.method Type of similarity method to use
#' 
#' @details Assumed that the x-axis of all the functions will be the same 
#' d1 corresponds to first derivatives. Make sure the matrix Mat is that of
#' first derivatives and not just the function values
#' 
#' @author Subhrangshu Nandi; Statistics PhD student, UW Madison; snandi@wisc.edu or nands31@gmail.com
#' 
#' @usage getKmaSimilarityMatrix( Mat, Xaxis, similarity.method = c("d1.pearson", "d1.L2", "d1.L2.centered") )
#' 
#' @examples
#' data( growth, package = 'fda')
#' Mat1 <- growth[['hgtm']]
#' Arguments <- growth[['age']]
#' getKmaSimilarityMatrix( Mat = Mat1, Xaxis = Arguments, similarity.method = "d1.pearson" )
#' 
#' @return List of two objects
#' \item{ColCombs}{A data frame of pairwise comparison similarity index}
#' \item{AvgSimilarity}{Average of all similarities}
#' 

########################################################################## 
## To use kma.similarity for a matrix                           
########################################################################## 
getKmaSimilarityMatrix <- function( 
  Mat, 
  Xaxis,
  similarity.method = c( "d1.pearson", "d1.L2", "d1.L2.centered" )
 ){
  
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
