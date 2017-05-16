#' @importFrom seqinr s2c
#' @export
#' 
########################################################################
## This function takes in a DNA sequence "agggctattaaccctttaaa" and returns    
## smaller subsequence chunks based on user specified length. This    
## was adapted from a function chunks.2 found online
## This was fn_subseq in Project_CurveReg
########################################################################
getSubsequenceChunks <- function( 
  Sequence, 
  numChunks, 
  force.number.of.groups = TRUE
) { 
  ## The argument force.number.of.groups=TRUE means it will append the 
  ## remainder of the elements to the last chunk
  
  ## Sequence should be a list of characters instead of a long string
  if( length( Sequence ) == 1 ){
    Sequence <- seqinr::s2c( Sequence )
  }
  
  len <- length( Sequence )
  groups <- trunc( len/numChunks )
  overflow <- len%%numChunks 

  if( force.number.of.groups ) {
    f1 <- as.character( sort( rep( 1:numChunks, groups ) ) )
    f <- as.character( c( f1, rep( numChunks, overflow ) ) )
  } else {
    f1 <- as.character( sort( rep( 1:groups, numChunks ) ) )
    f <- as.character( c( f1, rep( "overflow", overflow ) ) )
  }
  
  g <- split( Sequence, f )
  
  if( force.number.of.groups ) {
    g.names <- names( g )
    g.names.ordered <- as.character( sort( as.numeric( g.names ) ) )
  } else {
    g.names <- names( g[-length( g )] )
    g.names.ordered <- as.character( sort( as.numeric( g.names ) ) )
    g.names.ordered <- c( g.names.ordered, "overflow" )
  }
  
  return( g[ g.names.ordered ] )
}
########################################################################
#getSubsequenceChunks(Sequence = s2c( "agggctattaaccctttaaa" ), numChunks = 3, force.number.of.groups = T )
