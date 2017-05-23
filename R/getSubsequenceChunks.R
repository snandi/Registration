#' @title Get subsequence chunks from a long sequence
#' 
#' @description This function takes in a DNA sequence "agggctattaaccctttaaa" in the form of a list of 
#' characters and returns smaller subsequence chunks based on user specified number.
#' 
#' @param Sequence A DNA sequence in the form of a list of characters: "a" "g" "g" "c" "t", similar 
#' to the output of \code{seqinr::s2c}
#' @param numChunks Number of subsequence chunks desired
#' @param force.number.of.groups Defaults to \code{TRUE}. If \code{TRUE}, then it will append the 
#' remainder of the elements to the last chunk
#' 
#' @author Subhrangshu Nandi; PhD, Statistics; snandi@wisc.edu or nands31@gmail.com
#' 
#' @return Returns a list of subsequence chunks.
#' 
#' @examples
#' getSubsequenceChunks( 
#'   Sequence = "agggctattaaccctttaaa", 
#'   numChunks = 3, 
#'   force.number.of.groups = TRUE )
#'
#' @keywords DNA-sequence subsequence
#'  
#' @importFrom seqinr s2c
#' @export
#' 

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
