#' @title A dataframe of counts of k-mers
#' 
#' @description Return a dataframe of counts of k-mers of a user provided genomic sequence
#' 
#' @param seqData Sequence data of class \code{SeqFastadna}
#' @param maxNmers An integer between 1 and 5
#' @param seqWindow Defaults to 206. This is the subsequence window in which the counts of different 
#' k-mers is returned
#' 
#' @return A dataframe of counts of k-mers. Each column corresponds to a k-mer, and each row corresponds
#' to a 206 \code{seqWindow}
#' 
#' @author Subhrangshu Nandi, PhD Statistics, UW Madison; snandi@wisc.edu or nands31@gmail.com
#' 
#' @examples 
#' fpath <- system.file( "extdata", "chr22_Segment750.fa", package = "FscanStats" )
#' seqData <- seqinr::read.fasta( file = fpath )[[1]]
#' segmentSeqCompData <- getSeqCompCount( 
#'    seqData   = seqData,
#'    maxNmers  = 4
#' )
#' 
#' @importFrom seqinr count
#' @export

getSeqCompCount <- function( seqData, maxNmers = 5, seqWindow = 206 ){

  lengthSequence <- length( seqData )
  N_SubSegments <- lengthSequence / seqWindow

  SeqBy206bp <- getSubsequenceChunks( 
    Sequence               = seqData, 
    numChunks              = N_SubSegments, 
    force.number.of.groups = TRUE 
  )
  
  # SeqCompData <- segmentPixels  
  for( i in 1:maxNmers ){
    SplitSeqCount <- as.data.frame( do.call( what = rbind, lapply( X = SeqBy206bp, FUN = seqinr::count, wordsize = i ) ) )
    if( i == 1 ) { SeqCompData <- SplitSeqCount } else { SeqCompData <- cbind( SeqCompData, SplitSeqCount ) }
    rm( SplitSeqCount )
  }
  return( SeqCompData )
}
