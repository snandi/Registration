###########################################################################################################
# Gina-Maria Pomann
# 9/11/2012
#This script will test along each dimension of of the input data sets to see if they are different
#using the Anderson Darling Test
###########################################################################################################
#require( adk )

#'@import adk 
#'@export

test_AD_2sample <- function( datX, datY, p ){ 
  #p the dimension of the multivariate data
  #nmeansY the number of data sets in the Y array to be tested against 
  #those in the X array 
  
  #no transformation, just test 
  U <- datX
  V <- datY
  if( is.vector( U ) == T || is.vector( V ) == T ){ p = 1 }
  if( is.vector( U ) == F && is.vector( V ) == F ){ p = min( dim( U )[2],dim( V )[2] ) }
  #p=min( dim( U )[2],dim( V )[2] )   
  pvalvec <- rep( NA,p )
  
  
  if( p == 1 ){ k=1; pvalvec[k] <- adk.test( U,V )$adk[1,2] }
  if( p > 1 ){ 
    for( k in 1:p ){ 
      #get the p value from the anderson darling test
      pvalvec[k] <- adk.test( U[,k],V[,k] )$adk[1,2] 
    }
  }
  #now return the pvalue vector
  list( pval = pvalvec )
}

x <- list(c(1,3,2,5,7),c(2,8,1,6,9,4),c(12,5,7,9,11))
out <- adk.test(x) 
# or out <- adk.test(c(1,3,2,5,7),c(2,8,1,6,9,4), c(12,5,7,9,11))
## Examine the component names of out
names(out)

## Examine the matrix adk of out.
out$adk
