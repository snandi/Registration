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
  pvaluec <- rep( x = NA, times = p )
  testStat <- rep( x = NA, times = p )
  
  if( p == 1 ){ 
    k = 1
    testOutput <- adk.test( U, V )
    pvaluec[k] <- testOutput$adk[1,2]
    testStat[k] <- testOutput$adk[1,1]
  }
  if( p > 1 ){ 
    for( k in 1:p ){ 
      #get the p value from the anderson darling test
      testOutput <- adk.test( U[,k], V[,k] )
      pvaluec[k] <- testOutput$adk[1,2] 
      testStat[k] <- testOutput$adk[1,1]
      rm( testOutput )
    }
  }
  #now return the pvalue vector and the test statistic vector
  return(list( pval = pvaluec, testStat = testStat ))
}

# U1 <- rnorm(n = 20, mean = 2, sd = 1)
# U2 <- rnorm(n = 20, mean = 1, sd = 2)
# 
# V1 <- rnorm(n = 20, mean = 0, sd = 1)
# V2 <- rnorm(n = 20, mean = 1, sd = 2)
# 
# Mat1 <- cbind(U1, U2)
# Mat2 <- cbind(V1, V2)
# testOutput <- test_AD_2sample( datX = Mat1, datY = Mat2, p = 2 )

  