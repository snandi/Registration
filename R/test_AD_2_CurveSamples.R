## This function returns: test statistic values, p-values of Anderson-Darling
## goodness of fit test on the FPCA scores of 2 samples of curves

## This uses the test_AD_2sample function

## Reference: Pomann, G.M., Staicu, A.M., and Ghosh,S. ( 2016 ). A two-sample distribution-free test 
## for functional data with application to a diffusion tensor imaging study of multiple sclerosis. 
## Journal of the Royal Statistical Society: Series C ( Applied Statistics ).

#' @import refund
#' @import MASS

#' @export

test_AD_2_CurveSamples <- function( Mat1, Mat2, varpct = 0.95 ){ 

  ## This functions requires two matrices of same number of rows, 
  ## each column being a separate curve; 

  if( nrow( Mat1 ) != nrow( Mat2 ) ){
    stop( "Both datasets should have same numbers of observations per curve" )
  }
  L <- nrow( Mat1 )
  n1 <- ncol(Mat1)
  n2 <- ncol(Mat2)  

  Data_Combined <- cbind(Mat1, Mat2)

  FPCA.out <- refund::fpca.sc(Y = t(Data_Combined), var = TRUE, pve = varpct)
  
  scores1 <- FPCA.out$scores[1:n1,]
  scores2 <- FPCA.out$scores[(n1 + 1):(n1 + n2),]  
  ncomp <- dim(FPCA.out$efunctions)[2]
  
  testOutput <- test_AD_2sample( datX = scores1, datY = scores2, p = ncomp)
  return(testOutput)
	
}
  
# require(fda)
# require(mgcv)
# require(adk)
# require(refund)
# data(growth)
# Mat1 <- growth[['hgtm']]
# Mat2 <- growth[['hgtf']]
# testOutput <- test_AD_2_CurveSamples(Mat1 = Mat1, Mat2 = Mat2, varpct = 0.95)

