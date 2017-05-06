#'@import fdakma
#'@export

########################################################################## 
## To use kma.similarity for a matrix of curves with a template curve 
## ( median or mean )
########################################################################## 
getKmaSimilarityWithTemplate <- function( 
  Mat, 
  Template,
  Xaxis,
  similarity.method = c( "d1.pearson", "d1.L2", "d1.L2.centered" ),
  Deriv = FALSE
){
  ## Assume that similarity measure is to be estimated between columns
  ## of the matrix and the Template
  ## Assumed that the x-axis of all the functions will be the same
  ## Return the average of the similarity index, between all the columns
  ## d1 corresponds to first derivatives. Make sure the matrix Mat is that of
  ## first derivatives and not just the function values
  ## If data is not derivatives, then change to deriv
  
  if( Deriv == FALSE ){
    dx <- diff( Xaxis )[1]
    Curves.D1 <- round( diff( Mat )/dx, 4 )
    Template.D1 <- round( diff( Template )/dx, 4 )
    Xaxis <- Xaxis[-1]
  } else{
    Curves.D1 <- Mat
    Template.D1 <- Template
  }
  
  Sim <- c( )
  options( warn = -1 )
  for( Col in 1:ncol( Curves.D1 ) ){
    Sim <- c( Sim, kma.similarity( 
      x.f                = Xaxis, 
      y1.f               = Curves.D1[,Col], 
      x.g                = Xaxis, 
      y1.g               = Template.D1, 
      similarity.method  = similarity.method ) )
  }
  Sim <- round( Sim, 4 )
  options( warn = 0 )
  return( Sim )
}
########################################################################## 
