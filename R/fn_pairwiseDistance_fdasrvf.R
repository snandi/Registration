fn_pairwiseDistance_fdasrvf <- function(Mat, Xaxis){
  require(gtools) 
  require(fdasrv)
  ## Assume that similarity measure is to be estimated between columns
  ## of the matrix. 
  ## Assumed that the x-axis of all the functions will be the same
  Vector <- c(1:ncol(Mat))
  ColCombs <- as.data.frame(combinations(n=length(Vector), r=2, v=Vector))
  ColCombs$Dx <- 0 ## Phase distance
  ColCombs$Dy <- 0 ## Amplitude distance
  Row <- 1
  for(Row in 1:nrow(ColCombs)){
    Distance <- elastic.distance(f1=Mat[,ColCombs[Row,1]], f2=Mat[,ColCombs[Row,2]], time=Xaxis)
    ColCombs[Row,'Dx'] <- Distance$Dx
    ColCombs[Row,'Dy'] <- Distance$Dy
  }
  return(ColCombs)
}
################################################################## 
