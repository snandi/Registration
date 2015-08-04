####################################################################################
## This function library has simple R functions written by Subhrangshu Nandi
####################################################################################

################################################################## 
## Returns the SE of mean of each row of a dataset
##################################################################
rowSE <- function(Data){
  SE <- apply(X = Data, MARGIN = 1, FUN=function(Row){sd(Row)/sqrt(length(Row))})
  return(SE)
}
################################################################## 

################################################################## 
## Returns the SD of each row of a dataset
##################################################################
rowSD <- function(Data){
  SD <- apply(X = Data, MARGIN = 1, FUN=function(Row){sd(Row)})
  return(SD)
}

colSD <- function(Data){
  rowSD(t(Data))
}
################################################################## 

################################################################## 
## Returns the SD of each row of a dataset
##################################################################
rowVar <- function(Data){
  SD <- apply(X = Data, MARGIN = 1, FUN=function(Row){var(Row)})
  return(SD)
}

colVar <- function(Data){
  rowVar(t(Data))
}
################################################################## 

######################### Convert NAs to Zero ##########################
na.is.zero <- function(X)
{
  X1 <- X
  X1[is.na(X)] <- 0.0
  return(X1)
}
########################################################################

########################################################################
"%notin%" <- function(x, y){
  if(x %in% y){
    return(FALSE)
  } else{
    return(TRUE)
  }
}
########################################################################

########################################################################
"%w/o%" <- function(x, y){
  return(x[!x %in% y])
}
########################################################################

########################################################################
Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}
########################################################################

########################################################################
## Returns if x is inside range or not
########################################################################
inRange <- function(x, Range=NULL, Vec=NULL){
  if (is.null(Range) & is.null(Vec)) stop('Need either Vec or Range')
  if (!is.null(Vec)) Range <- range(Vec)
  if (Range[1] == Range[2]) stop('Range of vector ill defined')
  if(x >= Range[1] & x <= Range[2]) {
    Ans <- TRUE
  } else {
    Ans <- FALSE
  }
  return(Ans)
}
########################################################################

########################################################################
## Returns sup of T-type statistic between two matrices
########################################################################
fn_absTStat <- function(Mat1, Mat2){
  Mat1 <- as.matrix(Mat1)
  Mat2 <- as.matrix(Mat2)

  N1 <- ncol(Mat1)
  N2 <- ncol(Mat2)
  
  mean1 <- rowMeans(Mat1)
  mean2 <- rowMeans(Mat2)
  var1 <- apply(X=Mat1, MARGIN=1, FUN=var)/N1
  var2 <- apply(X=Mat2, MARGIN=1, FUN=var)/N2
  
  Tvals <- abs(mean1 - mean2)/sqrt(var1 + var2)
  TSup <- max(Tvals)
  TMedian <- median(Tvals)
  TMean <- mean(Tvals)
  
  return(list(TSup = TSup, TMedian = TMedian, TMean = TMean))
}
