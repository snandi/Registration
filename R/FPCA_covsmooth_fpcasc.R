############################################################################
## This was downloaded from Ana-Maria Staicu's webpage 
## http://www4.stat.ncsu.edu/~staicu/Research.html
## This is the R code for their paper: A Two Sample Distribution-Free Test 
## for Functional Data with Application to a Diffusion Tensor Imaging Study 
## of Multiple Sclerosis
############################################################################

#Gina-Maria Pomann
#9/29/12

#this code smooths the raw covariances first and then pools them. It uses the smoothing methodology 
#from the refund package. 

#Inputs:
#Y.1  : Data set 1-needs to be read in as nxd where n is the number of curves
#and d is the number of observations along the curve
# Each row is a separate curve
#Y.2 : data set 2-needs to be read in as nxd where n is the number of curves
#and d is the number of observations along the curve

#threshold (default .99): The proporton of variance explained by the eigenvectors
#that are used to explain the data and create the smoothed covaraince matrix. 

#Outputs:
#cov.dat.1 : smoothed covariance of data set 1
#cov.dat.2 : smoothed covariance of data set 2
#L = L_U : The number of eigenfunctions that will be used 
#to explain the data (this depends on the threshold)
#eigenfns ( Phi.U_eigenfn) : The eigenfunctions of the smoothed pooled covariance
#eigenvals1 (lambda_U1_eigenval) : The estimated eigenvalues from projecting onto 
#the space of pooled eigenfunctions
#eigenvals2(lambda_U2_eigenval) :  
#sigma.noise2.1 : estimated error variance of data set 1
#sigma.noise2.2 : estimated error variance of data set 2

# require(mgcv)
# require(refund)
# require(MASS)

#'@import mgcv
#'@import refund
#'@import MASS

#'@export
########################################################################################################
### Smoothing of an existent covariance matrix
########################################################################################################
#Takes as input the two data sets and the treshold for chosing the number of non zero eigenvalues
#then uses refund pacakge fpca.sc to get estimates of the 
#smoothed covariances and then the common eigenfunctions of the pooled covariance. 
###

FPCA_covsmooth_fpcasc <- function(Y.1, Y.2, threshold = 0.99){
  
  #run fpca using the fpca.sc function in refund
  fpca.Y1 <- fpca.sc(Y.1,pve = .999999,var = TRUE)
  fpca.Y2 <- fpca.sc(Y.2,pve = .999999,var = TRUE)
  
  #get the number of curves in each data set for pooled covariance computation
  n.1 <- dim(Y.1)[1]
  n.2 <- dim(Y.2)[1]
  
  #get the number of observations in each curve
  m.1 <- dim(Y.1)[2]
  m.2 <- dim(Y.2)[2]
  
  #get  Gy1 and Gy2
  if(dim(fpca.Y1$efunctions)[2]>1){   
  sigma.est1 = (fpca.Y1$efunctions*sqrt(m.1)) %*% (diag(fpca.Y1$evalues)/m.1) %*% (t(fpca.Y1$efunctions)*sqrt(m.1))
  }
  if(dim(fpca.Y1$efunctions)[2] == 1){
    sigma.est1 = (fpca.Y1$efunctions*sqrt(m.1)) %*% (fpca.Y1$evalues/m.1) %*% (t(fpca.Y1$efunctions)*sqrt(m.1))
  }
  if(dim(fpca.Y2$efunctions)[2]>1){
  sigma.est2 = (fpca.Y2$efunctions*sqrt(m.2)) %*% (diag(fpca.Y2$evalues)/m.2) %*% (t(fpca.Y2$efunctions)*sqrt(m.2))
  }
  if(dim(fpca.Y2$efunctions)[2]==1){
    sigma.est2 = (fpca.Y2$efunctions*sqrt(m.2)) %*% (fpca.Y2$evalues/m.2) %*% (t(fpca.Y2$efunctions)*sqrt(m.2))
  }
  
  #get the error variances of the two data sets 
  sig2.Y1 <- fpca.Y1$sigma2
  sig2.Y2 <- fpca.Y2$sigma2

  #get the smoothed pooled covariance without error
  sig.estP <- ((n.1-1)*sigma.est1+(n.2-1)*sigma.est2)/(n.1+n.2-2)
  N  <- ncol(sig.estP)
 
  ####################################################################################
  #Get the eigenfunctions of the pooled covariance so we can use this as the 
  #common basis functions
  ####################################################################################
  #pick the ones for which the eigenvalues are positive so that we get a pd covariance
  #recall pd is equivalent to having positive eigenvalues 
  
  eig_sm_Gw = eigen(sig.estP) #this will have norm 1 so later it is multiplied by sqrt(N)
  fit.lambda = eig_sm_Gw$values
  #pick the ones for which the eigenvalues are positive so that we get a pd covariance
  #recall pd is equivalent to having positive eigenvalues 
  LU_est = length(which(fit.lambda >0)  )
  fit.lambda_pos = fit.lambda[1:LU_est]
  #est_Gw=eig_sm_Gw$vectors[, 1:LU_est] %*% diag(fit.lambda_pos) %*% t(eig_sm_Gw$vectors[, 1:LU_est])
  
  #K_U=est_Gw; 
  L_U = which((cumsum(fit.lambda_pos)/sum(fit.lambda_pos))>threshold)[1]
  L_U = max(L_U, 2)
  Phi.U_eigenfn  = eig_sm_Gw$vectors[, 1:L_U] *sqrt(N)

  
  #get the estimates of the eigenvalues, ie. covariance of the scores
  
  lambda_U1_eigenval = diag(t(eig_sm_Gw$vectors[,1:L_U])%*%sigma.est1%*%(eig_sm_Gw$vectors[,1:L_U]))
  lambda_U1_eigenval = lambda_U1_eigenval*(lambda_U1_eigenval>0)/N
  lambda_U2_eigenval = diag(t(eig_sm_Gw$vectors[,1:L_U])%*%sigma.est2%*%(eig_sm_Gw$vectors[,1:L_U]))
  lambda_U2_eigenval = lambda_U2_eigenval*(lambda_U2_eigenval>0)/N
  
  #get the seperate covarainces of both of the observed data sets
  cov.dat.1 <- sigma.est1 + (diag(rep(1,N))*sig2.Y1)
  cov.dat.2 <- sigma.est2 + (diag(rep(1,N))*sig2.Y2)
  #cov.smooth = est_Gw
  output = list(
    cov.dat.1       = cov.dat.1,
    cov.dat.2       = cov.dat.2, 
    L               = L_U, 
    eigenfns        = Phi.U_eigenfn, 
    eigenvals1      = lambda_U1_eigenval, 
    eigenvals2      = lambda_U2_eigenval, 
    sigma.noise2.1  = sig2.Y1, 
    sigma.noise2.2  =  sig2.Y2
  )
  output
}
#so the estimate of the covariance I want to use is 
#c <- FPCA_covsmooth$Phi.U_eigenfn%*%diag(FPCA_covsmooth$eigenvals)%*%t(FPCA_covsmooth$Phi.U_eigenfn)+(FPCA_covsmooth$sigma.noise%*%diag(rep(1,dim(FPCA_covsmooth$Phi.U_eigenfn)[1])))

#To debug: 
# x <- rnorm(2000,0,3)
# y <- rnorm(2000,0,4)
# d1 <- matrix(x,20,100)
# d2 <- matrix(y,20,100)
# cd1 <- cov(d1)
# cd2 <- cov(d2)
# n.1 <- dim(d1)[1]
# n.2 <- dim(d2)[1]
# test <- FPCA_covsmooth_fpcasc(cd1,cd2, threshold = 0.96)

matplot(cd1, type = 'l')
matplot(cd2, type = 'l')

fpca.Y1 <- fpca.sc(Y = cd1, pve = .999999, var = TRUE)
fpca.Y2 <- fpca.sc(Y.2, pve = .99, var = TRUE)

library(ggplot2)
library(reshape2)
data(cd4)
