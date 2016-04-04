\name{FPCA_covsmooth_fpcasc}
\alias{FPCA_covsmooth_fpcasc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Rreturns the score estimates of covariance function for 2-sample testing
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Inputs two data sets and returns the score estimates of covariance function that can be used for 2-sample testing with the 'testAD_fun.R' function. 
}
\usage{
FPCA_covsmooth_fpcasc(Y.1, Y.2, threshold = 0.99)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Y.1}{
%%     ~~Describe \code{Y.1} here~~
Y1 needs to be read in as n x d where n is the number of curves and d is the number of observations along the curve; Each row is a separate curve
}
  \item{Y.2}{
%%     ~~Describe \code{Y.2} here~~
Y2 needs to be read in as n x d where n is the number of curves and d is the number of observations along the curve; Each row is a separate curve
}
  \item{threshold}{
%%     ~~Describe \code{threshold} here~~
The proporton of variance explained by the eigenvectors that are used to explain the data and create the smoothed covaraince matrix
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
This code smooths the raw covariances first and then pools them. It uses the smoothing methodology 
from the refund package. 
}
\value{
%%  ~Describe the value returned
The output is a list of following elements
%%  If it is a LIST, use
  \item{cov.dat.1 }{smoothed covariance of data set 1}
  \item{cov.dat.2 }{smoothed covariance of data set 2}
  \item{L }{L_U : The number of eigenfunctions that will be used to explain the data (this depends on the threshold)}
  \item{eigenfns }{(Phi.U_eigenfn): The eigenfunctions of the smoothed pooled covariance}
  \item{eigenvals1 }{(lambda_U1_eigenval) : The estimated eigenvalues from projecting onto the space of pooled eigenfunctions}
  \item{eigenvals2 }{(lambda_U2_eigenval) : The estimated eigenvalues from projecting onto the space of pooled eigenfunctions}
  \item{sigma.noise2.1}{ estimated error variance of data set 1}
  \item{sigma.noise2.2}{ estimated error variance of data set 2}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
Pomann, G.M., Staicu, A.M., and Ghosh,S. (2016). A two-sample distribution-free test for functional data with application to a diffusion tensor imaging study of multiple sclerosis. Journal of the Royal Statistical Society: Series C (Applied Statistics).
}
\author{
%%  ~~who you are~~
Author of the function: Gina-Maria Pomann 

Author of this package: Subhrangshu Nandi (snandi@wisc.edu or nands31@gmail.com)
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
x <- rnorm(2000,0,3)
y <- rnorm(2000,0,10)
d1 <- matrix(x,20,100)
d2 <- matrix(y,20,100)
cd1 <- cov(d1)
cd2 <- cov(d2)
n.1 <- dim(d1)[1]
n.2 <- dim(d2)[1]
test <- FPCA_covsmooth_fpcasc(cd1,cd2, threshold = 0.96)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ fpca }
\keyword{ testing }
\keyword{ functional data }% __ONLY ONE__ keyword per line
