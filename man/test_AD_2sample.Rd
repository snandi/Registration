\name{test_AD_2sample}
\alias{test_AD_2sample}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
2 sample Anderson Darling test
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Returns a vector of p-values and test statistics of a 2 sample Anderson-Darling goodness-of-fit test
}
\usage{
test_AD_2sample(datX, datY, p)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{datX}{
%%     ~~Describe \code{datX} here~~
Matrix of p columns to be tested with datY
}
  \item{datY}{
%%     ~~Describe \code{datY} here~~
Matrix of p columns
}
  \item{p}{
%%     ~~Describe \code{p} here~~
The number of columns of the matrices
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
This is a simple wrapper around the original adk.test function in the package adk, which is currently deprecated.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
  \item{pval }{Vector of p-values, of length p}
  \item{testStat }{Vector of test statistics, of length p}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
Pomann, G.M., Staicu, A.M., and Ghosh,S. (2016). A two-sample distribution-free test for functional data with application to a diffusion tensor imaging study of multiple sclerosis. Journal of the Royal Statistical Society: Series C (Applied Statistics).

Scholz, F. W., & Stephens, M. A. (1987). K-sample Anderson–Darling tests. Journal of the American Statistical Association, 82(399), 918-924.
}
\author{
%%  ~~who you are~~
Subhrangshu Nandi; Statistics PhD student, UW Madison; snandi@wisc.edu or nands31@gmail.com
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
adk.test
}
\examples{
  U1 <- stats::rnorm(n = 20, mean = 2, sd = 1)
  U2 <- stats::rnorm(n = 20, mean = 1, sd = 2)
  
  V1 <- stats::rnorm(n = 20, mean = 0, sd = 1)
  V2 <- stats::rnorm(n = 20, mean = 1, sd = 2)
  
  Mat1 <- cbind(U1, U2)
  Mat2 <- cbind(V1, V2)
  testOutput <- test_AD_2sample( datX = Mat1, datY = Mat2, p = 2 )
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Anderson-Darling }
\keyword{ testing }% __ONLY ONE__ keyword per line
