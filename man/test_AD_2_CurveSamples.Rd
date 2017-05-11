\name{test_AD_2_CurveSamples}
\alias{test_AD_2_CurveSamples}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Non parametric goodness-of-fit test of 2 samples of curves
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Returns a vector of p-values and test statistics of Anderson-Darling goodness-of-fit test of 2 samples of curves. It estimates fpca scores and conducts 2-sample AD test.
}
\usage{
test_AD_2_CurveSamples(Mat1, Mat2, varpct = 0.95)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Mat1}{
%%     ~~Describe \code{Mat1} here~~
Matrix 1 with each column being the discretized version of a separate curve 
}
  \item{Mat2}{
%%     ~~Describe \code{Mat1} here~~
Matrix 1 with each column being the discretized version of a separate curve 
}
  \item{varpct}{
%%     ~~Describe \code{varpct} here~~
Desired percentage of variance explained by the functional PCA scores
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
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
test_AD_2sample
}
\examples{
data( growth, package = 'fda' )
Mat1 <- growth[['hgtm']]
Mat2 <- growth[['hgtf']]
testOutput <- test_AD_2_CurveSamples(Mat1 = Mat1, Mat2 = Mat2, varpct = 0.95)
testOutput
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
