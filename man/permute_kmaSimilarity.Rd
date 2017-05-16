\name{permute_kmaSimilarity}
\alias{permute_kmaSimilarity}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Permutation test between two sets of curves
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Permutation test between two sets of curves using similarity index
}
\usage{
permute_kmaSimilarity(Mat1, Mat2, Nperm = 200, argvals, q = 0.05, D1 = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Mat1}{
%%     ~~Describe \code{Mat1} here~~
Matrix 1 with each column being the discretized version of a separate curve 
}
  \item{Mat2}{
%%     ~~Describe \code{Mat2} here~~
Matrix 2 with each column being the discretized version of a separate curve 
}
  \item{Nperm}{
%%     ~~Describe \code{Nperm} here~~
Number of permutations. Default 200
}
  \item{argvals}{
%%     ~~Describe \code{argvals} here~~
The abscissa, the x-axis or the time-axis. The length of this should be the same as the number of rows of Mat1 and Mat2
}
  \item{q}{
%%     ~~Describe \code{q} here~~
0.05 or an appropriate \eqn{\alpha} level 
}
  \item{D1}{
%%     ~~Describe \code{D1} here~~
\code{TRUE/FALSE} whether the matrix contains the curves of derivatives of the curves. \code{TRUE} if contains derivatives.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
Parodi, A., et al. "fdakma: Functional data analysis: K-mean alignment." R package version 1.1 (2014).

https://cran.r-project.org/web/packages/fdakma/index.html
}
\author{
%%  ~~who you are~~
Subhrangshu Nandi; Statistics PhD student, UW Madison; snandi@wisc.edu or nands31@gmail.com
}
\note{
%%  ~~further notes~~
This is an extension of the function fdakma::kma.similarity. This is using \eqn{\rho} to test similarity between two sets of curves.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
fdakma::kma.similarity
}
\examples{
data( growth, package = 'fda' )
Mat1 <- growth[['hgtm']]
Mat2 <- growth[['hgtf']]
Arguments <- growth[['age']]
PermTestResults <- permute_kmaSimilarity(
  Mat1 = Mat1, 
  Mat2 = Mat2, 
  Nperm = 200,
  argvals = Arguments,
  q = 0.05
)
names(PermTestResults)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ permutation }
\keyword{ hypothesis-testing }
\keyword{ similarity }% __ONLY ONE__ keyword per line
