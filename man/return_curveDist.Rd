\name{return_curveDist}
\alias{return_curveDist}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Returns similarity index between means of two samples of curves
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Returns the similarity index, introduced in fdakma, between means of two samples of curves
}
\usage{
return_curveDist(Mat1, Mat2, argvals, D1 = TRUE)
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
  \item{argvals}{
%%     ~~Describe \code{argvals} here~~
The abscissa, the x-axis or the time-axis. The length of this should be the same as the number of rows of Mat1 and Mat2
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
  \item{Dist_obs }{Distance (or similarity index) between means of the two curve samples}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
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
fdakma::kma.similarity, permute_curveDist
}
\examples{
data( growth, package = 'fda' )
Mat1 <- growth[['hgtm']]
Mat2 <- growth[['hgtf']]
Arguments <- growth[['age']]
SimIndex <- return_curveDist(
  Mat1 = Mat1, 
  Mat2 = Mat2, 
  argvals = Arguments
)
SimIndex
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ similarity }
\keyword{ distance }% __ONLY ONE__ keyword per line
