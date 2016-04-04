\name{return_supTstat}
\alias{return_supTstat}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Returns sup of T-type statistic between two matrices
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Returns sup of T-type statistic between two matrices, as described in fda::tperm.fd
}
\usage{
return_supTstat(Mat1, Mat2)
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
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
Returns a list of three elements
%%  If it is a LIST, use
  \item{TMean }{Mean of T-type statistics evaluated at each point (row) of the matrix}
  \item{TMedian }{Median of T-type statistics }
  \item{TSup }{Supremum of T-type statistics }
%% ...
}
\references{
%% ~put references to the literature/web site here ~
Graves, S., Hooker, G., & Ramsay, J. (2009). Functional data analysis with R and MATLAB. (Sec 10.5.1)
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
permute_pointwiseT, fda::tperm.fd
}
\examples{
require(fda)
data(growth)
Mat1 <- growth[['hgtm']]
Mat2 <- growth[['hgtf']]
TStats <- return_supTstat(Mat1, Mat2)
names(TStats)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Functional t-test }
\keyword{ testing }% __ONLY ONE__ keyword per line
