\name{return_supDiff}
\alias{return_supDiff}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Returns sup of difference between row means of two matrices
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Returns sup of difference between row means of two matrices. 
}
\usage{
return_supDiff(Mat1, Mat2)
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
  \item{DMean }{Mean of D-type statistics evaluated at each point (row) of the matrix}
  \item{DMedian }{Median of D-type statistics }
  \item{DSup }{Supremum of D-type statistics }
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
permute_pointwiseD, return_supTstat
}
\examples{
data( growth, package = 'fda' )
Mat1 <- growth[['hgtm']]
Mat2 <- growth[['hgtf']]
Diff <- return_supDiff(Mat1, Mat2)
names(Diff)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ hypothesis-testing }% __ONLY ONE__ keyword per line
