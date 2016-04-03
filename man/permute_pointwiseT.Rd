\name{permute_pointwiseT}
\alias{permute_pointwiseT}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Permutation test between two sets of curves
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Permutation test between two sets of curves
}
\usage{
permute_pointwiseT(Mat1, Mat2, Nperm = 200, argvals, q = 0.05, 
                   returnPlot = TRUE, TitleText = "")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Mat1}{
%%     ~~Describe \code{Mat1} here~~
Matrix with each column being the discretized version of a separate curve 
}
  \item{Mat2}{
%%     ~~Describe \code{Mat2} here~~
Matrix with each column being the discretized version of a separate curve 
}
  \item{Nperm}{
%%     ~~Describe \code{Nperm} here~~
Number of permutations. Default 200. 
}
  \item{argvals}{
%%     ~~Describe \code{argvals} here~~
The abscissa, the x-axis or the time-axis. The length of this should be the same as the number of rows of Mat1 and Mat2
}
  \item{q}{
%%     ~~Describe \code{q} here~~
0.05 or an appropriate \eqn{\alpha} level 
}
  \item{returnPlot}{
%%     ~~Describe \code{returnPlot} here~~
Whether to return the permutation distribution plot
}
  \item{TitleText}{
%%     ~~Describe \code{TitleText} here~~
A text intended to be appended to ggtitle
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
This is the permutation test, testing the hypothesis \eqn{H_0: F_{Y_1} = F_{Y_2}}
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
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
This is an adaptation of the function fda::tperm.fd. The theory was first described in book mentioned in the references.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
fda::tperm.fd
}
\examples{
library(fda)
data(growth)
Mat1 <- growth[['hgtm']]
Mat2 <- growth[['hgtf']]
Arguments <- growth[['age']]
PermTestResults <- permute_pointwiseT(
  Mat1 = Mat1, 
  Mat2 = Mat2, 
  Nperm = 200,
  argvals = Arguments,
  q = 0.05,
  returnPlot = TRUE,
  TitleText = 'Growths of boys and girls (Berkeley growth data)'
)
names(PermTestResults)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ permutation }
\keyword{ testing }
\keyword{ curves }% __ONLY ONE__ keyword per line
