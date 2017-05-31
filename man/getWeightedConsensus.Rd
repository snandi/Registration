\name{getWeightedConsensus}
\alias{getWeightedConsensus}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Estimate a weighted mean/median of a sample of Curves
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Instead of a simple rowMeans, this functions estimates a Weighted consensus of a sample of Curves
}
\usage{
getWeightedConsensus(Curves, abscissa)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Curves}{
%%     ~~Describe \code{Curves} here~~
A matrix of curves, with the columns being the curves
}
  \item{abscissa}{
%%     ~~Describe \code{abscissa} here~~
The x-axis
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
First an \code{L1-median} is estimated, then \code{kma.similarity} is estimated between the curves and the L1-median. This similarity is used as weights to estimate weighted consensus.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
\item{Mean_Weighted}{Weighted mean}
\item{Median_Weighted}{Weighted median}
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
\code{L1median}, \code{rowMeans}, \code{rowWeightedMedians}, \code{rowWeightedMeans}
}
\examples{
data( growth, package = 'fda' )
Mat1 <- growth[['hgtm']]
Arguments <- growth[['age']]
%getWeightedConsensus( Curves = Mat1, abscissa = Arguments )
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ consensus }% use one of  RShowDoc("KEYWORDS")
\keyword{ fda }% __ONLY ONE__ keyword per line
\keyword{ weighted mean }
\keyword{ weighted median }
\keyword{ L1median }


