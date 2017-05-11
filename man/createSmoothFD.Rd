\name{createSmoothFD}
\alias{createSmoothFD}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Create a smooth functional data object using bsplines
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
For a curve of a matrix this function creates a smooth functional data object as per the definition in R package fda, using bspline basis. 
}
\usage{
createSmoothFD(
  curvesToSmooth, 
  abscissa, 
  lambdas = exp(-5:5), 
  basisBreakFreq = 3, 
  basisOrder = 4, 
  pbasis = NULL
)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{curvesToSmooth}{
%%     ~~Describe \code{curvesToSmooth} here~~
A matrix or a vector to smooth
}
  \item{abscissa}{
%%     ~~Describe \code{abscissa} here~~
The x-axis. Should of of the same length as the vector curvesToSmooth or the number of rows in matrix curvesToSmooth
}
  \item{lambdas}{
%%     ~~Describe \code{lambdas} here~~
The different roughness penalty to try during smoothing
}
  \item{basisBreakFreq}{
%%     ~~Describe \code{basisBreakFreq} here~~
Defaults to 3. More frequent could overfit, less frequent could over-smooth.
}
  \item{basisOrder}{
%%     ~~Describe \code{basisOrder} here~~
Defaults to 4. If more noisy, then this can be higher.
}
  \item{pbasis}{
%%     ~~Describe \code{pbasis} here~~
User-defined pbasis if user-defined basisObject is desired. Else defaults to NULL.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
\item{curvesSmooth}{FDA object of smoothed curves}
\item{lambdaBest}{lambda associated with minimum GCV}
}
\references{
%% ~put references to the literature/web site here ~
Graves, S., Hooker, G., & Ramsay, J. (2009). Functional data analysis with R and MATLAB.
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
}
\examples{
data( growth, package = 'fda' )
Mat1 <- growth[['hgtm']]
Mat1Smooth <- createSmoothFD( 
  curvesToSmooth = Mat1, 
  abscissa       = 1:nrow( Mat1 ), 
  lambdas        = exp( -5:5 ), 
  basisBreakFreq = 3,
  basisOrder     = 4,
  pbasis         = NULL
)

plot( (Mat1Smooth$curvesSmooth)$fd )
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ fda }% use one of  RShowDoc("KEYWORDS")
\keyword{ create.basis }% __ONLY ONE__ keyword per line
\keyword{ create.bspline.basis }
