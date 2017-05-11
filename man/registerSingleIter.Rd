\name{registerSingleIter}
\alias{registerSingleIter}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Registers a sample of curves to its cross-sectional average
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Uses the \code{fda::register.fd} function to register a sample of curves, after smoothing them.
}
\usage{
registerSingleIter(
  dataToRegister, 
  Lambda_ConstrainedWarping, 
  abscissaFrom, 
  abscissaTo, 
  abscissaIncrement, 
  Convergence_Threshold = 1e-05,
  basisOrder = 4, 
  basisBreakFreq = 3, 
  Lambdas_Roughness = exp(-5:5)
)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataToRegister}{
%%     ~~Describe \code{dataToRegister} here~~
A matrix of curves to register. Each column is a new curve. 
}
  \item{Lambda_ConstrainedWarping}{
%%     ~~Describe \code{Lambda_ConstrainedWarping} here~~
The roughness penalty for estimating the warping function. Smaller values will allow more undulations in the warping function. 
}
  \item{abscissaFrom}{
%%     ~~Describe \code{abscissaFrom} here~~
Mininum value of abscissa
}
  \item{abscissaTo}{
%%     ~~Describe \code{abscissaTo} here~~
Maximum value of abscissa
}
  \item{abscissaIncrement}{
%%     ~~Describe \code{abscissaIncrement} here~~
Increment of abscissa
}
  \item{Convergence_Threshold}{
The value for argument \code{conv} in function \code{register.fd}  
}
  \item{basisOrder}{
%%     ~~Describe \code{basisOrder} here~~
Order of B-spline bases
}
  \item{basisBreakFreq}{
%%     ~~Describe \code{basisBreakFreq} here~~
Frequency of knots of basis functions
}
  \item{Lambdas_Roughness}{
%%     ~~Describe \code{Lambdas_Roughness} here~~
Roughness penalty for smoothing using B-splines
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
  \item{dataToRegister}{Original dataset to register}
  \item{Regfd_Final}{Final registered data object. Output of \code{register.fd}}
  \item{registeredCurves}{A matrix of registered curves extracted from the registered object}
  \item{registeredCurves.D1}{A matrix of first derivatives of registered curves extracted from the registered object}
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
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
Arguments <- growth[['age']]

#Mat1Registered <- registerSingleIter(
#  dataToRegister = Mat1, 
#  Lambda_ConstrainedWarping = 0.005, 
#  abscissaFrom = 1, 
#  abscissaTo = 16, 
#  abscissaIncrement = 0.5, 
#  Convergence_Threshold = 1e-05,
#  basisOrder = 5, 
#  basisBreakFreq = 3, 
#  Lambdas_Roughness = exp(-5:5), 
#)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Registration }% use one of  RShowDoc("KEYWORDS")
\keyword{ Phase variability }% __ONLY ONE__ keyword per line
