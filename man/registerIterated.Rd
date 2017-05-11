\name{registerIterated}
\alias{registerIterated}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Registers a sample of curves to its weighted average, in an iterated fashion
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Uses the \code{fda::register.fd} function to register a sample of curves, after smoothing them, iteratively. Checks for outliers before estimating the weighted average. The weights are estimated based on distance from a L1-median \code{robustX::L1median}
}
\usage{
registerIterated(
  dataToRegister, 
  Lambdas_ConstrainedWarping = c(0.001, 0.0001, 0.00005), 
  LambdaDefault = 0.00001,
  abscissaFrom, 
  abscissaTo, 
  abscissaIncrement, 
  basisOrder = 5, 
  basisBreakFreq = 3, 
  Lambdas_Roughness = exp(-5:5), 
  outlierTrimPct = 0.2, 
  N_Bootstrap_for_Outlier = 500,
  RE_REGISTER = FALSE,
  Convergence_Threshold = 1e-05,
  SimMeanDiff_Threshold = 0.0001, 
  MinSimilarityThreshold = 0.15, 
  MAX_ITERATION = 6
)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataToRegister}{
%%     ~~Describe \code{dataToRegister} here~~
A matrix of curves to register. Each column is a new curve. 
}
  \item{Lambdas_ConstrainedWarping}{
%%     ~~Describe \code{Lambdas_ConstrainedWarping} here~~
The roughness penalty for estimating the warping function. Smaller values will allow more undulations in the warping function. 
}
  \item{LambdaDefault}{
The default roughness penalty, set to 0.00001. The \eqn{\lambda} is set to this value if the function iterates longer than the list of \eqn{\lambda} provided in \code{Lambdas_ConstrainedWarping}
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
  \item{basisOrder}{
%%     ~~Describe \code{basisOrder} here~~
Order of B-spline basis functions
}
  \item{basisBreakFreq}{
%%     ~~Describe \code{basisBreakFreq} here~~
Frequency of knots of basis functions
}
  \item{Lambdas_Roughness}{
%%     ~~Describe \code{Lambdas_Roughness} here~~
Roughness penalty for smoothing using B-splines
}
  \item{outlierTrimPct}{
%%     ~~Describe \code{outlierTrimPct} here~~
Percentage to trim when detecting outliers
}
  \item{N_Bootstrap_for_Outlier}{
Number of bootstrap samples to draw, to estimate outlier curves  
}
  \item{RE_REGISTER}{
%%     ~~Describe \code{RE_REGISTER} here~~
Re-register the curves from previous iterations, or register the original noisy curves, at every iteration, to the updated consensus/template
}
  \item{Convergence_Threshold}{
The value for argument \code{conv} in function \code{register.fd}  
}
  \item{SimMeanDiff_Threshold}{
%%     ~~Describe \code{SimMeanDiff_Threshold} here~~
Defaults to 0.001, a criterion to stop the iterations
}
  \item{MinSimilarityThreshold}{
%%     ~~Describe \code{MinSimilarityThreshold} here~~
Minimum similarity between registered curve and template to estimate the final consensus
}
  \item{MAX_ITERATION}{
%%     ~~Describe \code{MAX_ITERATION} here~~
Maximum number of iterations
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
%% ...
  \item{dataToRegister}{Original dataset to register}
  \item{Regfd_Final}{Final registered data object. Output of \code{register.fd}}
  \item{registeredCurves}{A matrix of registered curves extracted from the registered object}
  \item{registeredCurves.D1}{A matrix of first derivatives of registered curves extracted from the registered object}
  \item{FinalConsensus}{Final consensus of registered curves}
  \item{curvesForConsensus}{Names of curves used to estimate the consensus. This excludes the detected outliers}
  \item{registeredCurvesAll}{A matrix of all registered curves, including outliers}
  \item{registeredCurvesAll.D1}{A matrix of first derivatives of all registered curves, including outliers}
  \item{Sim_toTemplate}{A matrix of similarities of curves to template}
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
Subhrangshu Nandi; snandi@wisc.edu or nands31@gmail.com
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
\code{register.fd}, \code{registerSingleIter}
}
\examples{
data( growth, package = 'fda' )
Mat1 <- growth[['hgtm']]
Arguments <- growth[['age']]

#Mat1Registered <- registerIterated(
#  dataToRegister = Mat1, 
#  Lambdas_ConstrainedWarping = c(0.001, 0.0001, 0.00005), 
#  LambdaDefault = 0.00001,
#  abscissaFrom = 1, 
#  abscissaTo = 16, 
#  abscissaIncrement = 0.5, 
#  basisOrder = 5, 
#  basisBreakFreq = 3, 
#  Lambdas_Roughness = exp(-5:5), 
#  outlierTrimPct = 0.05, 
#  N_Bootstrap_for_Outlier = 500,
#  RE_REGISTER = FALSE, 
#  Convergence_Threshold = 1e-05,
#  SimMeanDiff_Threshold = 0.001, 
#  MinSimilarityThreshold = 0.25, 
#  MAX_ITERATION = 4
#)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Registration }% use one of  RShowDoc("KEYWORDS")
\keyword{ Phase variability }% __ONLY ONE__ keyword per line
