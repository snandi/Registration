\name{getFunctionalOutliers}
\alias{getFunctionalOutliers}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Gets functional outliers from a sample of curves
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Gets functional outliers from a sample of curves using methods described in fda.usc
}
\usage{
getFunctionalOutliers( 
  Curves, 
  Xaxis, 
  Names = list(main = 'Main', xlab = 'PixelPosition', ylab = 'Intensity'),
  DepthType = c( 'FM', 'Mode', 'RTukey', 'RProj' ),
  N_Bootstrap = 500,
  Trim = c( 'Yes', 'No' ),
  TrimPct = 0.05
 )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Curves}{
%%     ~~Describe \code{Curves} here~~
A matrix (or dataframe) of curves, with each column being a separate curve
}
  \item{Xaxis}{
%%     ~~Describe \code{Xaxis} here~~
The abscissa
}
  \item{Names}{
%%     ~~Describe \code{Names} here~~
list(main = 'Main', xlab = 'PixelPosition', ylab = 'Intensity')
}
  \item{DepthType}{
%%     ~~Describe \code{DepthType} here~~
Either of the four \code{c( 'FM', 'Mode', 'RTukey', 'RProj' )}
}
  \item{N_Bootstrap}{
%%     ~~Describe \code{N_Bootstrap} here~~
Number of boostrap samples. Defaults to 500.
}
  \item{Trim}{
%%     ~~Describe \code{Trim} here~~
Whether to Trim the samples or not. 
}
  \item{TrimPct}{
%%     ~~Describe \code{TrimPct} here~~
How much to Trim the samples or not
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
Returns a vector of column names that are detected as outliers.
}
\references{
%% ~put references to the literature/web site here ~
Febrero-Bande, M., & Oviedo de la Fuente, M. (2012). Statistical computing in functional data analysis: the R package fda. usc. Journal of Statistical Software, 51(4), 1-28.

Cuevas A, Febrero M, Fraiman R. 2006. On the use of bootstrap for estimating functions with functional
data. Computational Statistics and Data Analysis 51: 1063-1074.

Febrero-Bande, M., Galeano, P., and Gonzalez-Manteiga, W. (2008). Outlier detection in functional
data by depth measures with application to identify abnormal NOx levels. Environmetrics 19, 4,
331-345.

Febrero-Bande, M., Galeano, P. and Gonzalez-Manteiga, W. (2007). A functional analysis of NOx
levels: location and scale estimation and outlier detection. Computational Statistics 22, 3, 411-427.

Febrero-Bande, M., Oviedo de la Fuente, M. (2012). Statistical Computing in Functional Data
Analysis: The R Package fda.usc. Journal of Statistical Software, 51(4), 1-28. 

http://www.jstatsoft.org/v51/i04/
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
\code{outliers.depth.trim}, \code{outliers.depth.pond}
}
\examples{
require(fda)
data(growth)
Mat1 <- growth[['hgtm']]
Arguments <- growth[['age']]

getFunctionalOutliers( 
  Curves = Mat1, 
  Xaxis = Arguments, 
  Names = list(main = 'Main', xlab = 'PixelPosition', ylab = 'Intensity'),
  DepthType = 'FM',
  N_Bootstrap = 500,
  Trim = 'Yes',
  TrimPct = 0.05
 )
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ functional data }% use one of  RShowDoc("KEYWORDS")
\keyword{ outliers }% __ONLY ONE__ keyword per line
\keyword{ depth }