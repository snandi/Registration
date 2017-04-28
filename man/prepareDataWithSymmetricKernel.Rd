\name{prepareDataWithSymmetricKernel}
\alias{prepareDataWithSymmetricKernel}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Prepares a design matrix with a symmetric kernel
}
\description{
Prepares a design matrix adding a symmetric kernel (either Gaussian or Uniform) with columns as counts of genomic k-mers
}
\usage{
prepareDataWithSymmetricKernel(
  SeqCompData, 
  kernelType = c("Uniform", "Gaussian"), 
  colnamesNonSeq = c("consensusMean", "Interval", "subInterval", "subIntervalPixels")
)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{SeqCompData}{
A design matrix \eqn{(n x p)} with each rows as responses and columns as covariates. The columns contains counts of k-mers
}
  \item{kernelType}{
Gaussian or Uniform
}
  \item{colnamesNonSeq}{
Names of columns that do not correspond to the Seq Comp covariates
}
}
\details{
The kernel is 5 units wide. So, 2 rows before and after a response row are added as additional columns. This was done because unlike regression methods random forest or gradient boosting methods do not have the ability to incorporate kernels. 
}
\value{
Returns a matrix \eqn{(n-4 x 3p)}
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
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
}
\examples{
data( SeqCompData )
segmentDataWithKernel <- prepareDataWithSymmetricKernel( 
    SeqCompData =  SeqCompData, 
    kernelType = KERNEL_TYPE, 
    colnamesNonSeq = c( 'Chr', 'Segment', 'segmentPixels' )
)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kernel }% use one of  RShowDoc("KEYWORDS")

