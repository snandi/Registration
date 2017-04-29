\name{getKmaSimilarityWithTemplate}
\alias{getKmaSimilarityWithTemplate}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
To get kma.similarity for a matrix of curves with a template curve 
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
To use kma.similarity from fdakma package, for a matrix of curves with a template curve and return a vector of \code{similarity}
}
\usage{
getKmaSimilarityWithTemplate(
  Mat, 
  Template, 
  Xaxis, 
  similarity.method = c("d1.pearson", "d1.L2", "d1.L2.centered"), 
  Deriv = FALSE
)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Mat}{
%%     ~~Describe \code{Mat} here~~
A matrix of curves, with the columns being the curves
}
  \item{Template}{
%%     ~~Describe \code{Template} here~~
A template curve with the same length as the curves
}
  \item{Xaxis}{
%%     ~~Describe \code{Xaxis} here~~
The abscissa
}
  \item{similarity.method}{
%%     ~~Describe \code{similarity.method} here~~
Could be one of the three \code{c( "d1.pearson", "d1.L2", "d1.L2.centered" )}
}
  \item{Deriv}{
%%     ~~Describe \code{Deriv} here~~
Whether the matrix contains Derivatives or not. If not, then derivatives are estimated before proceeding.
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
\item{Sim}{ A vector of \code{similarity} of each curve with the template curve}
}
\references{
%% ~put references to the literature/web site here ~
Parodi, A., et al. "fdakma: Functional data analysis: K-mean alignment." R package version 1.1 (2014).

https://cran.r-project.org/web/packages/fdakma/index.html
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
\code{fdakma::kma.similarity}
}
\examples{
data( growth, package = 'fda')
Mat1 <- growth[['hgtm']]
Arguments <- growth[['age']]
Template <- rowMeans( Mat1 )
getKmaSimilarityWithTemplate(
  Mat = Mat1, 
  Template = Template, 
  Xaxis = Arguments, 
  similarity.method = "d1.pearson", 
  Deriv = FALSE
)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ similarity }% use one of  RShowDoc("KEYWORDS")
\keyword{ functional data }% __ONLY ONE__ keyword per line
\keyword{ fda }
