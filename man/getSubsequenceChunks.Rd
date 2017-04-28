\name{getSubsequenceChunks}
\alias{getSubsequenceChunks}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Get subsequence chunks from a long sequence
}
\description{
This function takes in a DNA sequence "agggctattaaccctttaaa" in the form of an array of characters and returns smaller subsequence chunks based on user specified number
}
\usage{
getSubsequenceChunks( Sequence, numChunks, force.number.of.groups = TRUE )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Sequence}{
A DNA sequence in the form of an array of characters: "a" "g" "g" "c" "t"
}
  \item{numChunks}{
Number of subsequence chunks desired
}
  \item{force.number.of.groups}{
If \code{TRUE}, it will append the remainder of the elements to the last chunk
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Returns a list of subsequence chunks.
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
%library( seqinr )
getSubsequenceChunks( 
  Sequence = seqinr::s2c( "agggctattaaccctttaaa" ), 
  numChunks = 3, 
  force.number.of.groups = TRUE )
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~DNA sequence }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~subsequence }% __ONLY ONE__ keyword per line
