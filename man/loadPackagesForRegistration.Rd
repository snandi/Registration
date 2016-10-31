\name{loadPackagesForRegistration}
\alias{loadPackagesForRegistration}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Loads packages required by Registration. If not present, then installs them from user specified CRAN mirror.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
loadPackagesForRegistration(CRANMirror = 83, tryToInstall = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{CRANMirror}{
%%     ~~Describe \code{CRANMirror} here~~
\code{CRANMirror} is the Cran mirror from where the package will be installed if
}
  \item{tryToInstall}{
%%     ~~Describe \code{tryToInstall} here~~
\code{tryToInstall} is by default \code{FALSE}. If \code{TRUE}, it will try to install a package 
that is absent
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
\item{PackageList }{Returns a list of two vectors, 'Packages' and 'Packages_Par'. In this version, both are the same packages that are required by Registration}
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
\code{library}, \code{require}
}
\examples{
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ install }% use one of  RShowDoc("KEYWORDS")
\keyword{ library }% __ONLY ONE__ keyword per line
