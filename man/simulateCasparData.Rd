\name{simulateCasparData}
\alias{simulateCasparData}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Data simulator
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Data simulator for CaSpaR simulations.
}
\usage{
simulateCasparData(n, p, cl, sz, cor, sig, peaks, valleys, binom = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{
Number of data points.
%%     ~~Describe \code{n} here~~
}
  \item{p}{
Number of predictors in the data set.
%%     ~~Describe \code{p} here~~
}
  \item{cl}{
Number of clusters of nonzero entries in the true vector of linear coefficients (beta)
%%     ~~Describe \code{cl} here~~
}
  \item{sz}{
Size of each cluster of nonzero entries in beta.
%%     ~~Describe \code{sz} here~~
}
  \item{cor}{
Correlation between inputs with each cluster of nonzero entries.
%%     ~~Describe \code{cor} here~~
}
  \item{sig}{
Noise added.
%%     ~~Describe \code{sig} here~~
}
  \item{peaks}{
Magnitude of one of the nonzero entries in each cluster within beta.
%%     ~~Describe \code{peaks} here~~
}
  \item{valleys}{
Magnitude of the remaining entries within each clster in beta.
%%     ~~Describe \code{valleys} here~~
}
  \item{binom}{
Depreciated.
%%     ~~Describe \code{binom} here~~
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
n = 50
G = simulateCasparData(n,25,2,5,0,1,6,3)
X = G$X
Y=G$Y
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
