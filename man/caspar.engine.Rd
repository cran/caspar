\name{caspar.engine}
\alias{caspar.engine}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Internal function.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Internal function.
}
\usage{
caspar.engine(Y.train, X.train, Y.test = NULL, X.test = NULL, kp = 0, mix = 1, d.fun, trace = TRUE, eval = "cv")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Y.train}{
%%     ~~Describe \code{Y.train} here~~
Training data.
}
  \item{X.train}{
%%     ~~Describe \code{X.train} here~~
Training data.
}
  \item{Y.test}{
%%     ~~Describe \code{Y.test} here~~
Test data.
}
  \item{X.test}{
%%     ~~Describe \code{X.test} here~~
Test data.
}
  \item{kp}{
%%     ~~Describe \code{kp} here~~
Kernel parameter.
}
  \item{mix}{
%%     ~~Describe \code{mix} here~~
Mixing parameter.
}
  \item{d.fun}{
%%     ~~Describe \code{d.fun} here~~
Distance function.
}
  \item{trace}{
%%     ~~Describe \code{trace} here~~
Output?
}
  \item{eval}{
%%     ~~Describe \code{eval} here~~
CV, BIC, AIC.
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
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
