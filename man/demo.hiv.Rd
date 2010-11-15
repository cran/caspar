\name{demo.hiv}
\alias{demo.hiv}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
CaSpaR demo: HIV data.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Performs a demonstration on the HIV data presented in Structured, Sparse Regression With Application to HIV Drug Resistance.  AOAS, to appear.  See: http://arxiv.org/abs/1002.3128. 
}
\usage{
demo.hiv(res = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{res}{
Integer in the range 1 to 7 inclusive.  There are 7 HIV data sets, each corresponding to a different HIV PI (protease inhibitor) drug.
%%     ~~Describe \code{res} here~~
}
}
\details{
There are 7 HIV data sets.  These can be accessed directly using the data() command, with arguments corresponding to data sets 1-7 as follows: 1 == "HIVDAT1APV", 2 == "HIVDAT2ATV", 3 == "HIVDAT3IDV", 4 = "HIVDAT4LPV", 5 == "HIVDAT5NFV", 6 == "HIVDAT6RTV", 7 == "HIVDAT7SQV".  Once called, an object named hiv.dat will be in the workspace.  This object has the elements $X, $Y, representing the inputs and response.  The element $iX is the inputs coded as indicator variables.  The elements $geo and $grp define the sequence geometry for $iX.   $grp is a vector of integers giving the number of levels of each factor in $X, minus one.  $geo is a binary vector indicating which positions of the 99 length are excluded since these positions have no variation within the sample.
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Returns a CaSpaR fit of the appropriate HIV data set.
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
demo.hiv(res=2)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
