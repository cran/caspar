\name{caspar.main}
\alias{caspar.main}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Clustered and Sparse Regression
}
\description{
All-in-one caspar function.  Finds a sparse lienar model where the included covariates are clustered with respect to a distance function defined between covariates.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
caspar.main(X, Y, distanceFunction, kernelVector, mixVector, cvFolds, trace = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
%%     ~~Describe \code{X} here~~
Data matrix X
}
  \item{Y}{
%%     ~~Describe \code{Y} here~~
Response vector Y
}
  \item{distanceFunction}{
%%     ~~Describe \code{distanceFunction} here~~
There are two valid inputs for this argument.  (1) A function of the form function(ONE,TWO), where ONE and TWO are integers, and function(ONE,TWO) returns a distance between the pair of covariates indexed by ONE and TWO.  (2) A pairwise matrix between covariates.  That is, entry i,j of the matrix contains the distance between predictor i and j.
}
  \item{kernelVector}{
%%     ~~Describe \code{kernelVector} here~~
Vector of kernel bandwidth parameters to try.  CaSpaR uses the Epanechnikov kernel, defined as $f(u) = (3/4)(1-u^2) I_{|u| < 1}$.  For bandwidth parameter $h$, this kernel is redefined as $f_h(u) = f(u/h)$.  This kernel defines the clustering behavior of the resulting model.  CaSpaR will try each of these parameters, and return the model using the parameter with the best performance.
}
  \item{mixVector}{
%%     ~~Describe \code{mixVector} here~~
Vector of mixing parameters to try.  The full kernel is a mix between the Epanechnikov kernel and a uniform distribution.  That is, for mixing parameter $a$, we have $k(x) = a + (1-a)f(x)$, where $f(x)$ is the epanechnikov kernel.  CaSpaR will try each of these parameters, and return the model using the parameter with the best performance.
}
  \item{cvFolds}{
%%     ~~Describe \code{cvFolds} here~~
Integer, number of folds for cross validation.
}
  \item{trace}{
%%     ~~Describe \code{trace} here~~
Logical, should the program print its progress?
}
}
\details{
caspar.main performs a grid search, looking over all combinations of the inputs mixVector and kernelVector.  The program returns the model with the combination of parameters corresponding to the best performance.
%%  ~~ If necessary, more details than the description above ~~
}
\value{
\item{finalModel}{An lm() result of the final selected model.}
\item{grid}{The results of the grid search.  Contains two objects: cv.grid and cv.step.  Both are matrices of the results of the grid search, grid contains the cv scores, and step indicates the step of stepwise regression on which the procedure stopped.  For both matrices, the rows correspond to differing values of the bandwidth parameter and the columns to differing values of the mixing parameter.}
\item{call}{Call of this function.}
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
demo.caspar()

n = 45
{G = simulateCasparData(n,25,2,5,0,1,6,3); X = G$X; Y=G$Y}
d.fun.demo = function(ONE,TWO) abs(ONE-TWO)
kp.vec     = 2:4
mix.vec    = (1:10)/10
mod        = caspar.main(X,Y,d.fun.demo,kp.vec,mix.vec,cvFolds=5,trace=TRUE)

## Code for demo.caspar:
##  n = 45
##  {G = simulateCasparData(n,25,2,5,0,1,6,3); X = G$X; Y=G$Y}
##  nz = which(G$truth != 0)
##  d.fun.demo = function(ONE,TWO) abs(ONE-TWO)
##  ## d.fun.demo = abs(outer(1:ncol(X),1:ncol(X),"-"))
##  kp.vec     = 2:4
##  mix.vec    = (1:10)/10
##  out  = caspar.grid(X,Y,5,d.fun.demo,kp.vec = kp.vec, mix.vec = mix.vec, TRUE)
##  best = which(out$cv.grid == min(out$cv.grid),arr.ind = TRUE)[1,]

##  take = out$cv.step[best[1],best[2]]
##  out2 = caspar.engine(Y,X,NULL,NULL,kp = kp.vec[best[1]], mix = mix.vec[best[2]],d.fun.demo,trace=FALSE)
  
##  Covariate = X[,sort(out2$sequ[1:take])]
##  FL = function(vec)
##    vec[c(length(vec),1:(length(vec)-1))]
##  strLF <- function(x)
##    sapply(lapply(strsplit(x, NULL), FL), paste, collapse="")
##  colnames(Covariate) = strLF(sapply(sort(out2$sequ[1:take]),paste,"X ",sep=""))
##  print(summary(lm(Y~Covariate)))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
