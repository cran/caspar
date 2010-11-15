## CaSpaR: Cluster sparse regression
## Distance function version
##
## Dancsi Percival
## Performs stepwise and CaSpaR methods w/
## stetson kernel -- epanechnikov kernel only right now...

####################################
## caspar.engine function INPUTS:
## Y.train: univariate response Y
## X.train: multivariate response X
## Y.test:  univariate response Y -- test data, if null no test set
## X.test:  multivariate response X -- test data, if null no test set
##          if both Y.test and X.test the cv score is not computed
##
## kp:      kernel parameter (bandwidth, h in paper)
## mix:     mixing parameter (\in (0,1), alpha in paper) 1 => stepwise regression
## d.fun:   either:
##          (a) function of the form function(ONE,TWO),
##               maps from 1:ncol(X) to \mathbb{R}
##          (b) pairwise distance matrix, can also be a dist() object
##         defines the geometry of your problem
## trace:  do you want a lot of things printed?
#######################################

caspar.engine = function(Y.train, X.train, Y.test = NULL, X.test = NULL, kp = 0, mix = 1,  d.fun, trace = TRUE,eval="cv"){
  eval = c("cv","AIC","BIC")[match(eval,c("cv","AIC","BIC"),1)]

  ## Training data
  X = X.train
  Y = Y.train
  n = nrow(X)
  p = ncol(X)
    
  ## Set up predictor distances
  ## locs: the pairwise distance between predictors
  if(class(d.fun) != "function")
    locs = as.matrix(d.fun)
  if(class(d.fun) == "function"){
    locs = matrix(0,nr=p,nc=p)
    for(i in 1:p)
      for(j in 1:p)
        locs[i,j] = d.fun(i,j)
  }  
  
  ## set up test set...
  XT = X.test
  YT = Y.test
  ## if no data provided, set evaluation to BIC
  ## and set up a phony test set
  if(is.null(Y.test) || is.null(X.test)){
    XT   = as.matrix(scale(X))
    YT   = as.matrix(scale(Y))
    X    = as.matrix(scale(X))
    Y    = as.matrix(scale(Y))
    if(eval == "cv") eval = "BIC"
  }

  model = casparCV(Y,X,YT,XT,kp,mix,locs = locs, verbose = trace)

  ## the first step is always the intercept...
  scores   = model[[12]][-1]
  add.seq  = model[[11]][-1]
  
  if(eval == "AIC")
    scores  = 2*(1:length(scores)) + n*log(n*scores)
  if(eval == "BIC")
    scores  = log(n)*(1:length(scores)) + n*log(n*scores)
  
  if(eval == "cv")  return(list(cv = scores, sequ = add.seq))
  if(eval != "cv") return(list(score = scores, sequ = add.seq))
}

################
## function caspar.grid
## the main function for CaSpaR
## Arguments:
##
## X:       Data matrix (n rows, p columns)
## Y:       Univariate continuous response (vector of length n)
## folds:   Number of folds for cross validation
## d.fun:   distance function (format: d.fun(ONE,TWO) maps 1:p to \mathbb{R}
## kp.vec:  vector of kernel bandwidths to try in grid search
## mix.vec: vector of mising parameters to try in grid search
## trace:   print the progress?
##
## returns:
## cv.grid: the cv scores for the grid of kp.vec (rows), mix.vec (cols)
## cv.step: the step with the best cv scores for the same grid
##
##

caspar.grid = function(X,Y,folds,d.fun,kp.vec,mix.vec,trace = TRUE){
  ## cv.grid: for the cv scores
  ## cv.step: for storing the step with the best cv score
  cv.grid           = matrix(0,nr=length(kp.vec),nc=length(mix.vec))
  rownames(cv.grid) = kp.vec
  colnames(cv.grid) = mix.vec
  cv.step           = cv.grid
  
  ## Parameters
  n        = nrow(X)
  p        = ncol(X)
  nam      = names(data.frame(X))
  X        = data.frame(scale(X))
  Y        = scale(Y)
  names(X) = nam
  X        = as.matrix(X)

  ## Set up predictor distances
  ## locs: the pairwise distance between predictors
  if(class(d.fun) != "function")
    locs = as.matrix(d.fun)
  if(class(d.fun) == "function"){
    locs = matrix(0,nr=p,nc=p)
    for(i in 1:p)
      for(j in i:p)
        locs[i,j] = d.fun(i,j)
  }

  ## valid: defines how many cv scores will be OK to take
  ## tempmat: stores all the cv scores
  valid   = min(p,length(Y) - (ceiling(length(Y)/folds)) - 1)
  tempmat = array(0,dim=c(folds,length(kp.vec),length(mix.vec),valid))

  groups = split(sample(1:n),rep(1:folds,length=n))

  for(ii in 1:folds){
    omit    = groups[[ii]]
    data    = X[-omit,]
    if(n-length(omit) == 1)
      data = as.matrix(data,nr = 1)
    resp    = Y[-omit]
    for(i in 1:length(kp.vec)){
      for(j in 1:length(mix.vec)){
        if(trace) cat("\nGRID: h=",kp.vec[i]," a=",mix.vec[j]," fold=",ii)

        ## compute the cv scores
        test.X = X[omit,]
        if(length(omit) == 1)
          test.X = t(as.matrix((test.X),nr=1))

        tempmat[ii,i,j,]  = caspar.engine(resp,data,Y[omit],test.X,kp.vec[i],mix.vec[j],d.fun,trace,eval="cv")$cv[1:valid]

        cv.grid[i,j] = cv.grid[i,j] + min(tempmat[ii,i,j,],na.rm=TRUE)/folds
      }
    }
  }
  for(i in 1:length(kp.vec))
    for(j in 1:length(mix.vec))
      cv.step[i,j] = which.min(apply(tempmat[,i,j,],2,mean))

  cat("\n")
  list(cv.grid = cv.grid, cv.step = cv.step)
}

## c interface function

casparCV = function(Y,X,Ytest,Xtest,h,alpha,locs,verbose = FALSE){

  ## initialize...
  n    = nrow(X)
  nt   = nrow(Xtest)
  p    = ncol(X)
  sequ = rep(0,min(n,p+1))
  bics = rep(0,min(n,p+1))

  ## Add in locs for the intercept
  locs = cbind(2*(h+1),locs)
  locs = rbind(2*(h+1),locs)

  ##dyn.load("caspar-i-cv-dmat.so")
  ##out = .C("caspar",as.double(cbind(1,X)),as.double(Y),as.double(cbind(1,Xtest)),as.double(Ytest),as.integer(n),as.integer(p+1),as.integer(nt),as.double(locs),as.double(alpha),as.double(h),as.integer(sequ),as.double(bics),as.integer(verbose))
  out = .C("caspar",as.double(cbind(1,X)),as.double(Y),as.double(cbind(1,Xtest)),as.double(Ytest),as.integer(n),as.integer(p+1),as.integer(nt),as.double(locs),as.double(alpha),as.double(h),as.integer(sequ),as.double(bics),as.integer(verbose),PACKAGE='caspar')
  ##dyn.unload("caspar-i-cv-dmat.so")
  
  return(out)
}

## demo
## demonstrates CaSpaR

demo.caspar = function(){
  cat("\n\nWelcome to a demo of Clustered and Sparse Regression (CaSpaR)")
  cat("\nWe are going to fit CaSpaR to a simple simulated Data set...\n")
  readline("<PRESS RETURN>")
  
  ## Generate data
  n = 45
  {G = simulateCasparData(n,25,2,5,0,1,6,3); X = G$X; Y=G$Y}
  nz = which(G$truth != 0)

  ##X = matrix(rnorm(14*35),nc=35)
  ##signs = sign(runif(6,-1,1))
  ##Y = X[,1:6] %*% matrix(c(3,3,6,3,3,3)*signs,nc=1) + rnorm(14)

  ## Note that either of the below two will work for d.fun
  d.fun.demo = function(ONE,TWO) abs(ONE-TWO)
  ## d.fun.demo = abs(outer(1:ncol(X),1:ncol(X),"-"))
  kp.vec     = 2:4
  mix.vec    = (1:10)/10

  ## Fitting the model and taking the best combination of parameters
  out  = caspar.grid(X,Y,5,d.fun.demo,kp.vec = kp.vec, mix.vec = mix.vec, TRUE)
  best = which(out$cv.grid == min(out$cv.grid),arr.ind = TRUE)[1,]

  cat("\nAbove we fit a model to a data set (n=150) where there are 35 candidate predictors")
  cat("\nWe know covariates ",nz," are truly nonzero, we are testing if CaSpaR can recover these.")
  cat("\nThe CaSpaR process is completed (above was the trace) -- here are the cv scores for the grid search:\n")
  print(out$cv.grid)
  readline("<PRESS RETURN>")
  cat("\n\nWe see that entry", best, "has the best cv score.")
  cat("\nThis corresponds to parameters (h =",kp.vec[best[1]],"), (alpha =",mix.vec[best[2]],")")
  cat("\nIn order to fit the model properly, we then check which step had the best score")
  cat("\nfor this combination of parameters (below is this for all combinations we tried):\n")
  print(out$cv.step)
  cat("\n\nWe now refit the model at that number of parameters...\n")
  readline("<PRESS RETURN>")

  ## Fit caspar for the best combination of parameters...
  take = out$cv.step[best[1],best[2]]
  out2 = caspar.engine(Y,X,NULL,NULL,kp = kp.vec[best[1]], mix = mix.vec[best[2]],d.fun.demo,trace=FALSE)
  
  cat("\nWe then take the first",take," selected predictors, as suggested above:\n")
  cat(out2$sequ[1:take])
  cat("\n\nWe now refit the model for these predictors with lm()...\n")
  readline("<PRESS RETURN>")

  ## Fit the final linear model...
  Covariate = X[,sort(out2$sequ[1:take])]
  FL = function(vec)
    vec[c(length(vec),1:(length(vec)-1))]
  strReverse <- function(x)
    sapply(lapply(strsplit(x, NULL), rev), paste, collapse="")
  strLF <- function(x)
    sapply(lapply(strsplit(x, NULL), FL), paste, collapse="")
  colnames(Covariate) = strLF(sapply(sort(out2$sequ[1:take]),paste,"X ",sep=""))
  print(summary(lm(Y~Covariate)))

  cat("Note that we have selected",length(intersect(out2$sequ[1:take],nz)),"correct predictors of", length(nz),".\n")
}

## Simulator
## demo data
## {G = simulate.data(n,250,7,5,0,1,6,3); X = G$X; Y=G$Y}
simulateCasparData = function(n,p,cl,sz,cor,sig,peaks, valleys, binom = FALSE){
  Y = rep(0, n)

  valid = TRUE  # Wait for separated clusters

  ## place the centers of the clusters
  ## we want separated clusters...
  while(valid){
    centers  = sample(1:p,cl,replace = FALSE)
    outlying = numeric()
    clusters = list()
    
    for(i in 1:cl){
      maxplace = min(centers[i],sz)
      minplace = max(centers[i]-p+sz,1)
            
      where = sample(rep(minplace:maxplace,each=2),1,replace = FALSE)
      outlying = c(outlying,centers[i]:(centers[i]+sz-1) - where + 1)
            
      clusters[[i]] = setdiff(centers[i]:(centers[i]+sz-1) - where + 1,centers[i])
    }

    if(length(outlying) == length(unique(outlying))) valid = FALSE
    else valid = TRUE
    
    outlying = unique(outlying)
    outlying = setdiff(outlying,centers)
  }

  truevec           = rep(0,p)
  truevec[centers]  = (2*rbinom(cl,1,.5)-1)*rep(peaks,cl)
  truevec[outlying] = (2*rbinom((sz-1)*cl,1,.5)-1)*rep(valleys,(sz-1)*cl)

  X = matrix(rnorm(n*p), nr = n,nc = p)
  the.covmat = matrix(cor,nr = sz, nc = sz) + diag(1-cor,sz)
  for(i in 1:cl) X[,c(centers[i],clusters[[i]])] = rmvnorm(n, rep(0,sz), the.covmat)
  
  colnames(X) = seq(1:p)
  Y = X %*% truevec + rnorm(n,0,sig)
  
  return(list(X= X, Y=Y, truth = truevec, centers = centers) )
}


## hiv demo
## performs HIV data analysis...

demo.hiv = function(res = 1){
  resNumber = match(res,1:7,1)
  resName   = c("APV","ATV","IDV","LPV","NFV","RTV","SQV")[resNumber]
  cat("\nData file for this drug: ",paste("HIVDAT",resNumber,resName,".rda",sep=""))

  ## load the data files...
  if(resNumber == 1){
    data(HIVDAT1APV)
    resData = HIVDAT1APV
  }
  if(resNumber == 2){
    data(HIVDAT2ATV)
    resData = HIVDAT2ATV
  }
  if(resNumber == 3){
    data(HIVDAT3IDV)
    resData = HIVDAT3IDV
  }
  if(resNumber == 4){
    data(HIVDAT4LPV)
    resData = HIVDAT4LPV
  }
  if(resNumber == 5){
    data(HIVDAT5NFV)
    resData = HIVDAT5NFV
  }
  if(resNumber == 6){
    data(HIVDAT6RTV)
    resData = HIVDAT6RTV
  }
  if(resNumber == 7){
    data(HIVDAT7SQV)
    resData = HIVDAT7SQV
  }
  geo       = resData$geo
  grp       = resData$grp
  X         = resData$iX
  Y         = resData$Y
  d.fun.hiv = function(ONE,TWO) return(abs(cov2grid(col2cov(ONE,grp),geo) - cov2grid(col2cov(TWO,grp),geo)))

  caspar.main(X,Y,d.fun.hiv,2,0.1,5,TRUE)
}


## Some HIV distance helper functions
## converts from column space to covariate space
col2cov = function(k, grp){
  j = min(which(cumsum(grp)>=k))
  return(j)
}

## converts from covariate to grid space
cov2grid = function(k,geo){
  if(k == 0) k = sum(1-geo)
  j = which(geo == 0)[k]
  return(j)
}



###################
## A function which does it all...

caspar.main = function(X,Y,distanceFunction,kernelVector,mixVector,cvFolds,trace=FALSE){
  ## argument checking
  call = match.call()
  ## CHECK FOR AIC, BIC
  gridOut = caspar.grid(X,Y,cvFolds,distanceFunction,kp.vec = kernelVector, mix.vec = mixVector, trace = trace)
  best    = which(gridOut$cv.grid == min(gridOut$cv.grid),arr.ind = TRUE)[1,]

  X = as.matrix(X)
  Y = unlist(Y)
  
  if(trace){
    cat("\n\nWe see that entry", best, "has the best cv score.")
    cat("\nThis corresponds to parameters (h =",kernelVector[best[1]],"), (alpha =",mixVector[best[2]],")")
  }

  ## fit caspar at the best parameters
  take    = gridOut$cv.step[best[1],best[2]]
  bestOut = caspar.engine(Y,X,NULL,NULL,kp = kernelVector[best[1]], mix = mixVector[best[2]],distanceFunction,trace=trace,eval="BIC")

  finalCovariate = X[,sort(bestOut$sequ[1:take])]

  FL = function(vec)
    vec[c(length(vec),1:(length(vec)-1))]
  strReverse <- function(x)
    sapply(lapply(strsplit(x, NULL), rev), paste, collapse="")
  strLF <- function(x)
    sapply(lapply(strsplit(x, NULL), FL), paste, collapse="")

  colnames(finalCovariate) = strLF(sapply(sort(bestOut$sequ[1:take]),paste,"X ",sep=""))

  finalLM = lm(Y~finalCovariate)
  if(trace) print(summary(lm(Y~finalCovariate)))

  list(finalModel = finalLM, grid = gridOut, call = call)
}
