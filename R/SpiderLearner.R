# Implementation file for SpiderLearner class
library(corrplot)
library(doParallel)
library(igraph)
library(foreach)
library(Rsolnp)
library(R6)

Candidate = R6::R6Class(
  private = list(
    .identifier = ""
  ), # end private
  public = list(

    initialize = function(identifier)
    {
      private$.identifier = identifier;
    },

    # gets and sets
    setIdentifier = function(identifier){
      private$.identifier = identifier;
      invisible(self)},
    getIdentifier = function(identifier){return(private$.identifier)},

    # fit function
    # in: n x p matrix with n samples, p predictors
    # returns the estimated precision matrix
    fit = function(data)
    {
      result = 0*diag(ncol(data)) # vanilla option returns empty precision matrix
      return(result)
    }
  ) # end public
)

SpiderLearner = R6::R6Class(
  private = list(
    .library = list(), # list of Candidates
    .K = 5,
    .nCores = 1,
    .calcLoss = function(alphas,nets,dataTest)
    {
      testMat = 0*diag(nrow(nets[[1]]))
      for(i in 1:length(nets))
      {
        testMat = testMat + alphas[i]*nets[[i]]
      }

      return(-1*private$.loglikLossfunction(testMat,dataTest))
    },
    .estimateKNetworks = function(data,nCores=private$.nCores)
    {

      n=nrow(data)
      K=private$.K

      # split the data into K folds and fit the model in each one
      folds = c(rep(1:K,floor(n/K)))
      if(n %% K != 0) folds = c(folds,c(1:(n-K*floor(n/K))))
      folds = sample(folds)
      foldsDF = data.frame("sample"=paste("Sample",1:n),"fold"=folds)

      doParallel::registerDoParallel(private$.nCores)
      foldsNets = list()
      foldsNets = foreach(k=1:K) %dopar%
        {
	  print(paste("[SpiderLearner] Estimating in fold",k))
          foreach(i=1:length(private$.library))%do%
            {
	      print(paste("[SpiderLearner] Fitting model"), private$.library[[i]].getIdentifier())
              tmp <- private$.library[[i]]$fit(data[foldsDF$fold !=k,])
	      tmp
            }

        }

      doParallel::stopImplicitCluster()

      # return the data with fold assignments as well as the estimated networks
      return(list(foldsDF,foldsNets))
    },
    .loglikLossfunction = function(thetaEst,dataTest)
    {
      n = nrow(dataTest)
      p = ncol(dataTest)
      firstTerm = -n/2*log(det(solve(thetaEst)))
      secondTerm= 0
      for(i in 1:n)
        secondTerm = secondTerm - 1/2*t(dataTest[i,]) %*% thetaEst %*% dataTest[i,]
      return(firstTerm + secondTerm)
    },
    .objectiveFunction = function(alphas,foldsNets,data,foldsDF)
    {
      loglikResults = matrix(rep(NA,private$.K),ncol=1)
      for(k in 1:private$.K)
      {
        dataTest = data[foldsDF$fold == k,]
        print(paste("[SpiderLearner] Calculating in fold", k))
        loglikResults[k] = 1/nrow(dataTest)*private$.calcLoss(alphas,foldsNets[[k]],
                                                              dataTest)
      } # end of inner k-fold loop

      return(mean(loglikResults[,1])) # mean cross-validated loss
    },

    .boundedObjectiveFunction = function(alphas,foldsNets,data,foldsDF)
    {
      expitLossResults = matrix(rep(NA,private$.K),ncol=1)
      for(k in 1:private$.K)
      {
        dataTest = data[foldsDF$fold == k,]
        print(paste("[SpiderLearner] Calculating in fold", k))
        negll = 1/ncol(dataTest)*1/nrow(dataTest)*private$.calcLoss(alphas,foldsNets[[k]],
                                                                    dataTest)
        expitLossResults[k] = exp(negll)/(1+exp(negll))
      } # end of inner k-fold loop

      return(mean(expitLossResults[,1])) # mean cross-validated loss
    }


  ), # end private

  public = list(

    initialize = function(){},
    addCandidate = function(candidate){
      private$.library[[length(private$.library)+1]] = candidate;
      invisible(self);
    },
    printLibrary = function()
    {
      print(length(private$.library))
      for(i in 1:length(private$.library))
        print(private$.library[[i]]$getIdentifier())
    },
    removeCandidate = function(identifier)
    {
      candidateNames = sapply(private$.library,function(x) x$getIdentifier());
      delIndex = which(candidateNames == identifier)
      private$.library[[delIndex]] = NULL
      invisible(self)
    },
    runBootstrap = function(data,K=5,standardize=T,nBoot=1,nBootCores=1)
    {
      bootArray = array(rep(NA,nBoot*ncol(data)*ncol(data)),
                        dim=c(nBoot,ncol(data),ncol(data)))
      bootWeights = matrix(rep(NA,nBoot*length(private$.library)),ncol=length(private$.library))

      # convert to parallel processing
      library(foreach)
      library(doParallel)

      doParallel::registerDoParallel(nBootCores)
      bootModels = foreach(b=1:nBoot) %dopar%
        {
          thisData = data[sample(nrow(data),replace=T),]
          thisEst = self$runSpiderLearner(data,K,standardize)
          #bootArray[b,,] = thisEst$optTheta
          #bootWeights[b,] = thisEst$weights
        }
      doParallel::stopImplicitCluster()

      for(b in 1:nBoot)
      {
        thisEst = bootModels[[b]]
        bootArray[b,,] = thisEst$optTheta
        bootWeights[b,] = thisEst$weights
      }
      return(list(bootArray, bootWeights))

    },
    runSpiderLearner = function(data,K=5,standardize=T,nCores=1,boundedLoss=F)
    {
      private$.K = K
      private$.nCores = nCores

      if(standardize)
      {
        data = apply(data,2,function(x){return((x-mean(x))/sd(x))})
      }

      foldEstimates = private$.estimateKNetworks(data,nCores)
      foldsNets = foldEstimates[[2]]

      # find optimal coefficients
      equal = function(alphas,foldsNets,data,foldsDF){return(sum(alphas))}
      nMod = length(private$.library)
      if(!boundedLoss)
      {
        alphaOpt = Rsolnp::solnp(rep(1/nMod,nMod),
                         fun=private$.objectiveFunction,
                         eqfun=equal,
                         eqB=1,
                         foldsNets=foldsNets,
                         data=data,
                         foldsDF=foldEstimates[[1]],
                         LB=rep(0,nMod),UB=rep(1,nMod))
      }

      if(boundedLoss)
      {
        alphaOpt = Rsolnp::solnp(rep(1/nMod,nMod),
                         fun=private$.boundedObjectiveFunction,
                         eqfun=equal,
                         eqB=1,
                         foldsNets=foldsNets,
                         data=data,
                         foldsDF=foldEstimates[[1]],
                         LB=rep(0,nMod),UB=rep(1,nMod))
      }

      fullModels = list()
      for(i in 1:nMod){
        fullModels[[i]] = private$.library[[i]]$fit(data)
      }

      thetaEstOpt = 0*diag(ncol(data))
      for(i in 1:nMod){
        thetaEstOpt = thetaEstOpt + alphaOpt$pars[i]*fullModels[[i]]
      }

      # calculate a simple mean model for benchmarking
      thetaSimpleMean = 1/nMod*Reduce('+', fullModels)

      return(list("optTheta"=thetaEstOpt,
                  "weights"=alphaOpt$pars,
                  "simpleMeanNetwork"=thetaSimpleMean,
                  "foldsNets"=foldsNets,
                  "fullModels"=fullModels))
    }
    # gets and sets

  ) # end public
)


