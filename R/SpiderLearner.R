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
    .identifiers = c(), # character vector of method names
    .K = 5,
    .nCores = 1,
    .result = NULL, # when runSpiderLearner is called, fill this in with the result
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
              print(paste("[SpiderLearner] Fitting model", private$.library[[i]]$getIdentifier()))
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
      private$.identifiers = c(private$.identifiers, candidate$getIdentifier())
      invisible(self);
    },
    printLibrary = function()
    {
      print(private$.identifiers)
    },
    getGGM = function() # return adjacency matrix of GGM
    {
      if(is.null(private$.result)) return(NULL)
      return(1/2*(-cov2cor(private$.result$optTheta) + t(-cov2cor(private$.result$optTheta))) + diag(ncol(private$.result$optTheta)))
    },
    getPrecMat = function()
    {
      if(is.null(private$.result)) return(NULL)
      return(private$.result$optTheta)
    },
    getResults = function()
    {
      return(private$.result)
    },
    getWeights = function()
    {
      if(is.null(private$.result)) return(NULL)
      methods = c()
      for(i in 1:length(private$.library))
        methods = c(methods,private$.library[[i]]$getIdentifier())
      return(data.frame("method"=methods,"weight"=private$.result$weights))

    },
    plotCandidate = function(identifier,vertex.color="white",vertex.size=20,vertex.label.cex=0.5,layout=igraph::layout_in_circle)
    {
      if(is.null(private$.result)) return(NULL)
      modelIndex = which(private$.identifiers == identifier)
      slGraph = igraph::graph_from_adjacency_matrix(-cov2cor(private$.result$fullModels[[modelIndex]]),
                                                    diag=F,
                                                    weighted=T,
                                                    mode="plus")
      igraph::E(slGraph)$weight = 1/2*igraph::E(slGraph)$weight
      if(length(igraph::E(slGraph)) == 0)
      {
        plot(slGraph,
             vertex.color = vertex.color,
             vertex.size = vertex.size,
             vertex.label.cex = vertex.label.cex,
             layout = layout)
      }

      else
      {
        plot(slGraph,
             edge.width = 5*abs(igraph::E(slGraph)$weight),
             edge.color = ifelse(igraph::E(slGraph)$weight > 0, "red","blue"),
             vertex.color = vertex.color,
             vertex.size = vertex.size,
             vertex.label.cex = vertex.label.cex,
             layout = layout)
      }
      title(identifier)
    },
    plotSpiderLearner = function(vertex.color="white",vertex.size=20,vertex.label.cex=0.5,layout=igraph::layout_in_circle)
    {
      if(is.null(private$.result)) return(NULL)
      slGraph = igraph::graph_from_adjacency_matrix(-cov2cor(private$.result$optTheta),
                                                    diag=F,
                                                    weighted=T,
                                                    mode="plus")
      igraph::E(slGraph)$weight = 1/2*igraph::E(slGraph)$weight
      if(length(igraph::E(slGraph)) == 0)
      {
        plot(slGraph,
             vertex.color = vertex.color,
             vertex.size = vertex.size,
             vertex.label.cex = vertex.label.cex,
             layout = layout)
      }

      else
      {
        plot(slGraph,
             edge.width = 5*abs(igraph::E(slGraph)$weight),
             edge.color = ifelse(igraph::E(slGraph)$weight > 0, "red","blue"),
             vertex.color = vertex.color,
             vertex.size = vertex.size,
             vertex.label.cex = vertex.label.cex,
             layout = layout)
      }

      title("SpiderLearner Ensemble")
    },
    removeCandidate = function(identifier)
    {
      if(identifier %in% private$.identifiers)
      {
	      delIndex = which(private$.identifiers == identifier)
      	private$.library[[delIndex]] = NULL
      	private$.identifiers = private$.identifiers[-delIndex]
      }

      else
      {
        print(paste("[SpiderLearner::removeCandidate] No such candidate:", identifier))
      }

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
    runSpiderLearner = function(data,K=5,standardize=T,nCores=1,boundedLoss=F, ...)
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
                         LB=rep(0,nMod),UB=rep(1,nMod),
                         ...)
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
                         LB=rep(0,nMod),UB=rep(1,nMod),
                         ...)
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
      result = list("optTheta"=thetaEstOpt,
                     "weights"=alphaOpt$pars,
                     "simpleMeanNetwork"=thetaSimpleMean,
                     "foldsNets"=foldsNets,
                     "fullModels"=fullModels)
      private$.result = result # store for later helper functions
      return(result)
    }

  ) # end public
)


