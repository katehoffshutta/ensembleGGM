# Each candidate inherits from the Candidate class defined in spiderLearnerOOP.R
library(clime)
library(hglasso)
library(huge)
library(qgraph)

HugeEBICCandidate = R6::R6Class(

  inherit = Candidate,

  private = list(
    .identifier = "",
    .gamma = 0.5
  ), # end private

  public = list(

    initialize = function(gamma)
    {
      private$.identifier = paste("ebic",gamma,sep="_");
      private$.gamma = gamma
    },

    fit = function(trainData,criterion)
    {
      mod = huge::huge(trainData,method = "glasso",verbose=F)
      modOpt = huge::huge.select(mod,
                           criterion = "ebic",
                           ebic.gamma = private$.gamma,
			   verbose=F)
      return(modOpt$opt.icov)
    }
  ) # end public
)

HugeRICCandidate = R6::R6Class(

  inherit = Candidate,

  private = list(
    .identifier = "ric"
  ), # end private

  public = list(

    initialize = function(){},

    fit = function(trainData,criterion)
    {
      mod = huge::huge(trainData,method = "glasso",verbose=F)
      modOpt = huge::huge.select(mod,
                           criterion = "ric",
			   verbose=F)
      # return(modOpt$opt.icov)
      # there were some issues with a null return here, so hand coding a refit of the model with optimal lambda
      modFin = huge::huge(trainData, method = "glasso",lambda = modOpt$opt.lambda,verbose=F)
      return(modFin$icov[[1]]) # There is only one lambda, so there will only be one entry in modFin$icov
    }
  ) # end public
)

HugeStARSCandidate = R6::R6Class(

  inherit = Candidate,

  private = list(
    .identifier = "",
    .thres = 0.1,
    .subsampleRatio = NULL,
    .defaultSubsampleRatio = function(n) #from the huge.select documentation
    {
      if(n<=144) return(0.8)
      return(10*sqrt(n)/n)
    }
  ), # end private

  public = list(

    initialize = function(thres,subsampleRatio=NULL)
    {
      private$.thres = thres;
      private$.identifier = paste("stars",private$.thres,sep="_");
      private$.subsampleRatio = subsampleRatio;
    },

    fit = function(trainData)
    {
      if(is.null(private$.subsampleRatio))
        private$.subsampleRatio = private$.defaultSubsampleRatio(nrow(trainData))

      mod = huge::huge(trainData,method = "glasso", verbose=F)
      modOpt = huge::huge.select(mod,
                           criterion = "stars",
                           stars.thres = private$.thres,
                           stars.subsample.ratio = private$.subsampleRatio,
			   verbose=F)
      #return(modOpt$opt.icov)
      # there were some issues with a null return here, so hand coding a refit of the model with optimal lambda
      modFin = huge::huge(trainData, method = "glasso",lambda = modOpt$opt.lambda,verbose=F)
      return(modFin$icov[[1]]) # There is only one lambda, so there will only be one entry in modFin$icov
    }
  ) # end public
)


HGlassoCandidate = R6::R6Class(
  private = list(
    .identifier = "hglasso",
    .selectTuningParams = function(data)
    {
      bestBIC = 100000

      lambda1opt = 0
      lambda2opt = 0
      lambda3opt = 0

      for(lambda1 in c(0.05,0.1,0.25,0.5))
      {
        for(lambda2 in c(0.0001,0.001,0.01,0.1,1))
        {
          for(lambda3 in c(1,2,3))
          {
            thisGraph = hglasso::hglasso(S = cov(data), lambda1 = lambda1, lambda2 = lambda2, lambda3)
            thisBIC = hglasso::hglassoBIC(thisGraph, cov(data), c=0.2)$BIC
            if(thisBIC < bestBIC)
            {
              bestBIC = thisBIC
              lambda1opt = lambda1
              lambda2opt = lambda2
              lambda3opt = lambda3
            }
          }
        }
        return(c(lambda1opt,lambda2opt, lambda3opt))
      }
    }
  ), # end private
  public = list(

    initialize = function(){},

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
      optParams = private$.selectTuningParams(data)
      optHub = hglasso::hglasso(S = cov(data), lambda1 = optParams[1], lambda2 = optParams[2], lambda3 = optParams[3])
      return(optHub$Theta)
      return(result)
    }
  ) # end public
)

MLECandidate = R6::R6Class(

  inherit = Candidate,

  private = list(
    .identifier = "mle"
  ), # end private

  public = list(

    initialize = function(){},

    fit = function(trainData)
    {
      return(solve(cov(trainData)))
    }
  ) # end public
)

QGraphEBICCandidate = R6::R6Class(

  inherit = Candidate,

  private = list(
    .identifier = "",
    .gamma = 0.5
  ), # end private

  public = list(

    initialize = function(gamma)
    {
      private$.identifier = paste("qgraph_ebic",gamma,sep="_");
      private$.gamma = gamma
    },

    fit = function(trainData)
    {
      mod = qgraph::EBICglasso(S=cov(trainData),n=nrow(trainData),
                               lambda.min.ratio=0.01,
                               gamma = private$.gamma,
                               returnAllResults=T, checkPD=FALSE,
			       verbose=F)
      # return(mod$optwi) This is occasionally null for some reason, so return the direct result instead
      return(mod$results$wi[,,which.min(mod$ebic)])
    }
  ) # end public
)

CLIMECandidate = R6::R6Class(

  inherit = Candidate,

  private = list(
    .identifier = ""
  ), # end private

  public = list(

    initialize = function(gamma)
    {
      private$.identifier = "clime";
    },

    fit = function(trainData)
    {
      climemat = clime::clime(trainData, standardize=F)
      re.cv= clime::cv.clime(climemat)
      re.clime.opt = clime::clime(trainData, standardize=F, lambda = re.cv$lambdaopt)
      return(re.clime.opt$Omegalist[[1]])
    }
  ) # end public
)
