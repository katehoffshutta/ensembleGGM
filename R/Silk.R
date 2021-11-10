# Wrappers and helper functions for the SpiderLearner framework
source("R/SpiderLearner.R")
source("R/SpiderLearnerCandidates.R")

MakeSpiderLearner = function(candidates = list(MLECandidate$new(),HugeEBICCandidate$new(gamma = 0.5)))
{
  sl = SpiderLearner$new()

  # Add candidate for each of the candidates in the list
  for(candidate in candidates)
  {
    sl$addCandidate(candidate)
  }

  return(sl)
}


