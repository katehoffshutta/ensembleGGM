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

PlotSpiderLearner = function(slResult)
{
  slGraph = igraph::graph_from_adjacency_matrix(-cov2cor(slResult$optTheta),
                                        diag=F,
                                        weighted=T,
                                        mode="undirected")
  if(length(E(slGraph)) == 0)
  {
       igraph::plot(slGraph,
       vertex.color = "white",
       vertex.size = 20,
       vertex.label.cex = 0.5,
       layout = layout_in_circle)
  }
	
  else
  { 				
       igraph::plot(slGraph,
       edge.width = 5*abs(E(slGraph)$weight),
       edge.color = ifelse(E(slGraph)$weight > 0, "red","blue"),
       vertex.color = "white",
       vertex.size = 20,
       vertex.label.cex = 0.5,
       layout = layout_in_circle)
 }
}

PlotCandidates = function(slResult,index)
{
  slGraph = graph_from_adjacency_matrix(-cov2cor(slResult$fullModels[[index]]),
                                        diag=F,
                                        weighted=T,
                                        mode="undirected")
   if(length(E(slGraph)) == 0)
  {
       igraph::plot(slGraph,
       vertex.color = "white",
       vertex.size = 20,
       vertex.label.cex = 0.5,
       layout = layout_in_circle)
  }

  else
  {
       igraph::plot(slGraph,
       edge.width = 5*abs(E(slGraph)$weight),
       edge.color = ifelse(E(slGraph)$weight > 0, "red","blue"),
       vertex.color = "white",
       vertex.size = 20,
       vertex.label.cex = 0.5,
       layout = layout_in_circle)
 }
}
