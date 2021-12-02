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
  if(length(igraph::E(slGraph)) == 0)
  {
       plot(slGraph,
       vertex.color = "white",
       vertex.size = 20,
       vertex.label.cex = 0.5,
       layout = layout_in_circle)
  }
	
  else
  { 				
       plot(slGraph,
       edge.width = 5*abs(igraph::E(slGraph)$weight),
       edge.color = ifelse(igraph::E(slGraph)$weight > 0, "red","blue"),
       vertex.color = "white",
       vertex.size = 20,
       vertex.label.cex = 0.5,
       layout = layout_in_circle)
 }
}

PlotCandidates = function(slResult,index)
{
  slGraph = igraph::graph_from_adjacency_matrix(-cov2cor(slResult$fullModels[[index]]),
                                        diag=F,
                                        weighted=T,
                                        mode="undirected")
   if(length(igraph::E(slGraph)) == 0)
  {
       plot(slGraph,
       vertex.color = "white",
       vertex.size = 20,
       vertex.label.cex = 0.5,
       layout = layout_in_circle)
  }

  else
  {
       plot(slGraph,
       edge.width = 5*abs(igraph::E(slGraph)$weight),
       edge.color = ifelse(igraph::E(slGraph)$weight > 0, "red","blue"),
       vertex.color = "white",
       vertex.size = 20,
       vertex.label.cex = 0.5,
       layout = layout_in_circle)
 }
}
