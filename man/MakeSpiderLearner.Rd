\name{MakeSpiderLearner}
\alias{MakeSpiderLearner}
\title{Instantiate SpiderLearner Object}
\description{
This function instantiates a SpiderLearner object and populates it with a library of input Candidate objects. The result can then be used to fit ensemble GGMs on input data.
}
\usage{
MakeSpiderLearner(inputCandidates)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{inputCandidates}{A list of Candidate objects.
}
}
\details{Currently supported Candidate objects are:

+ HugeEBICCandidate

+ HugeRICCandidate

+ HugeStARSCandidate

+ HGlassoCandidate

+ MLECandidate

+ QGraphEBICCandidate

The default candidates for this function are MLECandidate and HugeEBICCandidate with gamma = 0.5. Additional Candidates can be instantiated and added to the SpiderLearner object using the syntax shown in the example.

%%  ~~ If necessary, more details than the description above ~~

}
\value{
A SpiderLearner object.
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
# Fit an ensemble GGM to the judge rating in the USJudgeRatings dataset

data(USJudgeRatings)
mySL = MakeSpiderLearner()
judgeMod = mySL$runSpiderLearner(data = USJudgeRatings,K=2,boundedLoss = T,standardize = T,nCores = 1)

# corrplot of precision matrix
corrplot(judgeMod$optTheta, is.corr=F, main = "Judge Ratings Precision Matrix")

# corrplot of adjacency matrix of GGM
corrplot(-cov2cor(judgeMod$optTheta) + diag(12), main = "Judge Ratings GGM Adjacency Matrix")

# Add more candidates

ebic_gamma_0_candidate = HugeEBICCandidate$new(gamma = 0)
qgraph_gamma_0_candidate = QGraphEBICCandidate$new(gamma = 0)

mySL$addCandidate(ebic_gamma_0_candidate)
mySL$addCandidate(qgraph_gamma_0_candidate)

mySL$printLibrary()

judgeMod2 = mySL$runSpiderLearner(data = USJudgeRatings,K=2,boundedLoss = T,standardize = T,nCores = 1)

judgeMod2$weights

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
