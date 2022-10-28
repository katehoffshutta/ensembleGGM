\name{runSpiderLearner}
\alias{runSpiderLearner}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to run the SpiderLearner ensemble model.
}
\description{
This function is a member function of the SpiderLearner class. To use it, you first need to instantiate a SpiderLearner object using the constructor or the wrapper function `MakeSpiderLearner`. See `?MakeSpiderLearner` for details.
}
\usage{

# For a SpiderLearner object s:

s$runSpiderLearner(data,K=5,standardize=T,nCores=1,boundedLoss=F, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{Data frame with rows corresponding to samples and columns corresponding to predictors.}
  \item{K}{number of folds. Default=5}
  \item{standardize}{Standardize predictors to mean zero and standard deviation 1. Default=T}
  \item{nCores}{Number of cores for the estimation. Default=1}
  \item{boundedLoss}{Use bounded loss function. Default=F}
  \item{...}{Additional arguments passed to `Rsolnp::solnp`, the optimization function.}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
\item{optTheta}{Estimated precision matrix}
\item{weights}{Weights of each candidate in the ensemble}
\item{simpleMeanNetwork}{For comparision, a precision matrix estimated as the simple mean of the candidate estimates}
\item{foldsNets}{Precision matrices estimated within each fold of the cross-validation.}
\item{fullModels}{Precision matrices estimated with each candidate method on the full dataset, for comparision.}
}
\references{\url{https://doi.org/10.1101/2021.07.13.452248}
}
\author{
Kate Hoff Shutta (kshutta@hsph.harvard.edu)
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

# Adjust control parameters for Rsolnp::solnp
judgeMod2 = mySL$runSpiderLearner(data = USJudgeRatings,K=2,boundedLoss = T,standardize = T,nCores = 1,control=list(delta = 1e-9,tol = 1e-10))

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
