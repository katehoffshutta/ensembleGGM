\name{PlotCandidates}
\alias{PlotCandidates}
\title{Plot SpiderLearner Library Candidates}
\description{
This function makes a simple graph of the GGM estimated by a candidate in the SpiderLearner library.
}
\usage{
PlotSpiderLearner(slResult,index)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{slResult}{A list returned by the runSpiderLearner() member function of the SpiderLearner object.
  }
  \item{index}{The index of the candidate model to plot.
  }
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A plot.
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

# Look at the MLE
PlotCandidates(judgeMod,1)

# Look at the huge-ebic-0.5 model
PlotCandidates(judgeMod,2)

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
