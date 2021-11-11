\name{PlotSpiderLearner}
\alias{PlotSpiderLearner}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot SpiderLearner Results
}
\description{
This function generates a simple igraph plot of the ensemble GGM, using the results from the runSpiderLearner() member function of the SpiderLearner object.
}
\usage{
PlotSpiderLearner(slResult)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{slResult}{
A list returned by the runSpiderLearner() member function of the SpiderLearner object.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
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
data(USJudgeRatings)
mySL = MakeSpiderLearner()
judgeMod = mySL$runSpiderLearner(data = USJudgeRatings,K=2,boundedLoss = T,standardize = T,nCores = 1)

PlotSpiderLearner(judgeMod)
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
