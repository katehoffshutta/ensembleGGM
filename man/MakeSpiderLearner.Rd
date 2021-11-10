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
\details{
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

mySL = MakeSpiderLearner()
lawyerMod = mySL$runSpiderLearner(data = USJudgeRatings,K=2,boundedLoss = T,standardize = T,nCores = 1)
corrplot(-cov2cor(lawyerMod$optTheta))

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
