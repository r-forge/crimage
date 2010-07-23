\name{createClassifier}
\alias{createClassifier}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  Creates a classifier for a training set.
}
\usage{
createClassifier(trainingData, cross = FALSE, topo = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{trainingData}{ A table, created by segmentImage with manually added classes. }
  \item{cross}{ Does 10-fold cross validation to test the classifiers performance.}
  \item{topo}{ Use topological features.}
}
\details{
  Topological features include the density of cells and the size of the surrounding cytoplasma of a cell. These features depend of the size of the image. If training image and the image to classify have different size, these features can fool the classification.
}
\value{
  A List containing:
  \item{comp1 }{The classifier}
  \item{comp2 }{cross validation performance}
  ...
}
\author{ Henrik Failmezger }
\note{ ~~further notes~~ 

 ~Make other sections like Warning with \section{Warning }{....} ~
}
\seealso{ 'createTrainingSet','classifyCells' }
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

f = system.file("data", "trainingData.txt", package="CRImage")
#read training data
trainingData=read.table(f,header=TRUE)
#create classifier
classifier=createClassifier(trainingData,topo=FALSE)[[1]]

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{misc}
