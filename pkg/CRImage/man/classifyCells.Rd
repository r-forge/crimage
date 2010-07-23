\name{classifyCells}
\alias{classifyCells}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  The function classifies cells and paints the different class types in the image.
}
\usage{
classifyCells(classifier, filename = "", image = NA, segmentedImage = NA, featuresObjects = NA, paint = TRUE, KS = FALSE,cancerIdentifier = NA)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{classifier}{ A Support Vector Machine created by createClassifier or directly by the package e1071 }
  \item{filename}{ A path to an image file. }
  \item{image}{An 'Image' object or an array.}
  \item{segmentedImage}{An 'Image' object or an array.The corresponding segmented image (created by segmentImage)}
  \item{featuresObjects}{Cell feature file of the segmentedImage (created by segmentImage)}
  \item{paint}{ If true, the classified cells are painted with different colors in the image}
  \item{KS}{ Use Kernel Smoohter in classification. }
  \item{cancerIdentifier}{A string which describes, how the cancer class is named. }
}
\details{
  The kernels smoother improves the classification for cells which are likely to occur in clusters, like tumour cells. The kernel smoothing method can only be applied for 	two classes. If there are more classes only the normal svm without kernel smoothing is applied.
  Different classes get differen colors in the image.
}
\value{
  A list with
  \item{comp1 }{classes}
  \item{comp2 }{Classes, painted in the image, if paint was true}
}

\author{Henrik Failmezger }
\note{ ~~further notes~~ 

 ~Make other sections like Warning with \section{Warning }{....} ~
}
\seealso{ ~~objects to See Also as \code{\link{help}}, ~~~ }
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

t = system.file("data", "trainingData.txt", package="CRImage")
#read training data
trainingData=read.table(t,header=TRUE)
#create classifier
classifier=createClassifier(trainingData,topo=FALSE)[[1]]
#classify cells
f = system.file("data", "exImg.jpg", package="CRImage")
classesValues=classifyCells(classifier,filename=f,KS=TRUE)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{misc}
