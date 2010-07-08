\name{calculateCellularity}
\alias{calculateCellularity}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
calculateCellularity(filename = "", image = NA, classifier, cancerIdentifier, KS=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{filename}{A path to an image file.}
  \item{image}{If filename is undefined, an Image object }
  \item{classifier}{A SVM object, created with createClassifier or directly with the package e1071}
  \item{cancerIdentifier}{ A string which describes, how the cancer class is named.}
  \item{KS}{Apply kernel smoother?}
}
\details{
  The method calculates tumour cellularity of an image. The cells of the image are classified and the cellularity is: numTumourCells/numPixel.
  Furthermore the number of cells of the different classes are counted. A heatmap of cellularity is created. The image is divided in 16 subwindows and cellularity is  calculated for every subwindow. Green in the heatmaps indicates strong cellularity, white low cellularity.
}
\value{
  A list containing
  \item{cellulairty values }{a vector, the n first values indicate the n numbers of cells in the n classes, the n+1 value indicates the tumour cellularity}
  \item{comp2 }{a 'Image' object}
  ...
}
\references{ ~put references to the literature/web site here ~ }
\author{ Henrik Failmezger }
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
#calculation of cellularity
f = system.file("data", "exImg.jpg", package="CRImage")
exImg=readImage(f)
cellularityValues=calculateCellularity(f,classifier=classifier,cancerIdentifier="1")





## The function is currently defined as
function(filename="",image=NA,classifier,cancerIdentifier){
	print(filename)
	print(cancerIdentifier)
    if(filename !=""){
		imageData=try(segmentImage(filename=filename,image=NA))
	}else{
		imageData=try(segmentImage(filename="",image=image))
	}
	print(inherits(imageData,"try-error"))
	img=imageData[[1]]
	imgW=imageData[[2]]
	cellData=imageData[[3]]
	indexWhitePixel=imageData[[4]]
	print(indexWhitePixel)
	classes=classifyCells(classifier,"",img,imW,cellData,FALSE)
	cellularity=try(determineCellularity(classes[[1]],cellData,dim(imgW),img,imgW,indexWhitePixel,cancerIdentifier))
	cellularity
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ misc }