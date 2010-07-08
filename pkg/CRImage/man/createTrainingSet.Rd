\name{createTrainingSet}
\alias{createTrainingSet}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  Creates a training set for cell classification.
}
\usage{
createTrainingSet(filename = "", image = NA)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{filename}{ Path to an image file. }
  \item{image}{An 'Image' object, if filename is not specified. }
}
\details{
  The image is segmented. A new image is created with numbers for every cell. A Table with the features of the cells is created,too.
  To create the training set, the table with the cell featurs has to be opened for example with a spreadsheet program. Class values for the cells
  have to be inserted in the column 'class'. The corresponding cell in the image can be identified by the column 'index' (numbers in column index correspond to numbers
  in the image). Class values for different classes can be numbers or strings.
  Be careful, this function does not work on MacOsX.
}
\value{
  A List containing:
  \item{comp1 }{Image with numbers.}
  \item{comp2 }{Table of the cell features.}
}
\author{ Henrik Failmezger }

}
\seealso{ 'createClassifier' }
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function(filename="",image=NA){
	if(filename != ""){
		imageData=try(segmentImage(filename=filename))
	}else{
		imageData=try(segmentImage(image=image))
	}
	print(inherits(imageData,"try-error"))
	featuresObjects=imageData[[3]]
	image=imageData[[1]]
	segmentedImage=imageData[[2]]
	res=paintObjects(segmentedImage,image,col="green")
	nuctext=featuresObjects[,"index"]
	xy=featuresObjects[,c("g.x","g.y")]
	font = drawfont(weight=600, size=10)
	text=drawtext(res,xy=as.matrix(xy),labels=as.character(nuctext),font=font,col="yellow")
	l=list(text,featuresObjects)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ misc }
