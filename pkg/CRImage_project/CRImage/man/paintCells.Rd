\name{paintCells}
\alias{paintCells}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
paintCells(imgWT, img, classes, index, classValues)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{imgWT}{ ~~Describe \code{imgWT} here~~ }
  \item{img}{ ~~Describe \code{img} here~~ }
  \item{classes}{ ~~Describe \code{classes} here~~ }
  \item{index}{ ~~Describe \code{index} here~~ }
  \item{classValues}{ ~~Describe \code{classValues} here~~ }
}
\details{
  ~~ If necessary, more details than the description above ~~
}
\value{
  ~Describe the value returned
  If it is a LIST, use
  \item{comp1 }{Description of 'comp1'}
  \item{comp2 }{Description of 'comp2'}
  ...
}
\references{ ~put references to the literature/web site here ~ }
\author{ ~~who you are~~ }
\note{ ~~further notes~~ 

 ~Make other sections like Warning with \section{Warning }{....} ~
}
\seealso{ ~~objects to See Also as \code{\link{help}}, ~~~ }
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function(imgWT,img,classes,index,classValues){
	indexClass=cbind(index,classes)
	indexClass=rep("n",max(imgWT))
	indexClass[index]=classes
	indexClass=c(indexClass,"n")
	imgOld=img
	cols=c("red","green","yellow","orange","blue")
	counter=1
	for (c in classValues){
		imgTC=imgWT
		imgTC[imgTC==0]=length(indexClass)
		a=array(indexClass[imgTC] !=c ,dim(imgWT))
		imgTC[a]=0
		imgTC[imgTC==length(indexClass)]=0
		imgOld = paintObjects(imgTC,imgOld,col=cols[counter])
		counter=counter+1
	}
	
	imgOld
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{misc}
