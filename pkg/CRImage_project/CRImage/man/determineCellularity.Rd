\name{determineCellularity}
\alias{determineCellularity}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
determineCellularity(classes, classifiedCells, dimImg, img, imgW, indexWhitePixel, cancerIdentifier)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{classes}{ ~~Describe \code{classes} here~~ }
  \item{classifiedCells}{ ~~Describe \code{classifiedCells} here~~ }
  \item{dimImg}{ ~~Describe \code{dimImg} here~~ }
  \item{img}{ ~~Describe \code{img} here~~ }
  \item{imgW}{ ~~Describe \code{imgW} here~~ }
  \item{indexWhitePixel}{ ~~Describe \code{indexWhitePixel} here~~ }
  \item{cancerIdentifier}{ ~~Describe \code{cancerIdentifier} here~~ }
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
function(classes,classifiedCells,dimImg,img,imgW,indexWhitePixel,cancerIdentifier){
	print(classes)
	classValues = unique(classes)
	print(classValues)
	wholeCellDensityImage=img
	print("determine cell density")
	imgTC=imgW
	imgW[img[,,1]==2]=-1
	rm(img)
	imgW[indexWhitePixel]=-1
	predictedClassesN=c(as.character(classes),cancerIdentifier)
	imgTC[imgTC==0]=length(predictedClassesN)
	a=array(predictedClassesN[imgTC] != cancerIdentifier,dim(imgW))
	imgTC[a]=0
	imgTC[imgTC==length(predictedClassesN)]=0
	
	print("imgTC created")
	classifiedCells[,"g.x"]=as.numeric(as.character(classifiedCells[,"g.x"]))
	classifiedCells[,"g.y"]=as.numeric(as.character(classifiedCells[,"g.y"]))
	numWindows=4
	xStepSize=dimImg[1]/numWindows
	yStepSize=dimImg[2]/numWindows
	xl=1
	xr=xStepSize
	yo=1
	yu=yStepSize
	cellMatrix=data.frame(t(c(rep(NA,length(classValues)),NA,NA,NA,NA,NA,NA)))
	colnames(cellMatrix)=c(classValues,"positivePixels","cancerPixels","ratioCancerPixels","ratioCancerCellPixel","densityCells","numberNeighbors")
	
	
	print("values for small Image")
	
	hColors=col2rgb(heat.colors(50))
	hColors=hColors[,dim(hColors)[2]:1]
	rColors=hColors[1,]/256
	gColors=hColors[2,]/256
	bColors=hColors[3,]/256
	
	
	
	for (i in 1:numWindows){
		xl=1
		xr=xStepSize
		print("yu yo")
		for (j in 1:numWindows){
			windowX=(xr-xStepSize)/2
			windowY=(yu-yStepSize)/2
		
			ptm1 <- proc.time()
			imgSubCells=imgTC[xl:xr,yo:yu]
		
			t1=proc.time()-ptm1
			
			ptm2 <- proc.time()
			imgSub1=imgW[xl:xr,yo:yu]
			t2=proc.time()-ptm2
			
			ptm3 <- proc.time()
			cellPixelImage=imgSubCells
			cellPixelImage[,]=0
			numPositivePixels=length(imgSub1[imgSub1 != -1])
			t3=proc.time()-ptm3
			
			ptm4 <- proc.time()
			numCancerPixels=length(imgSubCells[imgSubCells != 0])
			
			ratioCancerPixel=numCancerPixels/numPositivePixels
			if(is.na(ratioCancerPixel)){
				ratioCancerPixel=0
			}
			cellsWindow=c()
			for (classValue in classValues){
				cells=subset(classifiedCells,classifiedCells[,"g.x"]<=xr & classifiedCells[,"g.x"]>=xl & classifiedCells[,"g.y"]<=yu & classifiedCells[,"g.y"]>=yo & classes==classValue)
				print(classValue)
				if(classValue == cancerIdentifier){
					cancerCells=cells
				}
				cellsWindow=c(cellsWindow,dim(cells)[1])
			}
			
			
			
			ratioCancerCellPixel=dim(cancerCells)[1]/numPositivePixels
			if(is.na(ratioCancerCellPixel)){
				ratioCancerCellPixel=0
			}
			colorRatioCancerCellPixel=(ratioCancerCellPixel*500)*50+1
			if(colorRatioCancerCellPixel>length(rColors)){
				colorRatioCancerCellPixel=length(rColors)
			}
			if(is.na(ratioCancerCellPixel)){
				wholeCellDensityImage[xl:xr,yo:yu,1]=1
				wholeCellDensityImage[xl:xr,yo:yu,2]=1
				wholeCellDensityImage[xl:xr,yo:yu,3]=1
			}else{
				wholeCellDensityImage[xl:xr,yo:yu,1]=gColors[colorRatioCancerCellPixel]
				wholeCellDensityImage[xl:xr,yo:yu,2]=rColors[colorRatioCancerCellPixel]
				wholeCellDensityImage[xl:xr,yo:yu,3]=bColors[colorRatioCancerCellPixel]
			}
		
			cancerDensity=mean(cancerCells[,"densityValues"])
			cancerNeighbors=mean(cancerCells[,"numberNeighbors"])
			
			
			if(is.na(cancerDensity)){
				cancerDensity=0
			}
			if(is.na(cancerNeighbors)){
				cancerNeighbors=0
			}
			
			
			cellsWindow=c(cellsWindow,numPositivePixels,numCancerPixels,numCancerPixels/numPositivePixels,dim(cancerCells)[1]/numPositivePixels,cancerDensity,cancerNeighbors,dim(imgTC)[1],dim(imgTC[2]))
			cellMatrix=rbind(cellMatrix,cellsWindow)
			xl=xl+xStepSize
			xr=xr+xStepSize
			
		}
		yo=yo+yStepSize
		yu=yu+yStepSize
	}

	print("createList")
	wholeCellDensityImage=Image(wholeCellDensityImage)
	colorMode(wholeCellDensityImage)=Color
	l=list(cellMatrix,wholeCellDensityImage)
	

  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ misc }
