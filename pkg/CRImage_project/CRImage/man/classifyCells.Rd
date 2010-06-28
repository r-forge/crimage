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

## The function is currently defined as
function(classifier,filename="",image=NA,segmentedImage=NA,featuresObjects=NA,paint=TRUE,KS=FALSE){
	print("classify cell")
	if(filename != ""){
		imageData=try(segmentImage(filename))
		print(inherits(imageData,"try-error"))
		featuresObjects=imageData[[3]]
		image=imageData[[1]]
		segmentedImage=imageData[[2]]
	}
	indexCells=featuresObjects[,c("index")]
	
	print("try error")
	print(inherits(imageData,"try-error"))
	cellCoordinates=featuresObjects[,c("g.x","g.y")]
	if(length(classifier$x.scale[[1]])==94){
			print("classify without topological features")
			print(length(classifier$x.scale[[1]]))
			index=NULL
			densityValues=NULL
			sizeCytoplasma=NULL
			class=NULL
			g.x=NULL
			g.y=NULL
			g.edge=NULL
			featuresObjectsS=subset(featuresObjects, select = -c(index,densityValues,sizeCytoplasma,class,g.x,g.y,g.edge))
			predictedClasses=try(predict(classifier,featuresObjectsS, probability = TRUE))
	}else{
			index=NULL
			class=NULL
			g.x=NULL
			g.y=NULL
			g.edge=NULL
			featuresObjectsS=subset(featuresObjects, select = -c(index,class,g.x,g.y,g.edge))
			print(colnames(featuresObjects))
			predictedClasses=try(predict(classifier,featuresObjects, probability = TRUE))
	}
	classSVM=predictedClasses[1:dim(featuresObjects)[1]]
	
	if(KS==TRUE){
	classLabels=unique(as.character(classSVM))
	if(length(classLabels) >2){
		print("More than two classes, Kernel Smoother is not applied.")
	}else{
		        prob=attr(predictedClasses, "prob")
				probCancer=prob[,1]
				cellCoordinatesN=cellCoordinates
				cellCoordinatesN=cbind(indexCells,cellCoordinates)
				cellCoordinatesN=cbind(cellCoordinatesN,1:dim(cellCoordinates)[1])
				imgW=dim(segmentedImage)[1]
				imgH=dim(segmentedImage)[2]
				xs=imgW/2
				xPos=imgW/2
				ys=imgH/2
				yPos=imgH/2
				allNewClasses=data.frame()
				windowCounter=0
				for (iW in 1:2){
					xPos=xs
					for (jW in 1:2){
						if(length(cellCoordinatesN)>0){
						actualCoordinates=subset(cellCoordinatesN,cellCoordinatesN[,2]<=xPos & cellCoordinatesN[,3]<=yPos)
						indexActualCoordinates=which(cellCoordinatesN[,2]<=xPos & cellCoordinatesN[,3]<=yPos)
						cellCoordinatesN=subset(cellCoordinatesN, !(cellCoordinatesN[,4] \%in\% indexActualCoordinates))
						cellCoordinatesN[,4]=1:dim(cellCoordinatesN)[1]
						if(dim(actualCoordinates)[1]>0){
							windowCounter=windowCounter+1
							distMatrix=as.matrix(dist(actualCoordinates[,2:3]))
							
							f=c()
							for (i in 1:dim(actualCoordinates)[1]){
								probCancerW=probCancer[actualCoordinates[,"indexCells"]]
								t=distMatrix[i,]/50
								indexZero=which(t>1)
								k=(1-t^3)^3
								k[indexZero]=0
								f=c(f,sum(k*probCancerW)/sum(k))
							}
							index=actualCoordinates[,"indexCells"]
							classValues=f
							newClassesKS=cbind(index,classValues)
							if(windowCounter==1){
								allNewClasses=newClassesKS
							}else{
								allNewClasses=rbind(allNewClasses,newClassesKS)
							}
						}
						}
						xPos=xPos+xs
					}
					yPos=yPos+ys
				}
			
	
			allNewClasses=as.data.frame(allNewClasses)
			print(allNewClasses)
			allNewClasses$classValues[allNewClasses$classValues<0.5]=0
			allNewClasses$classValues[allNewClasses$classValues>=0.5]=2
			allNewClasses$classValues[allNewClasses$classValues==0]=1
			allNewClasses$classValues[allNewClasses$classValues==1]=classLabels[1]
			allNewClasses$classValues[allNewClasses$classValues==2]=classLabels[2]
	
			realIndex=data.frame(indexCells)
			colnames(realIndex)="index"
			colnames(allNewClasses)=c("index","classValues")
			classesRealIndices=merge(realIndex,allNewClasses,all.x=TRUE,all.y=TRUE)
			allClasses=classesRealIndices[,"classValues"]
			print(allClasses)
			classSVM=allClasses
	 }
	}
	
	
	classValues=unique(as.character(classSVM))
	if(paint==TRUE){
		paintedCells=paintCells(segmentedImage,image,as.character(classSVM),indexCells,classValues)
		list(as.character(classSVM),paintedCells)
	}else{
		list(as.character(classSVM))
	}
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{misc}
