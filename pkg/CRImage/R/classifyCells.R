`classifyCells` <-
function(classifier,filename="",image=NA,segmentedImage=NA,featuresObjects=NA,paint=TRUE,KS=FALSE,cancerIdentifier=NA){
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
			print("classify with topological features")
			index=NULL
			class=NULL
			g.x=NULL
			g.y=NULL
			g.edge=NULL
			featuresObjectsS=subset(featuresObjects, select = -c(index,class,g.x,g.y,g.edge))
			predictedClasses=try(predict(classifier,featuresObjectsS, probability = TRUE))
	}
	classSVM=predictedClasses[1:dim(featuresObjects)[1]]
	if(KS==TRUE){
		classLabels=unique(as.character(classSVM))
		if(length(classLabels) >2 && is.na(cancerIdentifier)){
			classValues=classifier$levels
			print("More than two classes, Kernel Smoother is not applied.")
		}else if(length(classLabels) >2 && is.na(cancerIdentifier)==FALSE){
			print("Apply Kernel Smoother")
			classSVM=kernelSmoother(predictedClasses,cellCoordinates,segmentedImage,cancerIdentifier,indexCells)
			classValues=sort(unique(classSVM))
		}else{
			print("Apply Kernel Smoother")
			classSVM=kernelSmoother(predictedClasses,cellCoordinates,segmentedImage,1,indexCells)
			classValues=sort(unique(classSVM))
		}
	}else{
		classValues=classifier$levels
	}
		print(classValues)
		if(paint==TRUE){
			paintedCells=paintCells(segmentedImage,image,as.character(classSVM),indexCells,classValues)
			list(as.character(classSVM),paintedCells,cbind(indexCells,cellCoordinates,as.character(classSVM)))
		}else{
			list(as.character(classSVM),cbind(indexCells,cellCoordinates,as.character(classSVM)))
		}

}

