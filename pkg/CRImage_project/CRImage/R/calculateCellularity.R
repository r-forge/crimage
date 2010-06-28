`calculateCellularity` <-
function(filename="",image=NA,classifier=NA,cancerIdentifier=NA,KS=FALSE){
	print(filename)
	print(cancerIdentifier)
	print(classifier)
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
	classes=try(classifyCells(classifier,"",img,imgW,cellData,paint=TRUE,KS=KS))
	print(inherits(classes,"try-error"))
	print(classes)
	cellularity=try(determineCellularity(classes[[1]],cellData,dim(imgW),img,imgW,indexWhitePixel,cancerIdentifier))
	list(cellularity[[1]],cellularity[[2]],classes[[2]])
}

