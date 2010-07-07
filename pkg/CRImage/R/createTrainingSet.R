`createTrainingSet` <-
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

