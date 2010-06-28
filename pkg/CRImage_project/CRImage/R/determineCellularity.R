`determineCellularity` <-
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
	print("values for small Image")
	hColors=col2rgb(heat.colors(50))
	hColors=hColors[,dim(hColors)[2]:1]
	rColors=hColors[1,]/256
	gColors=hColors[2,]/256
	bColors=hColors[3,]/256
	
	print("num postive pixels")
	print(classes)
	cellularity=length(classes[classes==cancerIdentifier])/length(imgW[imgW != -1])
	cellularityValues=c()
	numClassCells=c()
	cancerCells=0
	for (classValue in classValues){
				numClassCells=c(numClassCells,length(classes[classes==classValue]))
				if(classValue == cancerIdentifier){
					cancerCells=length(classes[classes==classValue])
				}
	}
	indexValues=c(classValues,"cellularity","ratioTumourCells")
	print(numClassCells)
	
	ratioCancerCells=cancerCells/sum(numClassCells)
	cellularityValues=c(cellularityValues,numClassCells,cellularity,ratioCancerCells)
	names(cellularityValues)=indexValues
	cancerCells=c()
	for (i in 1:numWindows){
		xl=1
		xr=xStepSize
		print("yu yo")
		for (j in 1:numWindows){
			imgSub1=imgW[xl:xr,yo:yu]
			numPositivePixels=length(imgSub1[imgSub1 != -1])						
			cellsWindow=c()
			cancerCells=c()
			for (classValue in classValues){
				cells=subset(classifiedCells,classifiedCells[,"g.x"]<=xr & classifiedCells[,"g.x"]>=xl & classifiedCells[,"g.y"]<=yu & classifiedCells[,"g.y"]>=yo & classes==classValue)
				print(classValue)
				if(classValue == cancerIdentifier){
					cancerCells=cells
				}
			}
			print(dim(cancerCells))
			print(cancerCells)
			if(length(cancerCells)==0){
				ratioCancerCellPixel=0
			}else{
				ratioCancerCellPixel=dim(cancerCells)[1]/numPositivePixels
			}
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
			xl=xl+xStepSize
			xr=xr+xStepSize
			
		}
		yo=yo+yStepSize
		yu=yu+yStepSize
	}

	print("createList")
	wholeCellDensityImage=Image(wholeCellDensityImage)
	colorMode(wholeCellDensityImage)=Color
	l=list(cellularityValues,wholeCellDensityImage)
	

}

