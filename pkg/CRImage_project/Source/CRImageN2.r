require(EBImage)
require(e1071)
require(MASS)

foo=function(){
print("hallo")
}

#calculates the automatic Otsu Threshold 
calculateThreshold=function(allGreyValues){
	gTrans=allGreyValues*255
	
	gTrans=round(gTrans)
	"h=hist(as.vector(gTrans))"
	t=tabulate(gTrans)
	t=t+0.001
	numPixel=length(allGreyValues)
	ni=t
	pi=ni/numPixel
	uT=function(p,L){
		res=0
		for (i in 1:L){
			res = res+i*p[i]
		}
		res
	}
	wK = function(p,k){
		res=0
		for (i in 1:k){
			res = res+p[i]
		}
		res
	}
	uK = function(p,k){
		res=0
		for (i in 1:k){
			res = res+i*p[i]
		}
		res
	}
	maxT=function(p,L,k){
		(uT(p,L)*wK(p,k)-uK(p,k))^2/(wK(p,k)*(1-wK(p,k)))
	}
	o=optimize(maxT,c(1,length(t)),maximum=TRUE,L=length(t),p=pi)
	threshold=o$maximum/255
}
segmentImage=function(filename="",image=NA,maxShape=NA,minShape=NA,failureRegion=NA){
	if(filename!=""){
		img = readImage(filename)
	}else{
		img=image
	}
	colorMode(img)="color"
	if(length(dim(img))>2){
		imgT=img[,,1]
		indexWhitePixel=which(img[,,1]>0.85 & img[,,2]>0.85 & img[,,3]>0.85)
	}else{
		print("only grayscale")
		imgT=img
		indexWhitePixel=which(img>0.85)
	}
	numWhitePixel=length(imgT[indexWhitePixel])
	numPixel=length(as.vector(imgT))
	if((numWhitePixel/numPixel)<0.9){
	print("start segmentation")
	
	print("delete dark regions")
	#search large dark regions and color them white
	imgG=channel(img,"grey")
	imgB=imgG
	imgB[imgB<0.2]=1
	imgB[imgB !=1]=0
	imgS=imgB
		
	imgS=opening(imgS,makeBrush(5, shape='disc'))
	imgS= bwlabel(imgS)
	imgSd=as.vector(imageData(imgS)+1)
	hc=tabulate(imgSd)
	print("c")
	imgSdN=imageData(imgS)+1
	rm(imgS)
	#set the failure region
	if(!is.na(failureRegion)){
		print("delete failure regions")
		a=array(hc[imgSdN]<failureRegion,dim(imgSdN))
		imgSdN[a]=-1
	}
	imgSdN=imgSdN-1
	imgR=img[,,1]
	imgGreen=img[,,2]
	imgBlue=img[,,3]
	print("mark failure regions")
	failures=which(imgSdN==-1)
	imgR[failures]=2
	imgGreen[failures]=2
	imgBlue[failures]=2
	img[,,1]=imgR
	img[,,2]=imgGreen
	img[,,3]=imgBlue
	rm(imgR)
	rm(imgBlue)
	rm(imgGreen)
	print("failure regions")
	#failure regions in thre greyscale image are colored white	
	imgG[failures]=1
	imgG[indexWhitePixel]=1
	positivePixels=dim(img[,,1])[1]*dim(img[,,1])[2]-length(indexWhitePixel)-length(failures)
		
	#local thresholding
	imgB=imgG
	numWindows=2
	xWs=round(dim(img)[1]/numWindows)
	yWs=round(dim(img)[2]/numWindows)
	xlw=1
	xrw=round(dim(img)[1]/numWindows)
	yow=1
	yuw=round(dim(img)[2]/numWindows)
	for(iW in 1:numWindows){
			xlw=1
			xrw=round(dim(img)[1]/numWindows)
			for(jW in 1:numWindows){
				if(xrw>dim(img)[1]){xrw=dim(img)[1]}
				imgGL=imgG[xlw:xrw,yow:yuw]
				if(length(as.vector(imgGL[imgGL != 1]))>0){
					globalThreshold=calculateThreshold(as.vector(imgGL[imgGL != 1]))
				}else{
					globalThreshold=0
				}
				print(paste("adaptiveThreshold",globalThreshold))
				imgBL=imgB[xlw:xrw,yow:yuw]
				imgBL[imgGL<globalThreshold]=-1
				imgBL[imgBL != -1]=0
				imgBL[imgBL==-1]=1
				imgB[xlw:xrw,yow:yuw]=imgBL
				rm(imgBL)
				rm(imgGL)
				xlw=xrw+1
				xrw=xrw+xWs
				
			}
			yow=yuw+1
			yuw=yuw+yWs
			if(yuw>dim(img)[2]){yuw=dim(img)[2]}
		}	
		#morphological opening to smooth the shape
		imgB=opening(imgB,makeBrush(5, shape='disc'))
		imgS=bwlabel(imgB)
		rm(imgB)
		
		
		#segment the image
		imgSdN=imageData(imgS)+1
		numSeq=tabulate(imageData(imgS)+1)
		#set minShape
		if(!is.na(minShape)){
			print("delete min cell nuclei")
			a=array(numSeq[imgSdN]<minShape,dim(imgSdN))
			imgSdN[a]=1
		}
		imgSdN=imgSdN-1 
		#watershed segmentation
		imgT=distmap(imgSdN)
		print("watershed")
		imgW=watershed(imgT)
		rm(imgT)
		imgSd=as.vector(imageData(imgW)+1)
		hc=tabulate(imgSd)
		rm(imgSd)
		imgSdN=imageData(imgW)+1
		if(!is.na(maxShape)){
			print("delete max cell nuclei")
			a=array(hc[imgSdN]<maxShape,dim(imgSdN))
			imgSdN[a]=1
		}
		imgSdN=imgSdN-1
			
		print("resegment large segments")
		#segment large segments another time
		imgW[imgSdN >0]=0
			
		imgLS=imgG
		imgLS[imgSdN==0]=0
		imgSdN=as.vector(imgSdN)
		imgLS=as.vector(imgLS)
		imgB_ls=imgSdN
		allSegments=unique(imgSdN)
		allSegments=allSegments[allSegments !=0]
			
		print("first segments calculated")
			
		sapply(allSegments,function(i){
				   greyPixel=imgLS[imgSdN==i]
				   indexGreyPixel=which(imgSdN==i)
				   localThreshold=calculateThreshold(as.vector(greyPixel[greyPixel !=1]))
				   thresholdedPixel=greyPixel
				   thresholdedPixel[thresholdedPixel<localThreshold]=-1
				   thresholdedPixel[thresholdedPixel!=-1]=0
				   thresholdedPixel[thresholdedPixel==-1]=1
				   imgB_ls[indexGreyPixel]<<-thresholdedPixel
				   })
		imgB_ls=array(imgB_ls,dim(imgG))
		print("opening")
		imgB_ls=opening(imgB_ls)
		print("distmap")
		imgT_ls=distmap(imgB_ls)
		print("watershed")
		imgW_ls=watershed(imgT_ls)
		rm(imgT_ls)
		maxSegment=max(imgW)
		imgW_ls=imgW_ls+maxSegment
		imgW_ls[imgW_ls==maxSegment]=0
		print("large segments segmented")
		imgW=imgW+imgW_ls
		rm(imgB_ls)
		imgSd=as.vector(imageData(imgW)+1)
		hc=tabulate(imgSd)
		imgSdN=imageData(imgW)+1
		if(!is.na(minShape)){
			a=array(hc[imgSdN]<minShape,dim(imgSdN))
			imgSdN[a]=1
		}
		imgSdN=imgSdN-1
		imgW=imgSdN
		rm(imgSdN)
		print("relabel regions")
		regions=unique(as.vector(imgW))
		regions=regions[regions != 0]
		if(length(regions>0)){
				#label the nuclei in the right way
				changeLabel=cbind(regions,1:length(regions))
				colnames(changeLabel)=c("oldL","newL")	
				imgWv=as.vector(imgW)
				imgWv_nZ=imgWv[imgWv>0]
				index_imgWv_nZ=which(imgWv>0)
				nV=rep(0,max(imgWv_nZ))
				nV[changeLabel[,1]]=changeLabel[,2]
				imgWv_nZ=nV[imgWv_nZ]
				imgWv[index_imgWv_nZ]=imgWv_nZ
				imgWv=array(imgWv,dim(imgW))
				imgW=imgWv
				rm(imgWv)
				print("calculateFeatures")
				#calculate the features
				hG=hullFeatures(imgW)
				hF=hullFeatures(imgW)
				allM=moments(imgW,imgG)
				zM=zernikeMoments(imgW,imgG)
				hFgrey=haralickFeatures(imgW,imgG)
				allFeatures=data.frame()
				class=c()
				class=rep(NA,dim(hF)[1])
				
				allFeatures=cbind(class,hF,allM,zM,hFgrey)
				allFeatures=allFeatures[allFeatures[,"g.s"]>0,]
				index=1:dim(allFeatures)[1]
				allFeatures=cbind(index,allFeatures)
				print("number neighbors")
				
				cellCoordinates=allFeatures[,c("g.x","g.y")]
				cellCoordinates[,1]=as.numeric(format(cellCoordinates[,1],digits=4))
				cellCoordinates[,2]=as.numeric(format(cellCoordinates[,2],digits=4))
				write.table(allFeatures[,c("g.x","g.y")],"cellCoordinates.txt")
				densityNeighbors=kde2d(cellCoordinates[,"g.x"],cellCoordinates[,"g.y"],n=100)
				dV=densityNeighbors$z
				dimnames(dV)=list(densityNeighbors$x,densityNeighbors$y)
				write.table(dV,"dV.txt")
				densityValues=interpolate(cellCoordinates,dV)				
				allFeatures=cbind(allFeatures,densityValues)
				
				print("size cytoplasma")
				#segment the cytoplasma
				colorMode(img)="gray"
				imgGn=as.vector(img[,,1])
				imgGn[img[,,1]==2]=-1
				imgGn[imgGn==1]=-1
				imgGn[imgW>1]=-1
				imgGn=imgGn[-indexWhitePixel]
				t=calculateThreshold(as.vector(imgGn[imgGn !=-1]))
				imgBC=imgG
				imgBC[img[,,1]<=t]=-1
				imgBC[imgBC != -1]=0
				imgBC[imgBC == -1]=1
				imgBCo=opening(imgBC,makeBrush(5, shape='disc'))
				imgBCs=bwlabel(imgBCo)
				imgBCsD=imageData(imgBCs)
				colorMode(img)="color"	
				"imgCyto=paintObjects(imgBCs,img,col=\"yellow\")"
				featuresCytoplasma=hullFeatures(imageData(imgBCs))[,"g.s"]
				imgBCs=imageData(imgBCs)
				f=function(index,imgWt,imgBCsDt){
					
					x=round(as.vector(hF[index,"g.x"]))
					y=round(as.vector(hF[index,"g.y"]))
					if(x>0 & y>0){
						segCyto=imgBCs[x,y]
						if(segCyto >0){
							s=featuresCytoplasma[segCyto]
						}else{
							0
						}
					}else{
						0
					}
				}
				
				sizeCytoplasma=sapply(index,f, imgWt=imgW,imgBCsD=as.vector(imgBCsD),simplify=TRUE)
				
				rm(imgGn)
				rm(imgBC)
				rm(imgBCo)
				
				print("number neighbors")
				
				xs=dim(img)[1]/2
				xPos=dim(img)[1]/2
				ys=dim(img)[2]/2
				yPos=dim(img)[2]/2
				cellCoordinatesN=cellCoordinates
				cellCoordinatesN=cbind(allFeatures[,"index"],cellCoordinates)
				cellCoordinatesN=cbind(cellCoordinatesN,1:dim(cellCoordinates)[1])
				print("class cell coordinates")
				#calculate the number of neighbors of every nuclei
				indexNeighbors=data.frame()
				for (i in 1:2){
					xPos=xs
					for (j in 1:2){
						if(length(cellCoordinatesN)>0){
							actualCoordinates=subset(cellCoordinatesN,cellCoordinatesN[,2]<=xPos & cellCoordinatesN[,3]<=yPos)
							indexActualCoordinates=which(cellCoordinatesN[,2]<=xPos & cellCoordinatesN[,3]<=yPos)
							cellCoordinatesN=subset(cellCoordinatesN, !(cellCoordinatesN[,4] %in% indexActualCoordinates))
							cellCoordinatesN[,4]=1:dim(cellCoordinatesN)[1]
							if(dim(actualCoordinates)[1]>0){
								distMatrix=as.matrix(dist(actualCoordinates[,2:3]))
								numberNeighbors=c()
								f=c()
								for (m in 1:dim(actualCoordinates)[1]){
									sortedCells=sort(distMatrix[m,])
									neighbors=length(sortedCells[sortedCells<50])
									numberNeighbors=c(numberNeighbors,neighbors)
								}
								numberNeighbors=cbind(actualCoordinates[,1],numberNeighbors)
								indexNeighbors=rbind(indexNeighbors,numberNeighbors)
							}
						}
						xPos=xPos+xs
					}
					yPos=yPos+ys
				}
				print("neighbors calculated")
				realIndex=data.frame(allFeatures[,"index"])
				colnames(realIndex)="index"
				colnames(indexNeighbors)=c("index","neighbors")
				neighborsRealIndices=merge(realIndex,indexNeighbors,all.x=TRUE,all.y=TRUE)
				numberNeighbors=neighborsRealIndices$neighbors
				allFeatures=cbind(allFeatures,numberNeighbors)
				allFeatures=cbind(allFeatures,sizeCytoplasma)
				
				allFeatures=allFeatures[allFeatures[,"g.s"]>0,]
				print("segmentation ended")
				l=list(img,imgW,allFeatures,indexWhitePixel,TRUE)
				
			}else{
				print("no cells")
				allFeatures=data.frame()
				l=list(img,img,allFeatures,indexWhitePixel,FALSE)
		}
		}else{
		print("almost white")
		allFeatures=data.frame()
		l=list(img,img,allFeatures,indexWhitePixel,FALSE)
	}
	
}

createTrainingSet=function(filename="",image=NA){
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

createClassifier=function(trainingData,cross=FALSE,topo=TRUE){
	trainingData = trainingData[is.na(trainingData$class)==FALSE,]
	classes=trainingData$class
	if(topo==FALSE){
		index=NULL
		densityValues=NULL
		sizeCytoplasma=NULL
		class=NULL
		g.x=NULL
		g.y=NULL
		g.edge=NULL
		trainingData=subset(trainingData, select = -c(index,densityValues,sizeCytoplasma,class,g.x,g.y,g.edge))
	}else{
	index=NULL
	class=NULL
	g.x=NULL
	g.y=NULL
	g.edge=NULL
	trainingData=subset(trainingData, select = -c(index,class,g.x,g.y,g.edge))
	} 
	model = svm(trainingData,as.character(classes),type='C',kernel='radial',probability=TRUE)
	allCrossValues=c()
	if(cross==TRUE){
		for (i in 1:10){
			crossValue = svm(trainingData,as.character(classes),type='C',kernel='radial',probability=TRUE,cross=10)
			allCrossValues=c(allCrossValues,crossValue$a)
		}
		l=list(model,mean(allCrossValues))
	}else{
		l=list(model)
	}
	
}

kernelSmoother=function(predictedClasses,cellCoordinates,segmentedImage,cancerIdentifier,indexCells){
				#options(stringsAsFactors = FALSE)
				prob=attr(predictedClasses, "prob")
				probCancer=prob[,cancerIdentifier]
				print("prob cancer")
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
						cellCoordinatesN=subset(cellCoordinatesN, !(cellCoordinatesN[,4] %in% indexActualCoordinates))
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
							newClassesKS=data.frame(index,classValues)
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
			allNewClasses$classValues[allNewClasses$classValues<0.5]=0
			allNewClasses$classValues[allNewClasses$classValues>=0.5]=2
			allNewClasses$classValues[allNewClasses$classValues==0]=1
			allNewClasses$classValues[allNewClasses$classValues==1]="o"
			allNewClasses$classValues[allNewClasses$classValues==2]=cancerIdentifier
			
			realIndex=data.frame(indexCells)
			colnames(realIndex)="index"
			colnames(allNewClasses)=c("index","classValues")
			classesRealIndices=merge(realIndex,allNewClasses,all.x=TRUE,all.y=TRUE)
			allClasses=classesRealIndices[,"classValues"]
			classSVM=allClasses
			classSVM
}
classifyCells=function(classifier,filename="",image=NA,segmentedImage=NA,featuresObjects=NA,paint=TRUE,KS=FALSE,cancerIdentifier=NA){
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
paintCells=function(imgWT,img,classes,index,classValues){
	if(length(classValues)<=10){
		indexClass=cbind(index,classes)
		indexClass=rep("n",max(imgWT))
		indexClass[index]=classes
		indexClass=c(indexClass,"n")
		imgOld=img
		cols=c("green","red","yellow","orange","blue","brown","cyan","gray","purple","violet")
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
	}else{
		print("Too many classes to paint. Do not use more than 10 classes.")
	}
}

calculateCellularity=function(filename="",image=NA,classifier=NA,cancerIdentifier=NA,KS=FALSE){
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
	classes=try(classifyCells(classifier,"",img,imgW,cellData,paint=TRUE,KS=KS))
	print(inherits(classes,"try-error"))
	classValues=classifier$levels
	cellularity=try(determineCellularity(classes[[1]],cellData,dim(imgW),img,imgW,indexWhitePixel,cancerIdentifier,classValues))
	list(cellularity[[1]],cellularity[[2]],classes[[2]])
}
determineCellularity=function(classes,classifiedCells,dimImg,img,imgW,indexWhitePixel,cancerIdentifier,classValues){
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
				if(classValue == cancerIdentifier){
					cancerCells=cells
				}
			}
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
# finds the correct section for every subimage
findSlices=function(imgFolder,pathToOutputFolder,numSlides){
	smallImage=readImage(paste(imgFolder,"SlideThumb.jpg",sep="/"))
	
	#segment the thumbnail to find the sections in the image
	pathToImgFolder=imgFolder
	imgG=channel(smallImage,"gray")
	globalThreshold=calculateThreshold(as.vector(imgG))
	imgG[imgG<globalThreshold]=-1
	imgG[imgG !=-1]=0
	imgG[imgG==-1]=1
	imgG=fillHull(imgG)
	imgG=closing(imgG,makeBrush(10,shape="diamond"))
	imgG=opening(imgG,makeBrush(10,shape="diamond"))
	imgS=bwlabel(imgG)
	hF=hullFeatures(imgS)[,c("g.x","g.y")]
	#if the image is found two times (should not happen)
	if(length(numSlides)>1){
		numSlides=numSlides[1]
	}
	print(paste("numSlides",numSlides))
	print("centers")
	
	if(numSlides==1){
		if(length(hF)==0){
			centers=matrix(c(1,1))
		}else{
		centers=matrix(c(mean(hF[,1])),c(1))
		}
	}else{
		#cluster the segments to the number of sections
		if(nrow(hF)>numSlides){
			"centers=kmeans(hF[,1],as.numeric(numSlides),iter.max=2000)$centers"
			centersIndex=cutree(hclust(dist(hF[,1])),k=as.numeric(numSlides))
			centers=matrix(0,as.numeric(numSlides))
			for(i in 1:as.numeric(numSlides)){
				print(hF[centersIndex==i,1])
				centers[i]=mean(as.numeric(as.character(hF[centersIndex==i,1])))
			}
		}else{
			if(nrow(hF) == numSlides){
				centers=matrix(as.numeric(as.character(hF[,1])),as.numeric(numSlides))
			}
		}
	}
	sliceColors=col2rgb(c("red","blue","green","yellow","orange","black","brown","purple4","darkolivegreen1","darkred","gray" ))
	
	
	finalScanInFile=paste(pathToImgFolder,"FinalScan.ini",sep="/")
	print(finalScanInFile)
	finalScanOutfileFile=paste(pathToImgFolder,"FinalScanModified.ini",sep="/")
	print(finalScanOutfileFile)
	pathToScripts = system.file("exec", "parseFinalScan.py", package="CRImage")
	print(pathToScripts)
	com=paste("python",pathToScripts,finalScanInFile,finalScanOutfileFile)
	print(com)
	system(com)
	blockPositions=read.delim(finalScanOutfileFile,header=FALSE)
	widthO=as.numeric(as.character(blockPositions[1,2]))
	heightO=as.numeric(as.character(blockPositions[1,3]))
	centerOX=widthO/2
	centerOY=heightO/2
	widthSmallI=dim(smallImage)[1]
	heightSmallI=dim(smallImage)[2]
	#draw the centers of the different clusters in the thumbnail
	for(i in 1:as.numeric(numSlides)){
		centerXSmall=centers[i,1]
		centerYSmall=round(heightSmallI/2)
		if(centerXSmall-10 >0 & centerXSmall+10<dim(smallImage)[1] & centerYSmall-10 >0 & centerYSmall+10<dim(smallImage)[2]){
		smallImage[(centerXSmall-10):(centerXSmall+10),(centerYSmall-10):(centerYSmall+10),1]=sliceColors[1,i]/255
		smallImage[(centerXSmall-10):(centerXSmall+10),(centerYSmall-10):(centerYSmall+10),2]=sliceColors[2,i]/255
		smallImage[(centerXSmall-10):(centerXSmall+10),(centerYSmall-10):(centerYSmall+10),3]=sliceColors[3,i]/255
		}
	}
	
	print(centers)
	centers=centers[sort.list(centers),]

	print(centers)
	print("hFs")
	writeImage(smallImage,paste(pathToOutputFolder,"smallImage.jpg",sep="/"))
	
	
	#assign the subimages to the sections
	blockSlice=data.frame( 2:dim(blockPositions)[1], 2:dim(blockPositions)[1], 2:dim(blockPositions)[1], 2:dim(blockPositions)[1])
	for(i in 2:dim(blockPositions)[1]){
		name=as.character(blockPositions[i,1])
		x=as.numeric(as.character(blockPositions[i,2]))
		y=as.numeric(as.character(blockPositions[i,3]))
		
		cBlockXP=x/4
		cBlockYP=y/4
		posBlockX=abs(centerOX-cBlockXP)
		posBlockY=abs(centerOY-cBlockYP)
		
		posBlockXSmall=posBlockX/(widthO/widthSmallI)
		posBlockYSmall=posBlockY/(heightO/heightSmallI)
		dMin=sqrt(widthSmallI^2+heightSmallI^2)
		slice=1
		for (k in 1:as.numeric(numSlides)){
			m=rbind(c(posBlockXSmall),centers[k])
			d=dist(m)
			if(d<dMin){
				slice=k
				dMin=d
			}
		}
		print(c(name,slice,posBlockXSmall,posBlockYSmall))
		blockSlice[(i-1),]=c(name,slice,posBlockXSmall,posBlockYSmall)
	}
	colnames(blockSlice)=c("block","slice")
	sizeO=c(widthO,heightO)
	l=list(blockSlice,sizeO)
}
# classifies an image, calls image segmentation, the classification of the cells, calculation of the cellularity and the drawing of the nuclei
classificationAperio=function(fileLocation,filename,pathToOutputFolderImgDir,classifier,pathToOutputFolderImgDirFiles,pathToOutputFolderImgDirImages,pathToOutputFolderImgDirCellDensity,blockSlice,sliceColors,sizeO,index,cancerIdentifier,pathToOutputFolderImgDirCellCoordinates){
				widthO=as.numeric(as.character(sizeO[1]))
				heightO=as.numeric(as.character(sizeO[2]))
				pathToFile=paste(fileLocation,filename,sep="/")
				print("")
				print("segment image")
				# segmentation of the image
				imageData=try(segmentImage(filename=pathToFile,,maxShape=800,minShape=40,failureRegion=2000))
				print(inherits(imageData,"try-error"))
				#original image
				img=imageData[[1]]
				#segmented image
				imgW=imageData[[2]]
				#features of the nulcei
				cellData=imageData[[3]]
				
				
				#indices of all white pixel in the image
				indexWhitePixel=imageData[[4]]
				#classify the nuclei, or is the image white?
				classify=imageData[[5]]
				#is the image a failure image(no classification)
				nameFile=strsplit(filename,"\\.")[[1]][1]
				#section of the image
				slice=blockSlice[as.character(blockSlice$block)==nameFile,2]
				
				if(classify==TRUE){
				print("classify cells and paint image")
			    #classify the cells
				classValues=classifyCells(classifier,image=img,segmentedImage=imgW,featuresObjects=cellData,paint=TRUE,KS=TRUE,cancerIdentifier=cancerIdentifier)
				classes=classValues[[1]]
				writeImage(classValues[[2]],paste(pathToOutputFolderImgDir,filename,sep="/"))
				colnames(classValues[[3]])=c("index","g.x","g.y","class")
				write.table(classValues[[3]],paste(pathToOutputFolderImgDirCellCoordinates,paste(nameFile,".txt",sep=""),sep="/"),row.names=FALSE,quote=FALSE)
				classes=as.character(classes)
				classifiedCells=cellData
				
				print("cells sub window")
				#calculation of cellularity, drawing of the heatmap
				cellValues=determineCellularity(classes,classifiedCells,dim(imgW),img,imgW,indexWhitePixel,"c",sort(unique(classes)))
				print(inherits(cellValues,"try-error"))
				cellsSubImages=cellValues[[1]]
				cellsDensityImage=cellValues[[2]]
				print("write cellsSubImages")
				print(paste(pathToOutputFolderImgDirFiles,paste(paste("section",slice,sep="_"),nameFile,sep="/"),sep="/"))
				cellValuesResult=data.frame(cellsSubImages[1],cellsSubImages[2],cellsSubImages[3],cellsSubImages[4])
				colnames(cellValuesResult)=names(cellsSubImages)
				write.table(cellValuesResult,paste(pathToOutputFolderImgDirFiles,paste(paste("section",slice,sep="_"),nameFile,sep="/"),sep="/"),sep="\t",row.names = FALSE)
				
				print("img density")
				#draw the heatmap
				writeImage(cellsDensityImage,paste(pathToOutputFolderImgDirCellDensity,paste(nameFile,".jpg",sep=""),sep="/"))
				print("img neighbors")
				print("image saved")
				classesN=classes
				classesN[which(classes==cancerIdentifier)]=1
				classesN[which(classes!=cancerIdentifier)]=2
				classifiedCells[,"class"]=as.numeric(classesN)
				print("feature file written")
				classifiedCells
				}else{
					nameFile=strsplit(filename,"\\.")[[1]][1]
					print("no classification")
					print(nameFile)
					if(length(dim(img))>2){
					cellsDensityImage=img
					}else{
					cellsDensityImage=img
					}
					cellsDensityImage[,,1]=1
					cellsDensityImage[,,2]=1
					cellsDensityImage[,,3]=1
					writeImage(cellsDensityImage,paste(pathToOutputFolderImgDirCellDensity,paste(nameFile,".jpg",sep=""),sep="/"))

					print("sliceColors")
					print(sliceColors)
					print(paste("slice",slice))
					if(slice>dim(sliceColors)[2]){
						img[,,1]=1
						img[,,2]=1
						img[,,3]=1
					}else{
						print("new color")
						img[,,1]=sliceColors[1,as.numeric(slice)]/255
						img[,,2]=sliceColors[2,as.numeric(slice)]/255
						img[,,3]=sliceColors[3,as.numeric(slice)]/255
					}
					
				print(paste("values for small Image",index))	
				writeImage(img,paste(pathToOutputFolderImgDir,paste(nameFile,".jpg",sep=""),sep="/"))
					
				classifiedCells=data.frame()
				}
}
processAperio=function(classifier=classifier,inputFolder=inputFolder,outputFolder=outputFolder,identifier=identifier,numSlides=numSlides,cancerIdentifier=cancerIdentifier){
	options(stringsAsFactors = FALSE)
	pathToFolder=inputFolder
	pathToOutputFolder=outputFolder
	pathToOutputFolderTempFiles=paste(pathToOutputFolder,"tempFiles",sep="/")
	dir.create(pathToOutputFolderTempFiles)
	
	dir.create(pathToOutputFolder)
	pathToOutputFolderImgDir=paste(pathToOutputFolder,"classifiedImage",sep="")
	pathToOutputFolderImgDirFiles=paste(pathToOutputFolder,"Files",sep="")
	pathToOutputFolderImgDirCellDensity=paste(pathToOutputFolder,"tumourDensity",sep="")
	pathToOutputFolderImgDirCellCoordinates=paste(pathToOutputFolder,"cellCoordinates",sep="")
	dir.create(pathToOutputFolderImgDir)
	dir.create(pathToOutputFolderImgDirFiles)
	dir.create(pathToOutputFolderImgDirCellDensity)
	dir.create(pathToOutputFolderImgDirCellCoordinates)
	#finds the right section for every subimage
	sliceSizeList=try(findSlices(inputFolder,pathToOutputFolder,numSlides))
	print(inherits(sliceSizeList,"try-error"))
	print(sliceSizeList)
	blockSlice=sliceSizeList[[1]]
	sizeO=sliceSizeList[[2]]
	print(paste("sizeO",sizeO))
	print(inherits(blockSlice,"try-error"))
	numberSlices=length(unique(blockSlice[,2]))
	print(paste("numberSlices",numberSlices))
	sliceFolder=c()
	for (i in 1:numberSlices){
		dir.create(paste(pathToOutputFolderImgDirFiles,paste("section",i,sep="_"),sep="/"))
	}
	sliceColors=col2rgb(c("red","blue","green","yellow","orange"))
	dir.create(pathToOutputFolderImgDirCellDensity)
	filenames=list.files(path =pathToFolder ,pattern=identifier)
	print(filenames)
	allCells=foreach( i = 1:length(filenames), .combine=rbind) %do% {
			classificationError=try(classificationAperio(pathToFolder,filenames[i],pathToOutputFolderImgDir,classifier,pathToOutputFolderImgDirFiles,pathToOutputFolderImgDir,pathToOutputFolderImgDirCellDensity,blockSlice,sliceColors,sizeO,i,cancerIdentifier,pathToOutputFolderImgDirCellCoordinates))
			print(inherits(classificationError,"try-error"))
		}
		print("write all cells")
		write.table(allCells,paste(pathToOutputFolderImgDir,paste("result",".txt",sep=""),sep="/"),col.names=TRUE,sep="\t",row.names=FALSE)
}