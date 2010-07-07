`segmentImage` <-
function(filename="",image=NA,maxShape=NA,minShape=NA,failureRegion=NA){
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

