`kernelSmoother` <-
function(predictedClasses,cellCoordinates,segmentedImage,cancerIdentifier,indexCells){
				#options(stringsAsFactors = FALSE)
				prob=attr(predictedClasses, "prob")
				print(prob)
				probCancer=prob[,cancerIdentifier]
				print("prob cancer")
				print(probCancer)
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
							print(index)
							classValues=f
							newClassesKS=data.frame(index,classValues)
							print(allNewClasses)
							print(newClassesKS)
							if(windowCounter==1){
								allNewClasses=newClassesKS
							}else{
								print(names(allNewClasses))
								print(names(newClassesKS))
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
			allNewClasses$classValues[allNewClasses$classValues==1]="o"
			allNewClasses$classValues[allNewClasses$classValues==2]=cancerIdentifier
			
			realIndex=data.frame(indexCells)
			colnames(realIndex)="index"
			colnames(allNewClasses)=c("index","classValues")
			classesRealIndices=merge(realIndex,allNewClasses,all.x=TRUE,all.y=TRUE)
			allClasses=classesRealIndices[,"classValues"]
			print(allClasses)
			classSVM=allClasses
			classSVM
}

