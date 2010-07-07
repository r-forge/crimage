`classificationAperio` <-
function(fileLocation,filename,pathToOutputFolderImgDir,classifier,pathToOutputFolderImgDirFiles,pathToOutputFolderImgDirImages,pathToOutputFolderImgDirCellDensity,blockSlice,sliceColors,sizeO,index,cancerIdentifier,pathToOutputFolderImgDirCellCoordinates){
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

