`processAperio` <-
function(classifier=classifier,inputFolder=inputFolder,outputFolder=outputFolder,identifier=identifier,numSlides=numSlides,cancerIdentifier=cancerIdentifier){
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

