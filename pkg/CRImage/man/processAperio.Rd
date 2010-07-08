\name{processAperio}
\alias{processAperio}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
processAperio(classifier = classifier, inputFolder = inputFolder, outputFolder = outputFolder, identifier = identifier, numSlides = numSlides, cancerIdentifier = cancerIdentifier)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{classifier}{ The classifier.}
  \item{inputFolder}{ The path to the image folder. }
  \item{outputFolder}{ The path to the output folder. }
  \item{identifier}{ The identifier of the files (Ss or Da) }
  \item{numSlides}{ The number of sections in the image. }
  \item{cancerIdentifier}{ The identifier of the cancer class}
}
\details{
  ~~ If necessary, more details than the description above ~~
}
\value{
  ~Describe the value returned
  If it is a LIST, use
  \item{comp1 }{Description of 'comp1'}
  \item{comp2 }{Description of 'comp2'}
  ...
}
\references{ ~put references to the literature/web site here ~ }
\author{ Henrik Failmezger }
\note{ ~~further notes~~ 

 ~Make other sections like Warning with \section{Warning }{....} ~
}
\seealso{ ~~objects to See Also as \code{\link{help}}, ~~~ }
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
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
	dir.create(pathToOutputFolderImgDir)
	dir.create(pathToOutputFolderImgDirFiles)
	dir.create(pathToOutputFolderImgDirCellDensity)
	#finds the right section for every subimage
	sliceSizeList=try(findSlices(inputFolder,pathToOutputFolder,numSlides,pathToScripts))
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
	allCells=foreach( i = 1:length(filenames), .combine=rbind) \%do\% {
			classificationError=try(classificationAperio(pathToFolder,filenames[i],pathToOutputFolderImgDir,classifier,pathToOutputFolderImgDirFiles,pathToOutputFolderImgDir,pathToOutputFolderImgDirCellDensity,blockSlice,sliceColors,sizeO,i,cancerIdentifier))
			print(inherits(classificationError,"try-error"))
		}
		print("write all cells")
		write.table(allCells,paste(pathToOutputFolderImgDir,paste("result",".txt",sep=""),sep="/"),col.names=T,sep="\t",row.names=F)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{misc}
