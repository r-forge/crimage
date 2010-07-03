source("/Users/failmezger/cambridge/tumor_IP/RSource/CRImageN2.r")
img=img=readImage("Ss46.jpg")
trainingData=read.table("/Users/failmezger/CRImage/data/trainingData.txt",header=T)
trainingDataClasses = read.delim("/Users/failmezger/cambridge/tumor_IP/workspace/trainingData/5207_Ss30_mC_myClassification_n2.txt",header=T)
trainingData[,"class"]=trainingDataClasses[,"class"]
trainingData[1458,"class"]=NA
trainingData=subset(trainingData, select = -c(hr.asm,hr.con,hr.cor,hr.var,hr.idm,hr.sav,hr.sva,hr.sen,hr.ent,hr.dva,hr.den,hr.f12,hr.f13,hg.asm,hg.con,hg.cor,hg.var,hg.idm,hg.sav,hg.sva,hg.sen,hg.ent,hg.dva,hg.den,hg.f12,hg.f13,hb.asm,hb.con,hb.cor,hb.var,hb.idm,hb.sav,hb.sva,hb.sen,hb.ent,hb.dva,hb.den,hb.f12,hb.f13))


write.table(trainingData,file="/Users/failmezger/CRImage/data/trainingData.txt",sep="\t",rownames=F)


f = system.file("data", "trainingData.txt", package="CRImage")
#read training data
trainingData=read.table(f,header=TRUE)
#create classifier
classifier=createClassifier(trainingData,topo=F)[[1]]

t = system.file("data", "trainingData.txt", package="CRImage")
#read training data
trainingData=read.table(t,header=TRUE)
#create classifier
classifier=createClassifier(trainingData,topo=FALSE)[[1]]
#classify cells
f = system.file("data", "exImg2.jpg", package="CRImage")
classes=classifyCells(classifier,filename=f,paint=TRUE,KS=FALSE)



cellularityValues=calculateCellularity("exImg.jpg",classifier=classifier,cancerIdentifier="1")

filename="Ss46.jpg"
f="/Users/failmezger/exImg.jpg"
exImg=readImage("exImg.jpg")

f = system.file("data", "exImg.png", package="CRImage")
exImgB=readImage(f)
#find white pixels and exclude them from thresholding(if white is background)
indexWhitePixel=which(exImgB[,,1]>0.85 &exImgB[,,2]>0.85 & exImgB[,,3]>0.85)
#convert to grayscale
exImgB=channel(exImgB,"gray")
#calculate threshold
t=calculateThreshold(as.vector(exImgB[-indexWhitePixel]))
#create binary image
exImgB[exImgB>t]=-1
exImgB[exImgB != -1]=0
exImgB[exImgB == -1]=1
#display
display(exImgB)


#segment image
f = system.file("data", "exImg.png", package="CRImage")
segmentationValues=segmentImage(f)
image=segmentationValues[[1]]
segmentedImage=segmentationValues[[2]]
imageFeatures=segmentationValues[[3]]

#classify cells
f = system.file("data", "exImg.png", package="CRImage")
classes=classifyCells(classifier,filename=f,KS=T)


f = system.file("data", "exImg.png", package="CRImage")
trainingValues=createTrainingSet(filename=f)


f = system.file("data", "trainingData.txt", package="CRImage")
#read training data
trainingData=read.table(f)
#create classifier
classifier=createClassifier(trainingData)[[1]]

#calculation of cellularity
f = system.file("data", "exImg2.png", package="CRImage")
exImg=readImage(f)
cellularityValues=calculateCellularity(filename=f,classifier=classifier,cancerIdentifier="1")




######
library(foreach)
library(EBImage)
source("/Users/failmezger/cambridge/tumor_IP/RSource/CRImageN2.r")
pathToImage="/Users/failmezger/package/5207_mod"
pathToOutput="/Users/failmezger/package/output"
load("/Users/failmezger/package/classifier")
classifier=model
processAperio(classifier,pathToImage,pathToOutput,"Da",2,"c")


processAperio(classifier,pathToImage,pathToOutput,"Da",3,"c")

package.skeleton(name="CRImage",list=c("calculateThreshold","segmentImage","createTrainingSet","createClassifier","classifyCells","paintCells","calculateCellularity","determineCellularity","processAperio","findSlices","classificationAperio","kernelSmoother"))