`createClassifier` <-
function(trainingData,cross=FALSE,topo=TRUE){
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

