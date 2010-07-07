`paintCells` <-
function(imgWT,img,classes,index,classValues){
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

