sampleRows <- function(dataFrame, num){
	dataFrame <- dataFrame[sample(nrow(dataFrame), num), ]
	return(dataFrame)
}

sampleColumns <- function(dataFrame, num){
	dataFrame <- dataFrame[,sample(nrow(dataFrame), num)]
	return(dataFrame)
}
