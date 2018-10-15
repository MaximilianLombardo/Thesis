
readRegevData <- function(file.path = "~/Documents/uva/master/data/RegevSingleCell/GSE102130_K27Mproject.RSEM.vh20170621.txt"){
  require(Matrix)
  if(grepl(pattern = "RSEM", x = file.path)){
    sc.data <- read.table(file.path, header = TRUE, row.names = 1)
  }else{
    sc.data <- read.table(file.path)
  }
  
  #convert from data frame to numeric matrix
  sc.data <- as.matrix(sc.data)
  #Convert single data matrix to sparse matrix
  sc.data <- Matrix(sc.data, sparse = TRUE)
  return(sc.data)
}

runSIMLRPipeline <- function(object, pca = TRUE, num.pcs = NULL){
  require(SIMLR)
  
 t(object@dr$pca@cell.embeddings)
  
  
  simlr.result <- SIMLR(X = t(object@dr$pca@cell.embeddings), c = 5, k = 25, cores.ratio = 0.5, normalize = FALSE)
  
}