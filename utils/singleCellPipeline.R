
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

runSIMLRPipeline <- function(object, pca = TRUE, num.pcs = NULL, smoothing.fxn.location = "~/Documents/uva/master/the/utils/knnSmoothing.R", n.sample){
  require(SIMLR)
  source(smoothing.fxn.location)
  
  smooth.exp <- knn_smoothing(mat = object@data[1:100, 1:100], k = 4, d = 10)
  
  object <- AddSmoothedScore(object)
  
  
 t(object@dr$pca@cell.embeddings)
  
  
 #This worked out, however, I'd like to try and avoid 
  idx <- sample(x = 1:ncol(data.matrix),size = n.sample)
  simlr.result <- SIMLR(X = t(object@dr$pca@cell.embeddings[1:1000,1:num.pcs]), c = 6, k = 25, cores.ratio = 0.5, normalize = FALSE)
  
}

calcNearestNeighbors <- function(data.matrix){
  require(RANN)
  require(dbscan)
  #data.matrix is a samples by features matrix (n x m)
  
  nn.result <- nn2(data = data.matrix, k = 10)
  
  nn.result$nn.dists
  
  kmeans.results <- kmeans(x = data.matrix, centers = 1000, iter.max = 1000)
  
  
  dbscan.result <- dbscan(x = data.matrix, eps = 1, borderPoints = FALSE, minPts = 3)
  
}


SNFSNNgraph <- function(){
  require(SNFtool)
  
  d1 <- dist2(dipg@dr$pca@cell.embeddings, dipg@dr$pca@cell.embeddings)
  W1 <- affinityMatrix(d1)
  
  resolutions <- seq(0.1, 2.8, 0.3)
  k.params <- seq(10, 100, 10)
  num.pcs <- (1:25)
  graphs <- lapply(num.pcs, FUN = function(pc){getSNNgraph(object, dims.use = 1:pc)})
  
  graph.1 <- as.matrix(BuildSNN(object = object, dims.use = 1:5)@snn)
  graph.2 <- as.matrix(BuildSNN(object = object, dims.use = 6:10)@snn)
  
  idx <- sort(sample(1:2458,500))
  
  graphs <- list(graph.1[idx, idx], graph.2[idx, idx])
  
  normalize <- function(X) {
    return(X/rowSums(X))
  }
  
  avg.graph <- (graph.1 + graph.2)/2
  W <- normalize(avg.graph)
  W <- W + t(W)
  W <- W * 1000
  temp.res <- SNF(graphs, K = 50, t = 10)
  
  temp.graphs <- lapply(graphs, FUN = function(graph){graph[1:1000, 1:1000]})
  
}

getSNNgraph <- function(object, k.param = 30, dims.use){
  require(Seurat)
  object <- BuildSNN(object, k.param = k.param, save.SNN = TRUE, force.recalc = TRUE, dims.use = dims.use)
  snn.graph <- as.matrix(object@snn)
  return(snn.graph)
}




