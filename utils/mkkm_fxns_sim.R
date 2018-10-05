#The Multiple kernel k means pipeline for simulated data
runMKKMPipeline <- function(Data1, Data2, truelabel,
                            C = 2,
                            scale.factor = 1/20,
                            kernel = "nn"){
  
  #Requirements  
  require(SNFtool)#for preprocessing of the data
  require(MKKC)

  Data1 <- Data1[, c("V1", "V2")]
  Data2 <- Data2[, c("V1", "V2")]
  data.views <- list(Data1, Data2)
  
  #Process the data set
  data.views <- lapply(data.views, FUN = function(view){as.matrix(view)})
  
  #Normalization
  data.views <- lapply(data.views, FUN = standardNormalization)
  
  #Create kernels for each data view
  data.features <- makeKernels1(data.views = data.views, scale.factor)#Make kernels
  data.features <- lapply(data.features, FUN = function(data.feature){StandardizeKernel(x = data.feature)}) #Standardize
  
  K1 <- data.features[[1]]
  K2 <- data.features[[2]]
  
  n.view = length(data.features)    # the number of views used
  
  K.array = array(as.numeric(unlist(data.features)), dim = c(nrow(data.features[[1]]), ncol(data.features[[1]]), n.view))
  
  #res <- mkkc(K = K, centers = 2)
  model <- runMKKM(kernel.array = K.array)
  K <- model$combined.kernel
  
  #Cluster
  group = model$state$clustering
  group.1 = runkkmeans(K1, C)
  group.2 = runkkmeans(K2, C)
  
  #"Clustering results of individual and fused data views"
  mkkmNMI = calNMI(group, truelabel)
  mkkmNMI.1 = calNMI(group.1$clustering, truelabel)
  mkkmNMI.2 = calNMI(group.2$clustering, truelabel)
  

  return(list(kernel.matrices = list(K1 = K1, K2 = K2, K = model$combined.kernel),
              nmi.values = list(nmi.1 = mkkmNMI.1, nmi.2 = mkkmNMI.2, nmi.fused = mkkmNMI),
              identity = list(ident.1 = group.1$clustering, ident.2 = group.2$clustering, ident.fused = group)))

}
##############################################################################################################33333
chooseClustersMI <- function(data.view, clusterings, num.pcs = 1:20, features.use = TRUE){
  require(parallel)

  if(!is.null(features.use)){
    features <- split(data.view[features.use,], row(data.view))
  }else{
    features <- split(data.view[,num.pcs], col(data.view))
  }
  
  
  mean.mut.info <- mclapply(clusterings, FUN = function(clustering){calcFeatureMI(features, clustering)}, mc.cores = (detectCores() - 1))
  
  return(mean.mut.info)
}
##############################################################################################################
calcFeatureMI <- function(features, clustering, num.pcs = NULL, genes = TRUE){
  require(mpmi)
  if(genes){
    mean.mut.info <- mean(unlist(lapply(features, FUN = function(feature){mmi.pw(feature, clustering)$mi})))
    
  }else{
    mean.mut.info <- mean(unlist(lapply(features[num.pcs], FUN = function(feature){mmi.pw(feature, clustering)$mi})))
    
  }
  
  
  return(mean.mut.info)
}
##############################################################################################################
chooseVariableFeatures <- function(data.view, feat.var.thresh = 6){
  feature.idx <- apply(data.view, 1, FUN = function(feature){any(feature > feat.var.thresh)})
  
  variable.features <- names(feature.idx[feature.idx])
  
  return(variable.features)
}
##############################################################################################################
calcInternalClusterQualityMeasures <- function(experimental.dist, models){
  require(clValid)
  quality.metrics <- list()
  clusterings <- lapply(c(1:length(models)), FUN = function(idx){models[[idx]]$state$clustering})
  quality.metrics$connectivity.scores <- unlist(lapply(clusterings, FUN = function(clustering){connectivity(distance = experimental.dist, clusters = clustering)}))
  quality.metrics$dunn <- unlist(lapply(clusterings, FUN = function(clustering){dunn(distance = experimental.dist, clusters = clustering)}))
  quality.metrics$silhouette <- unlist(lapply(clusterings, FUN = function(clustering){mean(silhouette(clustering, dist = experimental.dist)[,"sil_width"])}))
  
  return(quality.metrics)
}


  ##############################################################################################################
runMKKM <- function(kernel.array, clusters = 2, iter = 10){
  #source(fxn.location)
  require(MKKC)
  
  model <- list()
  
  #model$res <- mkkc(kernel.array, centers = clusters, iter.max = iter)
  ###########################################################################################
  #mehmet
  parameters <- list()
  parameters$cluster_count <- clusters
  parameters$iteration_count <- iter
  model$state <- mkkmeans_train(kernel.array, parameters)
  ###########################################################################################
  model$combined.kernel <- combineKernels(kernel.array, model$state$theta)
  
  return(model)
  
}
##############################################################################################################
combineKernels <- function(kernels, kernel.params){
  for(i in 1:length(kernel.params)){
    #print(i)
    kernels[, ,i] <- kernels[, ,i] * kernel.params[i]
  }
  return(apply(kernels, c(1,2), sum))
}
##############################################################################################################
mkkmeans_train <- function(Km, parameters) {
  require(Rmosek)
  state <- list()
  state$time <- system.time({
    P <- dim(Km)[3]
    theta <- rep(1 / P, P)
    K_theta <- matrix(0, nrow(Km), ncol(Km))
    for (m in 1:P) {
      K_theta <- K_theta + theta[m]^2 * Km[,,m]  
    }
    
    objective <- rep(0, parameters$iteration_count)
    for (iter in 1:parameters$iteration_count) {
      #print(sprintf("running iteration %d...", iter))
      H <- eigen(K_theta, symmetric = TRUE)$vectors[, 1:parameters$cluster_count]
      
      problem <- list()
      problem$sense <- "min"
      problem$c <- rep(0, P)
      problem$A <- Matrix(1, nrow = 1, ncol = P, sparse = TRUE)
      problem$bc <- rbind(blc = 1, buc = 1) 
      problem$bx <- rbind(blx = rep(0, P), bux = rep(1, P))
      problem$qobj <- list(i = 1:P, j = 1:P, v = sapply(1:P, function(m) {sum(diag(Km[,,m])) - sum(diag(t(H) %*% Km[,,m] %*% H))}))
      opts <- list()
      opts$verbose <- 0
      result <- mosek(problem, opts)
      theta <- result$sol$itr$xx
      K_theta <- matrix(0, nrow(Km), ncol(Km))
      for (m in 1:P) {
        K_theta <- K_theta + theta[m]^2 * Km[,,m]  
      }
      
      objective[iter] <- sum(diag(t(H) %*% K_theta %*% H)) - sum(diag(K_theta))
    }
    H_normalized <- H / matrix(sqrt(rowSums(H^2, 2)), nrow(H), parameters$cluster_count, byrow = FALSE)
    
    set.seed(NULL)
    state$clustering <- kmeans(H_normalized, centers = parameters$cluster_count, iter.max = 1000, nstart = 10)$cluster
    state$objective <- objective
    state$parameters <- parameters
    state$theta <- theta
  })
  return(state)
}

runkkmeans <- function(K, C){
  #initalize the parameters of the algorithm
  parameters <- list()
  
  #set the number of clusters
  parameters$cluster_count <- C
  
  #perform training
  state <- kkmeans_train(K, parameters)
  
  return(state)
}

kkmeans_train <- function(K, parameters) {
  state <- list()
  state$time <- system.time({
    H <- eigen(K, symmetric = TRUE)$vectors[, 1:parameters$cluster_count]
    objective <- sum(diag(t(H) %*% K %*% H)) - sum(diag(K))
    H_normalized <- H / matrix(sqrt(rowSums(H^2, 2)), nrow(H), parameters$cluster_count, byrow = FALSE)
    
    set.seed(NULL)
    state$clustering <- kmeans(H_normalized, centers = parameters$cluster_count, iter.max = 1000, nstart = 10)$cluster
    state$objective <- objective
    state$parameters <- parameters
  })    
  return(state)
}

##############################################################################################################
chooseKernel <- function(data.view, type = "nn", k = 30){
  require(vegan)
  
  if(type == "nn"){
    #Jacard similarity Kernel based on nearest neighbors
    neighbors <- getNearestNeighbors(data.view, k)
    jaccard.kernel <- jaccardKernel(neighbors)
    return(jaccard.kernel)
  }
}




##############################################################################################################
makeKernels1 <- function(data.views, scale.factor = 1){
  require(kernlab)
  require(parallel)
  kernel.functions <- lapply(data.views, FUN = function(data.view){rbfdot(sigma = (scale.factor * sqrt(nrow(data.view))))})
  data.features <- lapply(c(1:length(kernel.functions)),
                          FUN = function(idx){kernelMatrix(x = data.views[[idx]],
                                                           kernel = kernel.functions[[idx]])})
  return(data.features)
}

computeKernelMatrix <- function(data, kernelFunction){
  
  kernel <- matrix( ,nrow = ncol(data), ncol = ncol(data))
  
  row.names(kernel) <- colnames(data)
  colnames(kernel) <- colnames(data)
  
  for(i in c(1:ncol(data))){
    for(j in c(1:ncol(data))){
      
      kernel[i,j] <- kernelFunction(data[,i], data[,j])
      
    }
  }
  return(kernel)
}
##############################################################################################################
trace <- function(data)sum(diag(data))
##############################################################################################################
jaccardKernel <- function(neighbors){
  require(vegan)
  require(parallel)
  
  mat <- matrix(numeric(0), nrow(neighbors), nrow(neighbors))
  
  colnames(mat) <- rownames(neighbors)
  rownames(mat) <- rownames(neighbors)
  
  for(i in 1:nrow(neighbors)){
    for(j in 1:nrow(neighbors)){
      
      #Calc Jaccard distance/similarity
      mat[i,j] <- length(intersect(neighbors[i,], neighbors[j,]))/length(union(neighbors[i,], neighbors[j,]))
    }
  }
  
  return(mat)
}
##############################################################################################################
getNearestNeighbors <- function(data.view, k){
  require(FNN)
  
  neighbors <- get.knn(t(data.view), k = k)
  neighbors <- neighbors$nn.index
  rownames(neighbors) <- colnames(data.view)
  return(neighbors)
}
##############################################################################################################
normalizeKernel <- function(kernel){
  #Normalizing the kernel as per chapter 5 in kernel methods for pattern analysis
  D <- diag(1./sqrt(diag(kernel)))
  kernel <- D * kernel * D
  return(kernel)
}