runMKKMPipeline <- function(data.views, truelabel,
                            K = 20, alpha = 0.5, iter = 10,
                            down.pct = 1, C = 2,
                            simulation = FALSE,
                            kernel.normalization = TRUE,
                            choose.variable.genes = FALSE,
                            prcomp = TRUE,
                            kernel = "nn"){

  #Requirements  
  require(abind)
  require(SNFtool)#for preprocessing of the data
  require(MKKC)
  require(rlist)
  require(clValid)
  require(psych)
  
  #Process the data set -- in progress
  #data.views <- preProcessData(data.views)
  Data1
  
  
  #Calculate kernel/similarity matrices  --kernel parameter optimization --- 
  #standardize matrices for comparison --- transform into a matrix array
  #kernels <- lapply(c(1:length(data.views)), FUN = function(idx){chooseKernel(data.views[[idx]])}) -- Jaccard Kernel function
  kernels <- lapply(c(1:length(data.views)), FUN = function(idx){optimizeRBFKernel(data.views[[idx]])})#Optimized Kernels based on KSdist
  #kernels <- lapply(c(1:length(kernels)), FUN = function(idx){StandardizeKernel(kernels[[idx]])}) -- don't do normalization???
  kernel.array <- do.call(abind, c(kernels, list(along = 3)))
  
  mkkm.model <- parallel::mclapply(c(2:50), FUN = function(k){runMKKM(kernel.array, clusters = k)})
  
  combined.kernels <- lapply(mkkm.model, FUN = function(model){model$combined.kernel})
  model.clustering <- lapply(mkkm.model, FUN = function(model){model$state$clustering})
  model.kernel.coefficients <- lapply(mkkm.model, FUN = function(model){model$state$theta})
  
  #Normalized mutual information for when we know the labels
  mkkc.nmi <- lapply(model.clustering, FUN = function(clustering){NormMI(truelabel, model.clustering)})
  
  #Internal measures for when we do not know the true labels...
  #calculate distances
  
  data.views.dr <- lapply(data.views, FUN = function(data.view){prcomp(t(data.view), rank. = 100)$x})
  
  if(dist.type = "cor"){
    
    dists <- lapply(data.views, function(data.view){as.dist((1 - cor(data.view))/2)})
    
  }else{
    dists <- lapply(data.views.dr, FUN = function(data.view){dist(data.view, method = dist.type)})
  }
  
  quality.metrics <- lapply(dists, calcInternalClusterQualityMeasures, mkkm.model)
  
  plotInternalQulaitymetrics()

  #Feature/cluster id mutual information based cluster number optimization

  #PC based
  average.cluster.mi <- lapply(data.views.dr, FUN = function(data.view){chooseClustersMI(data.view, model.clustering)})

  #Gene based ------
  features.use <- lapply(data.views, FUN = function(data.view){chooseVariableFeatures(data.view)})
  average.cluster.mi <- lapply(data.views, FUN = function(data.view){chooseClustersMI(data.view, model.clustering, features.use)})
  
  
  
  
  #Choose the optimal clustering parameter setting based on dunn index and return relevant 
  
  which.min(quality.metrics[])

  return(list(affinity.matrices = list(W1 = W1, W2 = W2, W = W),
              nmi.values = list(nmi.1 = SNFNMI.1, nmi.2 = SNFNMI.2, nmi.fused = SNFNMI),
              identity = list(ident.1 = group.1, ident.2 = group.2, ident.fused = group)))

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

##############################################################################################################
preProcessData <- function(data.views, simulation, choose.variable.genes, prcomp){
  #Process the data set
  if(simulation){
    data.views <- lapply(c(1:length(data.views)), FUN = function(idx){as.matrix(data.views[[idx]][, c("V1", "V2")])})
    data.views <- lapply(c(length(data.views)), FUN = function(idx){t(data.views[[idx]])})
  }else{
    data.views <- lapply(c(1:length(data.views)), FUN = function(idx){as.matrix(data.views[[idx]])})
  }
  
  #Normalization section
  #Sample normalization to control for absolute size differences (ie num transcripts detected) between samples -- only for positive count data
  if(sample.normalize){
    data.views <- lapply(c(1:length(data.views)), FUN = function(idx){data.views[[idx]] / colSums(data.views[[idx]])})
  }
  
  
  #Choose variable genes -- TODO
  if(choose.variable.genes){
    data.views <- lapply(data.views, FUN = 
                         )
  
  
  #Principal components - DR
  
  if(prcomp){TODO principal components analysis}
  
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
  model$combined.kernel <- combineKernels(kernel.array, state$theta)
  
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
optimizeRBFKernel <- function(data.view){
  require(kernlab)
  require(parallel)
  require(FNN)
  
  param.choices.log <- c(1 * 10^(-8:3))
  kernels.log <- lapply(param.choices.log, FUN = rbfdot)
  
  log.features <- mclapply(kernels.log,
                           FUN = function(kernel){computeKernelMatrix(data.view, kernelFunction = kernel)},
                           mc.cores = (detectCores()-1))
  
  #Choose best kernel parameter based on KL divergence...resultoing distribution of values in kernel should be close to beta distribution
  #Centered around 0.5, with values between 0 and 1 (symmetric)
  
 test.results <- lapply(log.features, FUN = computeKSDist)
 min.idx <- which.min(test.results)
 
 #Test different parameter values on a linear scale based on the "optimal" value from above
 param.choices.linear <- unlist(lapply(param.choices.log[c((min.idx -1) : (min.idx + 1))],
                                         FUN = function(scale.value){scale.value * seq(1,10, 2)}))
 kernels.linear <- lapply(param.choices.linear, FUN = rbfdot)
 linear.features <- mclapply(kernels.linear,
                             FUN = function(kernel){computeKernelMatrix(data.view, kernelFunction = kernel)},
                             mc.cores = (detectCores()-1))
 #Calculate distances again
 test.results <- lapply(linear.features, FUN = computeKSDist)
 min.idx <- which.min(test.results)
 
 return(linear.features[[min.idx]])
 
}

##############################################################################################################
computeKSDist <- function(kern){
  test.result <- suppressWarnings(ks.test(x = c(kern), y = rbeta(length(kern), shape1 = 2, shape2 = 2)))
  return(test.result[["statistic"]])
}
##############################################################################################################
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