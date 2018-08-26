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
  require(MKKC)
  require(rlist)
  
  #Process the data set
  data.views <- preProcessData(data.views)
  
  
  #Calculate kernel/similarity matrices  --kernel parameter optimization --- 
  #standardize matrices for comparison --- transform into a matrix array
  kernels <- lapply(c(1:length(data.views)), FUN = function(idx){chooseKernel(data.views[[idx]])})
  kernels <- lapply(c(1:length(kernels)), FUN = function(idx){StandardizeKernel(kernels[[idx]])})
  kernel.array <- do.call(abind, c(kernels, list(along = 3)))
  
  mkkm.model <- parallel::mclapply(c(2:10), FUN = function(k){runMKKM(kernel.array, clusters = k)})
  
  combined.kernels <- lapply(mkkm.model, FUN = function(model){model$combined.kernel})
  model.clustering <- lapply(mkkm.model, FUN = function(model){model$state$clustering})
  model.kernel.coefficients <- lapply(mkkm.model, FUN = function(model){model$state$theta})
  
  mkkc.nmi <- lapply(model.clustering, FUN = function(clustering){NormMI(truelabel, model.clustering)})

  #################################################################################################################
  
  ## next, construct similarity graphs
  W1 = affinityMatrix(Dist1, K, alpha)
  W2 = affinityMatrix(Dist2, K, alpha)
  
  #Run SNF
  W = SNF(Wall = list(W1,W2), K = K, t = iter)
  
  
  #Cluster
  group = spectralClustering(W, C);
  group.1 = spectralClustering(W1, C);
  group.2 = spectralClustering(W2, C);
  
  #"Clustering results of individual and fused data views"
  SNFNMI = calNMI(group, truelabel)
  SNFNMI.1 = calNMI(group.1, truelabel)
  SNFNMI.2 = calNMI(group.2, truelabel)
  
  return(list(affinity.matrices = list(W1 = W1, W2 = W2, W = W),
              nmi.values = list(nmi.1 = SNFNMI.1, nmi.2 = SNFNMI.2, nmi.fused = SNFNMI),
              identity = list(ident.1 = group.1, ident.2 = group.2, ident.fused = group)))
}

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

combineKernels <- function(kernels, kernel.params){
  for(i in 1:length(kernel.params)){
    #print(i)
    kernels[, ,i] <- kernels[, ,i] * kernel.params[i]
  }
  return(apply(kernels, c(1,2), sum))
}

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


chooseKernel <- function(data.view, type = "nn", k = 30){
  require(vegan)
  
  if(type == "nn"){
    #Jacard similarity Kernel based on nearest neighbors
    neighbors <- getNearestNeighbors(data.view, k)
    jaccard.kernel <- jaccardKernel(neighbors)
    return(jaccard.kernel)
  }
}

################################
optimizeKernelParam <- function(data.view, params){
  
}


rbfKernel <- function(data.view){
  require(kernlab)
  
  param.choices.log <- c(1 * 10^(-2:3))
  kernels.log <- lapply(param.choices.log, FUN = rbfdot(sigma = )
  
  
}
################################
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


getNearestNeighbors <- function(data.view, k){
  require(FNN)
  
  neighbors <- get.knn(t(data.view), k = k)
  neighbors <- neighbors$nn.index
  rownames(neighbors) <- colnames(data.view)
  return(neighbors)
}

normalizeKernel <- function(kernel){
  #Normalizing the kernel as per chapter 5 in kernel methods for pattern analysis
  D <- diag(1./sqrt(diag(kernel)))
  kernel <- D * kernel * D
}
