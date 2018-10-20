runDreamPipeline <- function(object, bmtmkl.dir = "~/Documents/gitRepos/master/utils/", k = 5, C = 5){
  #Sourcing the functions needed for bayesian multitask multiple kernel learning
  source(paste0(bmtmkl.dir, "bayesian_multitask_multiple_kernel_learning_train.R"))
  source(paste0(bmtmkl.dir, "bayesian_multitask_multiple_kernel_learning_test.R"))
  require(abind)
  require(mlr)
  #####
  #Split up the task matrix (the target matrix) and the input data matrices
  drug.response <- object$Drug_Response_Training
  object <- object[-1]
  
  #Get the cell lines which are common to all experiments and for which we have drug response data
  cell.lines <- lapply(object, FUN = colnames)
  common.cell.lines <- sort(Reduce(intersect, cell.lines))
  object <- lapply(object, FUN = function(data.view){data.view[,common.cell.lines]})
  
  #Make the kernels for each data view, these are the complete kernels with both training and test data
  kernels <- lapply(object, FUN = function(data.view){getKernel(data.view, k)})
  kernels <- lapply(kernels, FUN = addNames, common.cell.lines)
  
  #Make an array (3 dimensions) with all the kernels
  kernel.array <- do.call(abind, c(kernels, list(along = 3)))
  mkkm.model <- runMKKM(kernel.array, clusters = C)#Model contains the complete combined kernel with all common cell lines
  combined.kernel <- mkkm.model$combined.kernel
  
  #Get the cell lines which are common to all data views..and which are common to 
  training.lines <- sort(intersect(rownames(drug.response), common.cell.lines))
 
  #subset the training data
  training.kernel.combined <-as.kernelMatrix(combined.kernel[training.lines,training.lines])
  drug.response.training <- drug.response[training.lines,]
  
  #empty object to store our svm models for each task
  drug.response.training <- drug.response.training[,-c(5, 12, 13, 21, 24, 26)]
  tasks <- colnames(drug.response.training)

  task.models <- list()
  task.preds <- list()
  for(i in 1:length(tasks)){
    print(task)
    #response <- as.matrix(drug.response.training[,task, drop = FALSE])
    response <- drug.response.training[,i]
    names(response) <- training.lines
    
    svm <- ksvm(x = training.kernel.combined, y = response, kernel = 'matrix', cross = 5)
    
    #index combined kernel matrix so that rows correspond to test cases and columns are the support vectors -- need to represent columns for
    #the test indices first so that indexes of support vectors are correct
    combined.kernel.test <- combined.kernel[setdiff(rownames(combined.kernel), training.lines), training.lines]
    combined.kernel.test <- combined.kernel.test[, SVindex(svm), drop = FALSE]
    combined.kernel.test <- as.kernelMatrix(combined.kernel.test)
    
    preds <- predict(object = svm, newdata = combined.kernel.test)
    
    task.preds[i] <- list(preds)
    task.models[i] <- svm
    
  }

  names(task.preds) <- tasks
  names(task.models) <- tasks

}

addNames <- function(kernel, common.cell.lines){
  rownames(kernel) <- common.cell.lines
  colnames(kernel) <- common.cell.lines
  
  return(kernel)
}

makeKernels <- function(object){
 
   kernel2Mat <- function(kernel){
    dim(kernel) <- c(sqrt(length(kernel)), sqrt(length(kernel)))
    return(kernel)
  }
  
  k <- 5
  kernels <- lapply(object, FUN = function(data.view){getKernel(data.view, k)})
  kernels <- lapply(kernels, FUN = as.numeric)
  kernels <- lapply(kernels, FUN = kernel2Mat)
  return(kernels)
}

getKernel <- function(data.view, k = 5){
  require(kernlab)
  require(parallel)
  require(FNN)
  
  #Calculate nearest neighbors
  neighbors <- FNN::get.knn(data = t(data.view), k = k)
  neighbor.strings <- apply(neighbors$nn.index,1,paste,collapse=" ")
  #Choose sigma 
  #sigma <- sigma
  #make kernel
  #kernel <- rbfdot(sigma = sigma)
  #kernel <- laplacedot()
  kernel  <- stringdot(type="string", length = 3)
  #transform data view into features
  #features <- kernelMatrix(kernel = kernel, x = t(data.view))
  #features <- kernelMatrix(kernel = kernel, x = neighbors$nn.index)
  #features <- kernelMatrix(kernel = kernel, x = neighbor.strings)
  features <- computeKernelMatrix(data = neighbor.strings, kernelFunction = kernel)
  
  return(features)
  
}



runMKKM <- function(kernel.array, clusters = 2, iter = 10){
  #source(fxn.location)
  require(MKKC)
  
  model <- list()
  
  #mehmet
  parameters <- list()
  parameters$cluster_count <- clusters
  parameters$iteration_count <- iter
  model$state <- mkkmeans_train(kernel.array, parameters)
  model$combined.kernel <- combineKernels(kernel.array, model$state$theta)
  
  return(model)
  
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
##################################################################################
combineKernels <- function(kernels, kernel.params){
  for(i in 1:length(kernel.params)){
    #print(i)
    kernels[, ,i] <- kernels[, ,i] * kernel.params[i]
  }
  return(apply(kernels, c(1,2), sum))
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