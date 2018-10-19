runDreamPipeline <- function(object, bmtmkl.dir = "~/Documents/gitRepos/master/utils/"){
  #Sourcing the functions needed for bayesian multitask multiple kernel learning
  source(paste0(bmtmkl.dir, "bayesian_multitask_multiple_kernel_learning_train.R"))
  source(paste0(bmtmkl.dir, "bayesian_multitask_multiple_kernel_learning_test.R"))
  require(abind)
  #####
  #Split up the task matrix (the target matrix) and the input data matrices
  drug.response <- object$Drug_Response_Training
  object <- object[-1]
  
  #For the data matrices, get the common cell lines and subset the data
  cell.lines <- lapply(object, FUN = colnames)
  common.cell.lines <- sort(Reduce(intersect, cell.lines))
  common.cell.lines <- sort(intersect(rownames(drug.response), common.cell.lines))
  
  drug.response <- drug.response[common.cell.lines,]
  object <- lapply(object, FUN = function(data.view){data.view[,common.cell.lines]})
  
  ##############################################################
  #Normalization is already done, but we will need to construct the kernel matrices and set as arrays
  k <- 5
  kernels <- lapply(object, FUN = function(data.view){getKernel(data.view, k)})
  
  #ar <- array(kernels, c(nrow(kernels[[1]]), ncol(kernels[[1]]), length(kernels)));
  kernel.array <- do.call(abind, c(kernels, list(along = 3)))
  
  
  #######################Training for bayesian multitask multiple kernel learning
 
  #initalize the parameters 
  parameters <- list()
  
  #set the hyperparameters of gamma prior used for sample weights
  parameters$alpha_lambda <- 1
  parameters$beta_lambda <- 1
  
  #set the hyperparameters of gamma prior used for intermediate noise
  parameters$alpha_upsilon <- 1
  parameters$beta_upsilon <- 1
  
  #set the hyperparameters of gamma prior used for bias
  parameters$alpha_gamma <- 1e-10
  parameters$beta_gamma <- 1e-10
  
  #set the hyperparameters of gamma prior used for kernel weights
  parameters$alpha_omega <- 1
  parameters$beta_omega <- 1
  
  #set the hyperparameters of gamma prior used for output noise
  parameters$alpha_epsilon <- 1
  parameters$beta_epsilon <- 1
  
  #set the number of iterations
  parameters$iteration <- 200
  
  #determine whether you want to calculate and store the lower bound values
  parameters$progress <- 0
  
  #set the seed for random number generator used to initalize random variables
  parameters$seed <- 1606
  
  #set the number of tasks (e.g., the number of compounds in Nature Biotechnology paper)
  tasks <- ncol(drug.response)
  #set the number of kernels (e.g., the number of views in Nature Biotechnology paper)
  P <- length(object)
    
  #initialize the kernels and outputs of each task for training
  Ktrain <- vector("list", tasks)
  ytrain <- vector("list", tasks)
  for (task in 1:tasks) {
    Ktrain[[task]] <- kernel.array #should be an Ntra x Ntra x P matrix containing similarity values between training samples of task t
      ytrain[[task]] <- drug.response[, task, drop = FALSE] #should be an Ntra x 1 matrix containing target outputs of task t
  }
  
  #perform training
  state <- bayesian_multitask_multiple_kernel_learning_train(Ktrain, ytrain, parameters)
  
  #display the kernel weights
  print(state$be$mu[(tasks+1):(tasks+P)])
  
  #initialize the kernels of each task for testing
  Ktest <- vector("list", tasks)
  for (t in 1:tasks) {
    Ktest[[t]] <- ?? #should be an Ntra x Ntest x P matrix containing similarity values between training and test samples of task t
  }
  
  #perform prediction
  prediction <- bayesian_multitask_multiple_kernel_learning_test(Ktest, state)
  
  
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
  features <- kernelMatrix(kernel = kernel, x = neighbor.strings)
  
  
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
