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
  
  kernels <- lapply(object, FUN = RBFKernel)
  
  ar <- array(kernels, c(nrow(kernels[[1]]), ncol(kernels[[1]]), length(kernels)));
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



RBFKernel <- function(data.view){
  require(kernlab)
  require(parallel)
  require(FNN)
  
  #Choose sigma 
  sigma <- sqrt(nrow(data.view))
  #make kernel
  kernel <- rbfdot(sigma = sigma)
  #transform data view into features
  features <- kernelMatrix(kernel = kernel, x = t(data.view))
  
  return(features)
  
}





