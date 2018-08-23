chooseDataType <- function(disease){
  
  #TODO ENTER BASE STRING AS PARAMETER
  
  if(disease == "GBM"){
    data.files <- list.files("~/Documents/gitRepos/master/data/SNF/GBM/", full.names = TRUE)
  }else if(disease == "BREAST"){
    data.files <- list.files("~/Documents/gitRepos/master/data/SNF/Breast/", full.names = TRUE)
  }else if(disease == "COLON"){
    data.files <- list.files("~/Documents/gitRepos/master/data/SNF/Colon/", full.names = TRUE)
  }else if (disease == "KIDNEY"){
    data.files <- list.files("~/Documents/gitRepos/master/data/SNF/Kidney/", full.names = TRUE)
  }else if(disease == "LUNG"){
    data.files <- list.files("~/Documents/gitRepos/master/data/SNF/Lung/", full.names = TRUE)
  }else{
    print("data for specified disease does not exist")
  }
  
  return(data.files)
  
}



computeKernelMatrix <- function(data, kernelFunction){
  
  kernel <- matrix(,nrow = ncol(data), ncol = ncol(data))
  
  row.names(kernel) <- colnames(data)
  colnames(kernel) <- colnames(data)
  
  for(i in c(1:ncol(data))){
    for(j in c(1:ncol(data))){
      
      kernel[i,j] <- kernelFunction(data[,i], data[,j])
      
    }
  }
  return(kernel)
}




computeKernelMatrix1D <- function(data){
  
  kernel <- matrix(,nrow = length(data$Survival), ncol = length(data$Survival))
  
  row.names(kernel) <- data$PatientID
  colnames(kernel) <- data$PatientID
  
  for(i in c(1:length(data$Survival))){
    for(j in c(1:length(data$Survival))){
      
      kernel[i,j] <- sqrt((data$Survival[i] - data$Survival[j])^2)
      
    }
  }
  # kernel <- abs(kernel - 1)
  return(1 - (kernel/max(kernel)))
}





kernel.distance <- function(K1, K2){
  K <- abs(K1 - K2)
  return(norm(K, 'f'))
}

return.best.kernel <- function(kernel.variations, target.kernel){
  distances <- unlist(lapply(kernel.variations, FUN = kernel.distance, target.kernel))
  return(kernel.variations[[which.min(distances)]])
}

#Generate kernels with different hyperparameters for selection...

generate.kernels <- function(data){
  kernels <- list()
  #sigmas <- c(1:2 %o% 10^-(0:5))
  sigmas <- c(1:9 %o% 10^-(4:5))
  for(i in 1:length(sigmas)){
    #print(i)
    rbf <- rbfdot(sigma = sigmas[i])
    kernels[[i]] <- computeKernelMatrix(data, rbf)
  }
  return(kernels)
}


kernelTargetAlignment <- function(data.view, survival.similarity){
  
  data.view.kernels <- generate.kernels(data.view)
  data.view.kern <- return.best.kernel(data.view.kernels, survival.similarity)
  data.view.kern.kl <- as.kernelMatrix(data.view.kern)
  
  return(data.view.kern.kl)
}