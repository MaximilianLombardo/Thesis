###########################################################################################################################
#Noisy Data Functions
simulateNoisyData <- function(num = 1000, gaussian.noise = 0.5, gamma.noise = 2){
  require(mlbench)

  data.modalities <- list()
  
  #Slightly overlapping data
  normals <- mlbench.2dnormals(n = num, cl = 2, r = 3, sd = 0.75)
  normals <- data.frame(V1 = normals$x[,1], V2=normals$x[,2], label=factor(normals$classes))
  
  noisy.normals.gaussian <- addGaussianNoise(normals, gaussian.noise)
  noisy.normals.gamma <- addGammaNoise(normals, gamma.noise)
  
  data.modalities$ground.truth <- normals
  data.modalities$gaussian.noise <- noisy.normals.gaussian
  data.modalities$gamma.noise <- noisy.normals.gamma
  data.modalities$true.label <- normals$label
    
  return(data.modalities)
}


addGaussianNoise <- function(pure.data, noise.level = 0.01){
  #2D only
  require(MASS)
  points <- nrow(pure.data)
  dims <- ncol(pure.data) - 1
  
  noise <- MASS::mvrnorm(points, mu = c(0,0), diag(x = noise.level, nrow = dims, ncol = dims))
  noise <- data.frame(V1 = noise[,1], V2 = noise[,2])
  
  noisy.data <- noise + pure.data[,-(ncol(pure.data))]
  noisy.data <- cbind(noisy.data, label = pure.data$label)
  
  return(noisy.data)
  
}

addGammaNoise <- function(pure.data, noise.level = 1){
  #2D only
  
  points <- nrow(pure.data)
  dims <- ncol(pure.data) - 1
  
  noise <- data.frame(V1 = rgamma(points, noise.level), V2 = rgamma(points, noise.level))
  
  noisy.data <- noise + pure.data[,-(ncol(pure.data))]
  noisy.data <- cbind(noisy.data, label = pure.data$label)
  
  return(noisy.data)
  
}
############################################################################################################################################
#Misclassification Simulation Data Functions
simulateMisclassificationData <- function(num = 1000, rad = 2.2){
  
  require(plyr)
  require(mlbench)
  data.modalities <- list()
  
  #Slightly overlapping data
  normals <- mlbench.2dnormals(n = num, cl = 2, r = 3, sd = 0.75)
  normals <- data.frame(V1 = normals$x[,1], V2=normals$x[,2], label=factor(normals$classes))
  
  inBoundary <- apply(normals, 1, checkBoundary, radius = rad)
  
  switch.labels <- normals[inBoundary,]
  keep.labels <- normals[!inBoundary,]
  
  #Keep track of real label
  real.label <- switch.labels$label
  
  #Biased misclassification of labels, always one class or the other...
  switch.labels.1 <- plyr::mapvalues(real.label, from = c(1), to = c(2))
  switch.labels.2 <- plyr::mapvalues(real.label, from = c(2), to = c(1))
  
  switch.labels.1 <- cbind(switch.labels[,c(1:2)], label = switch.labels.1)
  switch.labels.2 <- cbind(switch.labels[,c(1:2)], label = switch.labels.2)
  
  data.modalities$ground.truth <- normals
  data.modalities$misclassification.1 <- rbind(keep.labels, switch.labels.1)
  data.modalities$misclassification.2 <- rbind(keep.labels, switch.labels.2)
  data.modalities$true.label <- unlist(list(keep.labels$label, real.label))
  
  return(data.modalities)
}

checkBoundary <- function(point, radius){
  x <- as.numeric(point["V1"])
  y <- as.numeric(point["V2"])
  
  if(sqrt(x^2 + y^2) < radius){
    return(TRUE)
  }else{
    return(FALSE)
  }
}
###########################################################################################################
