processSNFData <- function(dat1, dat2, label){
  
  data.modalities <- list()
  
  colnames(dat1) <- c("V1", "V2")
  colnames(dat2) <- c("V1", "V2")
  
  dat1 <- data.frame(dat1, label = label)
  dat2 <- data.frame(dat2, label = label)
  
  data.modalities$classes <- dat1
  data.modalities$inverted.classes <- dat2
  data.modalities$true.label <- label
  
  return(data.modalities)
}

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

switchLabels <- function(point, switch){
}

plotSNFHeatComparison <- function(Data1, Data2, truelabel,
                                  K = 20, alpha = 0.5, t = 10){
  require(SNFtool)
  
  #Process the data set
  Data1 <- as.matrix(Data1[,c("V1", "V2")])
  Data2 <- as.matrix(Data2[,c("V1", "V2")])
  
  ## Here, the simulation data (Data1, Data2) has two data types. They are complementary to each other. And two data types have the same number of points. The first half data belongs to the first cluster; the rest belongs to the second cluster.
   ##the ground truth of the simulated data;
  
  ## Calculate distance matrices(here we calculate Euclidean Distance, you can use other distance, e.g,correlation)
  ## If the data are all continuous values, we recommend the users to perform standard normalization before using SNF, though it is optional depending on the data the users want to use.  
  Data1 = standardNormalization(Data1)
  Data2 = standardNormalization(Data2)
  
  ## Calculate the pair-wise distance; If the data is continuous, we recommend to use the function "dist2" as follows; if the data is discrete, we recommend the users to use ""
  Dist1 = dist2(as.matrix(Data1),as.matrix(Data1));
  Dist2 = dist2(as.matrix(Data2),as.matrix(Data2));
  
  ## next, construct similarity graphs
  W1 = affinityMatrix(Dist1, K, alpha)
  W2 = affinityMatrix(Dist2, K, alpha)
  
  ## These similarity graphs have complementary information about clusters.
  displayClusters(W1,truelabel, main.title = "Data View 1");
  displayClusters(W2,truelabel, main.title = "Data View 2");
  
  ## next, we fuse all the graphs
  ## then the overall matrix can be computed by similarity network fusion(SNF):
  W = SNF(list(W1,W2), K, t)
  
  C = 2 					# number of clusters
  group = spectralClustering(W, C); 	# the final subtypes information
  group.1 = spectralClustering(W1, C); 	# the final subtypes information
  group.2 = spectralClustering(W2, C); 	# the final subtypes information
  
  ## you can evaluate the goodness of the obtained clustering results by calculate Normalized mutual information (NMI): if NMI is close to 1, it indicates that the obtained clustering is very close to the "true" cluster information; if NMI is close to 0, it indicates the obtained clustering is not similar to the "true" cluster information.
  
  ###Similarity network fusion result....
  plot.new()
  par(mfrow = c(1,3), oma = c(10,0,10,0))
  displayClusters(W1, group.1, main.title = "Data View 1")
  displayClusters(W2, group.2, main.title = "Data View 2")
  displayClusters(W, group, main.title = "Fused Data Views")
  title(main = list("Effect of Data fusion on Simulated Data", cex = 2), outer=TRUE)
  
  #"Clustering results of individual and fused data views"
  
  SNFNMI = calNMI(group, truelabel)
  SNFNMI.1 = calNMI(group.1, truelabel)
  SNFNMI.2 = calNMI(group.2, truelabel)
  
  return(list(nmi.fused = SNFNMI, nmi.1 = SNFNMI.1, nmi.2 = SNFNMI.2))
  
}

plotSimDataScatter <- function(dat, lab, main.title){
  require(ggplot2)
  
  colnames(dat) <- c("V1", "V2")
  
  plt <- ggplot(data.frame(dat), aes(x = V1, y = V2, color = factor(lab))) + geom_point() 
  plt + ggtitle(main.title) + labs(colour = "True Identity")
}
