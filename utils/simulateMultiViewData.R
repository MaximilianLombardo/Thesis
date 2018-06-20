simulateMultiviewData(gaussian.noise = 0.5, gamma.noise = 2){
  require(mlbench)
  require(ggplot2)
  
  data.modalities <- list()
  
  #Slightly overlapping data
  normals <- mlbench.2dnormals(n = 1000, )
  normals <- data.frame(x = normals$x[,1], y=normals$x[,2], class=factor(normals$classes))
  
  noisy.normals.gaussian <- addGaussianNoise(normals, gaussian.noise)
  noisy.normals.gamma <- addGammaNoise(normals, gamma.noise)
  
  data.modalities$gaussian.noise <- noisy.normals.gaussian
  data.modalities$gamma.noise <- noisy.normals.gamma
  
  return(data.modalities)
}


addGaussianNoise <- function(pure.data, noise.level = 0.01){
  #2D only
  require(MASS)
  points <- nrow(pure.data)
  dims <- ncol(pure.data) - 1
  
  noise <- MASS::mvrnorm(points, mu = c(0,0), diag(x = noise.level, nrow = dims, ncol = dims))
  noise <- data.frame(x = noise[,1], y = noise[,2])
  
  noisy.data <- noise + pure.data[,-(ncol(pure.data))]
  noisy.data <- cbind(noisy.data, class = pure.data$class)
  
  return(noisy.data)
  
}

addGammaNoise <- function(pure.data, noise.level = 1){
  #2D only
  
  points <- nrow(pure.data)
  dims <- ncol(pure.data) - 1
  
  noise <- data.frame(x = rgamma(points, noise.level), y = rgamma(points, noise.level))
  
  noisy.data <- noise + pure.data[,-(ncol(pure.data))]
  noisy.data <- cbind(noisy.data, class = pure.data$class)
  
  return(noisy.data)
  
}

spiral = mlbench.spirals(1000,sd=0.1)
spiral = data.frame(x=spiral$x[,1],y=spiral$x[,2],class=factor(spiral$classes))

ggplot(spiral,aes(x,y,color=class)) + geom_point()


#Slightly overlapping data
normals <- mlbench.2dnormals(n = 1000)
normals <- data.frame(x = normals$x[,1], y=normals$x[,2], class=factor(normals$classes))


#Progressively noisier overlapping data...gaussian noise
ggplot(normals,aes(x,y,color=class)) + geom_point()
ggplot(addGaussianNoise(normals),aes(x,y,color=class)) + geom_point()
ggplot(addGaussianNoise(normals, noise.level = 0.5),aes(x,y,color=class)) + geom_point()

#Progressively noisier overlapping data...gamma noise
ggplot(normals,aes(x,y,color=class)) + geom_point()
ggplot(addGammaNoise(normals),aes(x,y,color=class)) + geom_point()
ggplot(addGammaNoise(normals, noise.level = 2),aes(x,y,color=class)) + geom_point()



#Show that PCA on the concatenated feature set yields garbage.....
concatViews <- cbind(addGaussianNoise(normals, noise.level = 0.5)[,-3], addGammaNoise(normals, noise.level = 2)[,-3])
concatViews.scale <- scale(concatViews, center = TRUE, scale = TRUE)
concat.pca <- prcomp(x = t(concatViews.scale))
concat.pca <- data.frame(x = concat.pca$rotation[,1], y = concat.pca$rotation[,3], class = normals$class)
ggplot(concat.pca,aes(x,y,color=class)) + geom_point()

#Do the same on the TSNE plot?
##############################################################################################3


library(SNFtool)


#Generate plots with simulation data


## First, set all the parameters:
K = 20;		# number of neighbors, usually (10~30)
alpha = 0.5;  	# hyperparameter, usually (0.3~0.8)
T = 10; 	# Number of Iterations, usually (10~20)

## Data1 is of size n x d_1, where n is the number of patients, d_1 is the number of genes, e.g.
## Data2 is of size n x d_2, where n is the number of patients, d_2 is the number of methylation, e.g.
data(Data1)
data(Data2)


plotSNFHeatComparison <- function(Data1, Data2, truelabel = c(matrix(1,100,1),matrix(2,100,1)),
                                  K = 20, alpha = 0.5, T = 10){
  
}

## Here, the simulation data (Data1, Data2) has two data types. They are complementary to each other. And two data types have the same number of points. The first half data belongs to the first cluster; the rest belongs to the second cluster.

truelabel = c(matrix(1,100,1),matrix(2,100,1)) ##the ground truth of the simulated data;


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
displayClusters(W1,truelabel);
displayClusters(W2,truelabel);

## next, we fuse all the graphs
## then the overall matrix can be computed by similarity network fusion(SNF):
W = SNF(list(W1,W2), K, T)




C = 2 					# number of clusters
group = spectralClustering(W, C); 	# the final subtypes information
group.1 = spectralClustering(W1, C); 	# the final subtypes information
group.2 = spectralClustering(W2, C); 	# the final subtypes information

## you can evaluate the goodness of the obtained clustering results by calculate Normalized mutual information (NMI): if NMI is close to 1, it indicates that the obtained clustering is very close to the "true" cluster information; if NMI is close to 0, it indicates the obtained clustering is not similar to the "true" cluster information.

###Similarity network fusion result....
plot.new()
par(mfrow = c(1,3), oma = c(10,0,10,0))
displayClusters(W1, group.1, main.title = 'Data View 1')
displayClusters(W2, group.2, main.title = 'Data View 2')
displayClusters(W, group, main.title = 'Fused Data Views')
title(main = list("Title", cex = 4), outer=TRUE)

"Clustering results of individual and fused data views"

SNFNMI = calNMI(group, truelabel)
SNFNMI.1 = calNMI(group.1, truelabel)
SNFNMI.2 = calNMI(group.2, truelabel)


plotSNFHeatComparison <- function(Data1, Data2, truelabel = c(matrix(1,100,1),matrix(2,100,1)),
                                  K = 20, alpha = 0.5, T = 10){
  
  ## Here, the simulation data (Data1, Data2) has two data types. They are complementary to each other. And two data types have the same number of points. The first half data belongs to the first cluster; the rest belongs to the second cluster.
  truelabel = c(matrix(1,100,1),matrix(2,100,1)) ##the ground truth of the simulated data;
  
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
  displayClusters(W1,truelabel);
  displayClusters(W2,truelabel);
  
  ## next, we fuse all the graphs
  ## then the overall matrix can be computed by similarity network fusion(SNF):
  W = SNF(list(W1,W2), K, T)
  
  C = 2 					# number of clusters
  group = spectralClustering(W, C); 	# the final subtypes information
  group.1 = spectralClustering(W1, C); 	# the final subtypes information
  group.2 = spectralClustering(W2, C); 	# the final subtypes information
  
  ## you can evaluate the goodness of the obtained clustering results by calculate Normalized mutual information (NMI): if NMI is close to 1, it indicates that the obtained clustering is very close to the "true" cluster information; if NMI is close to 0, it indicates the obtained clustering is not similar to the "true" cluster information.
  
  ###Similarity network fusion result....
  plot.new()
  par(mfrow = c(1,3), oma = c(10,0,10,0))
  displayClusters(W1, group.1, main.title = 'Data View 1')
  displayClusters(W2, group.2, main.title = 'Data View 2')
  displayClusters(W, group, main.title = 'Fused Data Views')
  title(main = list("Title", cex = 4), outer=TRUE)
  
  #"Clustering results of individual and fused data views"
  
  SNFNMI = calNMI(group, truelabel)
  SNFNMI.1 = calNMI(group.1, truelabel)
  SNFNMI.2 = calNMI(group.2, truelabel)
  
  return(list(nmi.fused = SNFNMI, nmi.1 = SNFNMI.1, nmi.2 = SNFNMI.2))
  
}



