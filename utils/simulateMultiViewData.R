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
