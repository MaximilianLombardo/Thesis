#Exploratory plots of simulated data


plotSimulatedData1 <- function(){
  require(SNFtool)
  require(ggplot2)
  require(cowplot)
  
  
  ####Plotting the 
  data("Data1")
  data("Data2")
  
  truelabel <- c(matrix(1,100,1),matrix(2,100,1))
  
  Data1 <- cbind(Data1, truelabel)
  Data2 <- cbind(Data2, truelabel)
  
  plot.new()
  p1 <- ggplot(Data1, aes(V1, V2)) + geom_point(aes(colour = factor(truelabel)))
  p1 <- p1 + scale_color_discrete(name ="Type", labels=c("Class 1", "Class 2"))
  
  p2 <- ggplot(Data2, aes(V3, V4)) + geom_point(aes(colour = factor(truelabel)))
  p2 <- p2 + scale_color_discrete(name ="Type", labels=c("Class 1", "Class 2"))
  
  plot_grid(p1, p2, labels = "AUTO")
}


plotSimulatedData2 <- function(){
  require(ggplot2)
  
  dataViews <- simulateMultiviewData(n.points = 200, st.dev = 0.25, radius = 1,
                                     gaussian.noise = 0.5, gamma.noise = 0.5)
  
  plot.new()
  
  p <- ggplot(dataViews$ground.truth, aes(x, y)) + geom_point(aes(colour = factor(class)))
  p <- p + scale_color_discrete(name ="Type", labels=c("Class 1", "Class 2"))
  
  p1 <- ggplot(dataViews$gaussian.noise, aes(x, y)) + geom_point(aes(colour = factor(class)))
  p1 <- p1 + scale_color_discrete(name ="Type", labels=c("Class 1", "Class 2"))
  
  p2 <- ggplot(dataViews$gamma.noise, aes(x, y)) + geom_point(aes(colour = factor(class)))
  p2 <- p2 + scale_color_discrete(name ="Type", labels=c("Class 1", "Class 2"))
  
  plot_grid(p, p1, p2, labels = "AUTO", nrow = 1)
  
  
}


simulateMultiviewData <- function(n.points = 200 ,gaussian.noise = 0.5, gamma.noise = 2, radius = sqrt(2), st.dev = 0.5){
  require(mlbench)
  require(ggplot2)
  
  data.modalities <- list()
  
  #Slightly overlapping data
  normals <- mlbench.2dnormals(n = n.points, cl=2, sd = st.dev, r = radius )
  normals <- data.frame(x = normals$x[,1], y=normals$x[,2], class=factor(normals$classes))
  
  noisy.normals.gaussian <- addGaussianNoise(normals, gaussian.noise)
  noisy.normals.gamma <- addGammaNoise(normals, gamma.noise)
  
  data.modalities$ground.truth <- normals
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