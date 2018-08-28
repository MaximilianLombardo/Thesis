############################################################################################################################################3
#Plotting functions

plotInternalQulaitymetrics <- function(quality.metrics, data.views.names){
  require(cowplot)
  
  sils <-lapply(quality.metrics, FUN = function(quality.metric){quality.metric$silhouette})
  names(sils) <- data.views.names
  dunns <-lapply(quality.metrics, FUN = function(quality.metric){quality.metric$dunn})
  names(dunns) <- data.views.names
  
  connectivities <-lapply(quality.metrics, FUN = function(quality.metric){quality.metric$connectivity})
  names(connectivities) <- data.views.names
  
  plot_grid(qualityTriPlot(sils, data.views.names),
               qualityTriPlot(dunns, data.views.names, y.lab = "Dunn Index"),
               qualityTriPlot(connectivities, data.views.names, y.lab = "Connectivity Score"), nrow = 1)
}

qualityTriPlot <- function(measures, data.views.names, x.lab = "Number of Clusters", y.lab = "Silhouette Score", main.title = "", max.ident = TRUE){
  require(ggplot2)
  require(data.table)
  require(dplyr)
  
  measure.dfs <- lapply(measures, FUN = function(measure){data.frame(measure, clusters = c(1:length(measure) + 1))})
  df <- rbindlist(measure.dfs, idcol = "Data.View")
  
  if(max.ident){
    extrema <- aggregate(measure ~ Data.View, data = df, max)
    extrema.pts <- which(df$measure %in% extrema$measure)
  }
  
  g <- ggplot(df, aes(x = clusters, y = measure, color = Data.View)) + geom_point() +
    geom_point(data=df[extrema.pts, ], aes(x=clusters, y=measure, color = Data.View), size=5, shape = 8)
  g <- g + xlab(x.lab) + ylab(y.lab) + ggtitle(main.title)

  g
  
  return(g)
}

plotSNFHeatComparison <- function(all.data, truelabel, down.pct = 1){
  require(RColorBrewer)
  
  #Breakout data
  W1 <- all.data[[1]]
  W2 <- all.data[[2]]
  W <- all.data[[3]]
  
  
  ## These similarity graphs have complementary information about clusters.
  new.palette=colorRampPalette(brewer.pal(9, "Spectral"),space="rgb")
  
  #Display Data View 1
  displayClusters(W1,truelabel, main.title = "Data View 1",
                  col = rev(new.palette(100)))
  #Display Data View 2
  displayClusters(W2,truelabel, main.title = "Data View 2",
                  col = rev(new.palette(100)))
  #Display fused kernel
  displayClusters(W, truelabel, main.title = "Fused Data View",
                  col = rev(new.palette(100)))
  #Display average kernel
  displayClusters((W1 + W2)/2, truelabel, main.title = "Average Data View",
                  col = rev(new.palette(100)))
  
  #Downsample points for big display
  samp.vec = c(1:nrow(W1))
  num.to.sample = floor(down.pct * nrow(W1))
  ind = sort(sample(samp.vec, num.to.sample))
  
  #Save old/default graphics parameters
  old.par <- par(mar = c(1, 1, 1, 1),
                 oma = c(10,0,10,0),
                 mfrow = c(1,1))
  
  #3 by 1 plot for figures
  plot.new()
  par(mfrow = c(1,3), oma = c(10,0,10,0))
  displayClusters(W1[ind,ind], truelabel[ind], main.title = "Data View 1",
                  col = rev(new.palette(100)))
  displayClusters(W2[ind,ind], truelabel[ind], main.title = "Data View 2",
                  col = rev(new.palette(100)))
  displayClusters(W[ind,ind], truelabel[ind], main.title = "Fused Data Views",
                  col = rev(new.palette(100)))
  title(main = list("Effect of Data fusion on Simulated Data", cex = 2), outer=TRUE)
  
  #Restore old graphics parameters
  par(old.par)
}

###Not Needed???
plotMKLKMHeatComparison <- function(Data1, Data2, truelabel,
                                    K = 20, alpha = 0.5, iter = 10, source.path){
  require(Rmosek)
  source(source.path)
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
  displayClusters(W1,truelabel, main.title = "Data View 1")
  displayClusters(W2,truelabel, main.title = "Data View 2")
  
  ## next, we fuse all the graphs
  ## then the overall matrix can be computed by similarity network fusion(SNF):
  W = SNF(Wall = list(W1,W2), K = K, t = iter)
  
  C = 2 					# number of clusters
  group = spectralClustering(W, C); 	# the final subtypes information
  group.1 = spectralClustering(W1, C); 	# the final subtypes information
  group.2 = spectralClustering(W2, C); 	# the final subtypes information
  
  ## you can evaluate the goodness of the obtained clustering results by calculate Normalized mutual information (NMI): if NMI is close to 1, it indicates that the obtained clustering is very close to the "true" cluster information; if NMI is close to 0, it indicates the obtained clustering is not similar to the "true" cluster information.
  
  ###Similarity network fusion result....
  plot.new()
  par(mfrow = c(1,3), oma = c(10,0,10,0))
  #displayClusters(W1, group.1, main.title = "Data View 1")
  #displayClusters(W2, group.2, main.title = "Data View 2")
  #displayClusters(W, group, main.title = "Fused Data Views")
  displayClusters(W1, truelabel, main.title = "Data View 1")
  displayClusters(W2, truelabel, main.title = "Data View 2")
  displayClusters(W, truelabel, main.title = "Fused Data Views")
  title(main = list("Effect of Data fusion on Simulated Data", cex = 2), outer=TRUE)
  
  #"Clustering results of individual and fused data views"
  
  SNFNMI = calNMI(group, truelabel)
  SNFNMI.1 = calNMI(group.1, truelabel)
  SNFNMI.2 = calNMI(group.2, truelabel)
  
  return(list(nmi.fused = SNFNMI, nmi.1 = SNFNMI.1, nmi.2 = SNFNMI.2))
  
}


####Scatter plots for class, color by class
plotSimDataScatter <- function(dat, lab, main.title){
  require(ggplot2)
  
  colnames(dat) <- c("V1", "V2")
  
  plt <- ggplot(data.frame(dat), aes(x = V1, y = V2, color = factor(lab))) + geom_point() 
  plt + ggtitle(main.title) + labs(colour = "True Identity")
}
