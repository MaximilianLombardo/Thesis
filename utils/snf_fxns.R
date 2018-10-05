##############################################################################################################3
runSNFPipeline <- function(Data1, Data2, truelabel,
                           K = 20, alpha = 0.5, iter = 10,
                           down.pct = 1, C = 2){
  require(SNFtool)
  require(RColorBrewer)
  #Process the data set
  Data1 <- as.matrix(Data1[,c("V1", "V2")])
  Data2 <- as.matrix(Data2[,c("V1", "V2")])
  
  #Normalization
  Data1 = standardNormalization(Data1)
  Data2 = standardNormalization(Data2)
  
  #Calculate distance
  Dist1 = dist2(as.matrix(Data1),as.matrix(Data1));
  Dist2 = dist2(as.matrix(Data2),as.matrix(Data2));
  
  ## next, construct similarity graphs
  W1 = affinityMatrix(Dist1, K, alpha)
  W2 = affinityMatrix(Dist2, K, alpha)
  
  #Run SNF
  W = SNF(Wall = list(W1,W2), K = K, t = iter)
  
  
  #Cluster
  group = spectralClustering(W, C);
  group.1 = spectralClustering(W1, C);
  group.2 = spectralClustering(W2, C);
  
  #"Clustering results of individual and fused data views"
  SNFNMI = calNMI(group, truelabel)
  SNFNMI.1 = calNMI(group.1, truelabel)
  SNFNMI.2 = calNMI(group.2, truelabel)
  
  return(list(affinity.matrices = list(W1 = W1, W2 = W2, W = W),
              nmi.values = list(nmi.1 = SNFNMI.1, nmi.2 = SNFNMI.2, nmi.fused = SNFNMI),
              identity = list(ident.1 = group.1, ident.2 = group.2, ident.fused = group)))
}


runSNFPipeline2 <- function(data.views, truelabel = NULL,
                           K = 20, alpha = 0.5, iter = 10,
                           down.pct = 1, C = 2){
  require(SNFtool)
  require(RColorBrewer)
  
  #Temporary
  data.views <- chooseDataType("GBM")
  data.views <- lapply(data.views, FUN = read.table)
  names(data.views) <- c("ge", "meth", "mirna", "survival")
  
  survival <- data.views$survival
  data.views <- data.views[-4]
  
  #Process the data set
  data.views <- lapply(data.views, FUN = function(view){t(as.matrix(view))})
  
  #Normalization
  data.views <- lapply(data.views, FUN = standardNormalization)
  
  #Calculate distance
  data.views <- lapply(data.views, FUN = function(view){dist2(as.matrix(view), as.matrix(view))})

  ## next, construct similarity graphs
  data.views <- lapply(data.views, FUN = affinityMatrix, K, alpha)
  
  #Run SNF
  data.views$W <- SNF(Wall = data.views, K = K, t = iter)
  
  #Cluster
  group = spectralClustering(W, C);
  group.1 = spectralClustering(W1, C);
  group.2 = spectralClustering(W2, C);
  
  #"Clustering results of individual and fused data views"
  SNFNMI = calNMI(group, truelabel)
  SNFNMI.1 = calNMI(group.1, truelabel)
  SNFNMI.2 = calNMI(group.2, truelabel)
  
  return(list(affinity.matrices = list(W1 = W1, W2 = W2, W = W),
              nmi.values = list(nmi.1 = SNFNMI.1, nmi.2 = SNFNMI.2, nmi.fused = SNFNMI),
              identity = list(ident.1 = group.1, ident.2 = group.2, ident.fused = group)))
}




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
