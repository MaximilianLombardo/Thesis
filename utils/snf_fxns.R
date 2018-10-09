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


runSNFPipeline2 <- function(disease = "GBM", opt.iter, truelabel = NULL,
                           K = 20, alpha = 0.5, iter = 10,
                           down.pct = 1, C = 2){
  require(SNFtool)
  require(RColorBrewer)
  require(cluster)
  #require(rBayesianOptimization)
  
  #######
  #Load up the appropriate disease data views file paths
  data.views <- chooseDataType(disease)
  #Separate survival data file path
  survival <- data.views[4]
  #Read in data views and survival
  data.views <- lapply(data.views[1:3], FUN = read.table)
  names(data.views) <- c("ge", "meth", "mirna")
  
  survival <- read.table(survival, header = TRUE)
  
  #Process the data set
  data.views <- lapply(data.views, FUN = function(view){t(as.matrix(view))})
  
  #Normalization
  data.views <- lapply(data.views, FUN = standardNormalization)
  #########
  #Random Search hyperparameter optimization
  optimized.params <- randomSearchParams(data.views = data.views, search.iterations = opt.iter)
  
  K <- optimized.params$model.params$K
  alpha <- optimized.params$model.params$alpha
  iter <- optimized.params$model.params$iter
  C <- optimized.params$model.params$C

  #Calculate distance
  data.views.dist <- lapply(data.views, FUN = function(view){dist2(as.matrix(view), as.matrix(view))})
  
  ## next, construct similarity graphs
  graph.views <- lapply(data.views.dist, FUN = affinityMatrix, K, alpha)#This One
  
  #Run SNF
  graph.views$W_fused <- SNF(Wall = graph.views, K = K, t = iter)#This One
  
  #Cluster
  groupings = lapply(graph.views, FUN = function(graph.view){spectralClustering(graph.view, C)})#This one
  
  
  #Cast the previously calculated distance matrics to type dist
  data.views.dist <- lapply(data.views.dist, FUN = function(view.dist){as.dist(view.dist)})
  
  #Calculate silhouette scores...
  #Matching up single data views with respectively defined clusterings
  sils.single.views <- lapply(names(data.views),
                              FUN = function(idx){silhouette(x = groupings[[idx]],
                                                             dist = data.views.dist[[idx]])})
  names(sils.single.views) <- names(data.views)
  
  sil.multiview <- lapply(data.views.dist,
                          FUN = function(data.view.dist){silhouette(x = groupings$W_fused,
                                                                    dist = data.view.dist)})
  
  
  return(list(affinity.matrices = graph.views,
              silhouette.values = list(single.view = sils.single.views, multi.view = sil.multiview),
              identity = groupings, optimized.params = optimized.params))
}



randomSearchParams <- function(data.views,search.iterations,
                               K.range = c(2:10), alpha.range = c(0.1, 5.0),
                               iter.range = 10, C.range = c(2:10), num.pcs = 100){
  params <- list()
  params$best.score <- numeric(0)
  
  i <- 1
  
  while(i < search.iterations){
    print(i)
    #Sample new set of params
    #K <- sample(x = K.range, size = 1)
    K = 3
    alpha = runif(n =1, min = alpha.range[1], max = alpha.range[2])
    iter = iter.range
    C <- sample(x = C.range, size = 1)
    
    print(list(K = K, alpha = alpha, iter = iter, C = C))
    
    model.res <- try(maximizeSilhouetteSNF(data.views = data.views,
                                       K = K, alpha = alpha,
                                       iter = iter, C = C))
    if(i == 1){
      #On first iterations
      params$best.score <- model.res$Score
      params$model.params <- list(K = K, alpha = alpha,
                                  iter = iter, C = C,
                                  num.pcs = num.pcs)
    }else{
      #The rest of the iters
      if(model.res$Score > params$best.score & is.list(model.res)){
        params$best.score <- model.res$Score
        params$model.params <- list(K = K, alpha = alpha,
                                    iter = iter, C = C,
                                    num.pcs = num.pcs)
      }
    }
    i = i + 1
  }
  
  return(params)
  
}

maximizeSNFGBM <- function(K, alpha, iter, C){
  #require(caret)
  #require(rpart)
  
  model.results <- list()
  #LoadData
  data.views <- chooseDataType(disease)
  #Separate survival data file path
  survival <- data.views[4]
  #Read in data views and survival
  data.views <- lapply(data.views[1:3], FUN = read.table)
  names(data.views) <- c("ge", "meth", "mirna")
  
  survival <- read.table(survival, header = TRUE)
  
  #Process the data set
  data.views <- lapply(data.views, FUN = function(view){t(as.matrix(view))})
  
  #Normalization
  data.views <- lapply(data.views, FUN = standardNormalization)
  
  #PCA
  #data.views.pca <- lapply(data.views, FUN = function(view){prcomp(t(view), rank. = num.pcs)$rotation})
  
  #Calculate distance
  data.views.dist <- lapply(data.views, FUN = function(view){dist2(as.matrix(view), as.matrix(view))})
  
  ## next, construct similarity graphs
  graph.views <- lapply(data.views.dist, FUN = affinityMatrix, K, alpha)#This One
  
  #Run SNF
  graph.views$W_fused <- try(SNF(Wall = graph.views, K = K, t = iter))#This One
  
  #Cluster
  groupings = lapply(graph.views, FUN = function(graph.view){spectralClustering(graph.view, C)})#This one
  
  #Cast the previously calculated distance matrics to type dist
  data.views.dist <- lapply(data.views.dist, FUN = function(view.dist){as.dist(view.dist)})
  
  #Calculate silhouette scores...
  #Matching up single data views with respectively defined clusterings
  sils.single.views <- lapply(names(data.views),
                              FUN = function(idx){silhouette(x = groupings[[idx]],
                                                             dist = data.views.dist[[idx]])})
  names(sils.single.views) <- names(data.views)
  
  sil.multiview <- lapply(data.views.dist,
                          FUN = function(data.view.dist){silhouette(x = groupings$W_fused,
                                                                    dist = data.view.dist)})
  
  #Cast the previously calculated distance matrics to type dist
  data.views.dist <- lapply(data.views.dist, FUN = function(view.dist){as.dist(view.dist)})
  
  #Calculate silhouette scores...
  #Matching up single data views with respectively defined clusterings
  sils.single.views <- lapply(names(data.views),
                              FUN = function(idx){silhouette(x = groupings[[idx]],
                                                             dist = data.views.dist[[idx]])})
  names(sils.single.views) <- names(data.views)
  
  sil.multiview <- lapply(data.views.dist,
                          FUN = function(data.view.dist){silhouette(x = groupings$W_fused,
                                                                    dist = data.view.dist)})
  avg.sils <- lapply(sil.multiview, FUN = function(sil.score){mean(sil.score[,"sil_width"])})
  
  model.results$Score <- mean(unlist(avg.sils))
  model.results$Pred <- 0#groupings$W_fused
  
  return(model.results)
}

maximizeSilhouetteSNF <- function(data.views, K, alpha, iter, C){
  #require(caret)
  #require(rpart)
  
  model.results <- list()
  
  #Calculate distance
  data.views.dist <- lapply(data.views, FUN = function(view){dist2(as.matrix(view), as.matrix(view))})
  
  ## next, construct similarity graphs
  graph.views <- lapply(data.views.dist, FUN = affinityMatrix, K, alpha)#This One
  
  #Run SNF
  graph.views$W_fused <- SNF(Wall = graph.views, K = K, t = iter)#This One
  
  #Cluster
  groupings = lapply(graph.views, FUN = function(graph.view){spectralClustering(graph.view, C)})#This one
  
  #Cast the previously calculated distance matrics to type dist
  data.views.dist <- lapply(data.views.dist, FUN = function(view.dist){as.dist(view.dist)})
  
  #Calculate silhouette scores...
  #Matching up single data views with respectively defined clusterings
  sils.single.views <- lapply(names(data.views),
                              FUN = function(idx){silhouette(x = groupings[[idx]],
                                                             dist = data.views.dist[[idx]])})
  names(sils.single.views) <- names(data.views)
  
  sil.multiview <- lapply(data.views.dist,
                          FUN = function(data.view.dist){silhouette(x = groupings$W_fused,
                                                                    dist = data.view.dist)})
  
  #Cast the previously calculated distance matrics to type dist
  data.views.dist <- lapply(data.views.dist, FUN = function(view.dist){as.dist(view.dist)})
  
  #Calculate silhouette scores...
  #Matching up single data views with respectively defined clusterings
  sils.single.views <- lapply(names(data.views),
                              FUN = function(idx){silhouette(x = groupings[[idx]],
                                                             dist = data.views.dist[[idx]])})
  names(sils.single.views) <- names(data.views)
  
  sil.multiview <- lapply(data.views.dist,
                          FUN = function(data.view.dist){silhouette(x = groupings$W_fused,
                                                                    dist = data.view.dist)})
  avg.sils <- lapply(sil.multiview, FUN = function(sil.score){mean(sil.score[,"sil_width"])})
  
  model.results$Score <- mean(unlist(avg.sils))
  model.results$Pred <- 0#groupings$W_fused
  
  return(model.results)
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
