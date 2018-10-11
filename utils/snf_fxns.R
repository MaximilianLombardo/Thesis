##############################################################################################################3
chooseDataType <- function(disease, root.dir = "~/Documents/gitRepos/master/data/SNF/"){
  
  if(disease == "GBM"){
    data.files <- list.files(paste0(root.dir, disease), full.names = TRUE)
  }else if(disease == "BREAST"){
    data.files <- list.files(paste0(root.dir, disease), full.names = TRUE)
  }else if(disease == "COLON"){
    data.files <- list.files(paste0(root.dir, disease), full.names = TRUE)
  }else if (disease == "KIDNEY"){
    data.files <- list.files(paste0(root.dir, disease), full.names = TRUE)
  }else if(disease == "LUNG"){
    data.files <- list.files(paste0(root.dir, disease), full.names = TRUE)
  }else{
    print("data for specified disease does not exist")
  }
  
  return(data.files)
  
}


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


processRealData <- function(disease){
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
  
  #TODO Find variable gene function use for pca....
  
  #PCA
  data.views.pca <- lapply(data.views, FUN = function(view){prcomp(t(view), rank. = 100)$rotation})
  return(list(data.views = data.views, views.pca = data.views.pca, survival = survival))
}

runSNFPipelineRealData <- function(object, opt.iter, truelabel = NULL,
                            K = 20, alpha = 0.5, iter = 10,
                            down.pct = 1, C = 2, num.pcs){
  require(SNFtool)
  require(RColorBrewer)
  require(cluster)
  require(rBayesianOptimization)
  
  data.views <- object$data.views
  views.pca <- object$views.pca
  survival <- object$survival
  
  ######
  #optimize individual data views
  
  opt.params <- list()
  
  for(idx in c(1:length(data.views))){
    view.pca <- views.pca[[idx]]
    
    #opt.params[[idx]] <- BayesianOptimization(FUN = maximizeSingleView,
    #                                          bounds = list(K = c(10L,50L), alpha = c(0.1, 5.0),
    #                                                        C = c(2L,15L),  iter = c(10,10),
    #                                                        num.pcs = c(5L,100L)),
    #                                          n_iter = 10, init_points = 10, acq = 'ucb',
    #                                          init_grid_dt = NULL,
    #                                          eps = 0.1)
    
    opt.params[[idx]] <- randomSearchParams(search.iterations = 500)
    
  }
  
  names(opt.params) <- names(views.pca)
  
  ######
  
  
  
  
  #Calculate distance
  data.views.dist <- lapply(data.views.pca, FUN = function(view){dist2(as.matrix(view), as.matrix(view))})
  
  ## next, construct similarity graphs
  graph.views <- lapply(data.views.dist, FUN = affinityMatrix, K, alpha)#This One
  
  #Run SNF
  graph.views$W_fused <- SNF(Wall = graph.views, K = K, t = iter)#This One
  
  #Cluster
  groupings = lapply(graph.views, FUN = function(graph.view){spectralClustering(affinity = graph.view, K = C)})#This one
  
  
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
              identity = groupings))
}

maximizeSingleView <- function(K, alpha, C, iter, num.pcs){
  require(SNFtool)
  require(fcd)
  require(anocva)
  
  model.results <- list()
  
  #Calculate distance
  view.dist <- dist2(as.matrix(view.pca[,1:num.pcs]),
                      as.matrix(view.pca[,1:num.pcs]))
  
  ## next, construct similarity graphs
  graph.view <- try(affinityMatrix(view.dist, K, alpha))
  
  #Cluster
  #grouping <- try(spectralClustering(affinity = graph.view, K = C, type = 1))
  #grouping <- try(spectral.clustering(A = graph.view, K = C))
  grouping <- anocva::spectralClustering(W = graph.view, k = C)
  
  if(is.integer(grouping)){
    #Cast the previously calculated distance matrics to type dist
    view.dist <- as.dist(view.dist)
    
    #Calculate silhouette scores...
    #Matching up single data views with respectively defined clusterings
    sil <- silhouette(x = grouping, dist = view.dist)
  
    
    model.results$Score <- mean(sil[,"sil_width"])
    model.results$Pred <- 0#groupings$W_fused
  
  return(model.results)
  }else{
    model.results$Score <- 0
    model.results$Pred <- 0
  }
  
}


discriminativeFeatureSelection <- function(data.views, ident){
  require(SNFtool)
  
  calFeatureNMI <- function(data.view, ident){
    #return(apply(data.view, 2, FUN = function(gene){calNMI(x = kmeans(gene, 10)$cluster, y = ident)}))
    features.mi <- apply(data.view, 2, FUN = function(gene){calNMI(x = round(gene), y = ident)})
    return(rev(sort(features.mi)))
  }
  
  data.views.sorted.feats <- lapply(data.views, FUN = calFeatureNMI, ident)
  
  return(data.views.sorted.feats)
}

runSNFPipeline2 <- function(disease = "GBM", opt.iter, truelabel = NULL,
                           K = 20, alpha = 0.5, iter = 10,
                           down.pct = 1, C = 2, num.pcs){
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
  
  
  #PCA
  data.views.pca <- lapply(data.views, FUN = function(view){prcomp(t(view), rank. = num.pcs)$rotation})
  
  #Calculate distance
  data.views.dist <- lapply(data.views.pca, FUN = function(view){dist2(as.matrix(view), as.matrix(view))})
  
  #data.views.dist <- lapply(data.views.pca,
  #                          FUN = function(data.view){as.matrix(dist(data.view,
  #                                                                   method = "manhattan"))})
  
  #########
  #Random Search hyperparameter optimization
  #optimized.params <- randomSearchParams(data.views = data.views, search.iterations = opt.iter)
  
  
  #optimized.params.bae <- BayesianOptimization(FUN = maximizeSNFGBM,
  #                                             bounds = list(K = c(10L,50L), alpha = c(0.1, 5.0),
  #                                                           C = c(2L,15L),  iter = c(10,10)),
  #                                             n_iter = 49, init_points = 150, acq = 'ei',
  #                                             eps = 0.1, kernel = list(type = "matern"))
  
  #maximizeSNFGBM()
  
  #K <- optimized.params$model.params$K
  #alpha <- optimized.params$model.params$alpha
  #iter <- optimized.params$model.params$iter
  #C <- optimized.params$model.params$C

  #Calculate distance
  #data.views.dist <- lapply(data.views, FUN = function(view){dist2(as.matrix(view), as.matrix(view))})
  
  ## next, construct similarity graphs
  graph.views <- lapply(data.views.dist, FUN = affinityMatrix, K, alpha)#This One
  
  #Run SNF
  graph.views$W_fused <- SNF(Wall = graph.views, K = K, t = iter)#This One
  
  #Cluster
  groupings = lapply(graph.views, FUN = function(graph.view){spectralClustering(affinity = graph.view, K = C)})#This one
  
  
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
              identity = groupings))
}



randomSearchParams <- function(search.iterations,
                               K.range = c(5:100), alpha.range = c(0.1, 5.0),
                               iter.range = 10, C.range = c(2:10), num.pcs.range = c(1:99)){
  params <- list()
  params$best.score <- numeric(0)
  
  i <- 1
  
  while(i < search.iterations){
    K <- sample(x = K.range, size = 1)
    alpha <- runif(n =1, min = alpha.range[1], max = alpha.range[2])
    iter <- iter.range
    C <- sample(x = C.range, size = 1)
    num.pcs <- sample(x = num.pcs.range, size = 1)
    
    model.res <- try(maximizeSingleView(K = K, alpha = alpha,
                                        iter = iter, C = C, num.pcs = num.pcs))
    if(i == 1){
      #On first iterations
      params$best.score <- model.res$Score
      params$model.params <- list(K = K, alpha = alpha,
                                  C = C, iter = iter, num.pcs)
    }else{
      #The rest of the iters
      if(model.res$Score > params$best.score & is.list(model.res)){
        print(paste0("round",i))
        print(paste0("score:", model.res$Score))
        
        params$best.score <- model.res$Score
        params$model.params <- list(K = K, alpha = alpha,
                                    iter = iter, C = C)
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
  
  ## next, construct similarity graphs
  graph.views <- lapply(data.views.dist, FUN = affinityMatrix, K, alpha)#This One
  
  #Run SNF
  graph.views$W_fused <- SNF(Wall = graph.views, K = K, t = iter)#This One
  
  #Cluster
  groupings = try(lapply(graph.views, FUN = function(graph.view){spectralClustering(graph.view, K = C, type = 2)}))#This one
  
  #Cast the previously calculated distance matrics to type dist
  data.views.dist <- lapply(data.views.dist, FUN = function(view.dist){as.dist(view.dist)})
  
  #Calculate silhouette scores...
  #Matching up single data views with respectively defined clusterings
  print(names(data.views))
  sils.single.views <- try(lapply(names(data.views),
                              FUN = function(idx){silhouette(x = groupings[[idx]],
                                                             dist = data.views.dist[[idx]])}))
  names(sils.single.views) <- names(data.views)
  
  sil.multiview <- try(lapply(data.views.dist,
                          FUN = function(data.view.dist){silhouette(x = groupings$W_fused,
                                                                    dist = data.view.dist)}))
  
  #Cast the previously calculated distance matrics to type dist
  data.views.dist <- lapply(data.views.dist, FUN = function(view.dist){as.dist(view.dist)})
  
  #Calculate silhouette scores...
  #Matching up single data views with respectively defined clusterings
  sils.single.views <- try(lapply(names(data.views),
                              FUN = function(idx){silhouette(x = groupings[[idx]],
                                                             dist = data.views.dist[[idx]])}))
  names(sils.single.views) <- names(data.views)
  
  sil.multiview <- try(lapply(data.views.dist,
                          FUN = function(data.view.dist){silhouette(x = groupings$W_fused,
                                                                    dist = data.view.dist)}))
  avg.sils <- lapply(sil.multiview, FUN = function(sil.score){mean(sil.score[,"sil_width"])})
  
  model.results$Score <- mean(unlist(avg.sils))
  model.results$Pred <- 0#groupings$W_fused
  
  return(model.results)
}

maximizeSilhouetteSNF <- function(data.views, K, alpha, iter, C, num.pcs = 100){
  #require(caret)
  #require(rpart)
  
  model.results <- list()
  
  #PCA
  data.views.pca <- lapply(data.views,
                           FUN = function(data.view){t(prcomp(x = data.view, rank. = num.pcs)$rotation)})
  
  #Calculate distance
  data.views.dist <- lapply(data.views.pca, FUN = function(view){dist2(as.matrix(view), as.matrix(view))})
  
  ## next, construct similarity graphs
  graph.views <- lapply(data.views.dist, FUN = affinityMatrix, K, alpha)#This One
  
  #Run SNF
  graph.views$W_fused <- SNF(Wall = graph.views, K = K, t = iter)#This One
  
  #Cluster
  groupings = try(lapply(graph.views,
                         FUN = function(graph.view){spectralClustering(affinity = graph.view,
                                                                       K = C, type = 1)}))#This one
  
  #Cast the previously calculated distance matrics to type dist
  data.views.dist <- lapply(data.views.dist, FUN = function(view.dist){as.dist(view.dist)})
  
  #Calculate silhouette scores...
  #Matching up single data views with respectively defined clusterings
  sils.single.views <-try(lapply(names(data.views),
                              FUN = function(idx){silhouette(x = groupings[[idx]],
                                                             dist = data.views.dist[[idx]])}))
  names(sils.single.views) <- names(data.views)
  
  sil.multiview <- try(lapply(data.views.dist,
                          FUN = function(data.view.dist){silhouette(x = groupings$W_fused,
                                                                    dist = data.view.dist)}))
  
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
  
  print(model.results$Score)
  
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
