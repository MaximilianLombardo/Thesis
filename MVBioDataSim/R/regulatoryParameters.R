.regulatoryParameters <- function(regnet)
{
  g <- Genes(regnet)
  m <- Mirnas(regnet)
  gCount <- regnet$genes
  mCount <- regnet$mirnas
  #gidx <- getIndex(regnet,g)
  #midx <- getIndex(regnet,m)
  
  parameters  <- list()
  degradation <- rnorm(gCount + mCount, 0.5, 0.1)
  degradation[degradation<0] <- 0
  degradation[degradation>1] <- 1

  production <- rnorm(gCount + mCount, 0.5, 0.1)
  production[production<0] <- 0
  production[production>1] <- 1

  
  names(degradation)  <- c(g,m)
  names(production)   <- c(g,m)
  
  h       <- Matrix(0, nrow=gCount + mCount, ncol=gCount + mCount, sparse=T)
  theta   <- Matrix(0, nrow=gCount + mCount, ncol=gCount + mCount, sparse=T)
  rownames(h) <- rownames(theta) <- rownames(regnet$interactions)
  colnames(h) <- colnames(theta) <- colnames(regnet$interactions)
  
  #The activation function maps the positive real line into the [0,1] interval.
  #But the expression values used are relative, so also input si bound to [0,1].
  #We assumed that when expression value is 1, then the activation must be at 
  #least 0.95. For a fixed threshold value theta, the minimum value for h such 
  #that activation is >= 0.95 when x=1 is h >= ln(0.05/0.95)/ln(theta)
  idx         <- regnet$interactions != 0
  count       <- sum(idx)
  theta[idx]  <- runif(count, 0, 1)
  theta[theta < 0] <- 0.06      #min value for which h > 1
  theta[theta > 1] <- 0.94
  lower       <- -2.944439/log(theta[idx])   #-2.944439 = log(0.05/0.95)
  h[idx]      <- sapply(lower, FUN=function(l){runif(1, l, 2*l)})
  
  
  gamma <- Matrix(0, nrow=mCount, ncol=gCount, sparse=T)

  rownames(gamma) <- rownames(regnet$interactions[m, g])
  colnames(gamma) <- colnames(regnet$interactions[m, g])
  idx <- regnet$interactions[m, g] != 0
  gamma[idx] <- runif(sum(idx), 0, 1)
  gamma[gamma < 0] <- 0
  gamma[gamma > 1] <- 1
  params <- list(production=round(production, 4), 
                 degradation=round(degradation,4), h=round(h,4),
                 theta=round(theta,4), gamma=round(gamma,4))
  return(params)
}