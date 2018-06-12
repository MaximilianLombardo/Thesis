#'  @title Generation of regulatory networks
#'  @description Generation of a transcriptional and post-trascriptional regulatory network
#'  @param genes number of genes
#'  @param mirnas number of mirnas
#'  @param genes.in_degree a numeric value that specifies the scale parameter of a power-law degree distribution for the in degree of genes, or a numeric vector specifying a custom degree distribution
#'  @param genes.out_degree a numeric value that specifies the scale parameter of a power-law degree distribution for the out degree of genes, or a numeric vector specifying a custom degree distribution
#'  @param mirnas.in_degree a numeric value that specifies the parameter lambda of an exponential degree distribution for the in degree of mirnas, or a numeric vector specifying a custom degree distribution
#'  @param mirnas.out_degree a numeric value that specifies the parameter lambda of an exponential degree distribution for the out degree of mirnas, or a numeric vector specifying a custom degree distribution
#'  @param clustering_coefficient target clustering coefficient of the network
#'  @param max_regulators the maximum number of regulators that a node can have
#'  @param max_size vector of integers that specifies the size of each motif, if not sure leave it as default
#'  @param pool_size size of the pool of generated motifs for each step of the procedure
#'  @param repression probability that regulator represses one of its regulated nodes
#'  @param control_genes number of signalling genes controlling the state of the network durign simulation
#'  @param cooperation probability that two or more regulating nodes cooperate
#'  @param progress the type of progress bar to be shown
#'  
RegulatoryNetwork <- function(genes=200, mirnas=100, genes.in_degree=2.2, 
                              genes.out_degree=2.2, mirnas.in_degree=2.2,
                              mirnas.out_degree=0.1, clustering_coefficient=0.2, 
															max_regulators=5, max_size=c(5,8,5,5), 
                              pool_size=30, repression=0.35, control_genes=5, 
                              cooperation=0.5, progress=c('none', 'text', 'gui'))
{
  precision <- 1e-6  #can be changed
  #input checking
  if(!is.numeric(genes))
      stop("genes must be integer!")
  if(genes <= 0)
      stop("genes must be positive!")
  
  if(!is.numeric(mirnas))
      stop("mirnas must be integer!")
  if(mirnas < 0)
      stop("mirnas must be non-negative!")
  
  if(length(genes.in_degree) == 1)
  {
    if(genes.in_degree <= 0)
        stop("genes.in_degree must be a positive integer between 2 and 3 or a distribution of probability!")
  } else if(abs(1 - sum(genes.in_degree)) > precision)
      stop("genes.in_degree must sum to 1!")
  
  if(length(genes.out_degree) == 1)
  {
    if(genes.out_degree <= 0)
      stop("genes.out_degree must be a positive integer between 2 and 3 or a distribution of probability!")
  } else if(abs(1 - sum(genes.out_degree)) > precision)
    stop("genes.out_degree must sum to 1!")
  
  if(length(mirnas.in_degree) == 1)
  {
    if(mirnas.in_degree <= 0)
      stop("mirnas.in_degree must be a positive integer between 2 and 3 or a distribution of probability!")
  } else if(abs(1 - sum(mirnas.in_degree)) > precision)
    stop("mirnas.in_degree must sum to 1!")
  
  if(length(mirnas.out_degree) == 1)
  {
    if(mirnas.out_degree <= 0)
      stop("mirnas.out_degree must be a positive integer between 2 and 3 or a distribution of probability!")
  } else if(abs(1 - sum(mirnas.out_degree)) > precision)
    stop("mirnas.out_degree must sum to 1!")
  
  if(clustering_coefficient < 0 || clustering_coefficient > 1)
    stop("clustering_coefficient must be a positive real between 0 and 1!")
  
  if(max_regulators <= 0)
    stop("max_regulators must be a positive integer!")
  
  if(!all(max_size > 0))
    stop("max_size must be a vector of 4 positive integers!")
  
  if(pool_size <= 0)
    stop("max_regulators must be a positive integer!")
  
  if(repression < 0 || repression > 1)
    stop('repression must be a positive real between 0 and 1!')
  
  if(control_genes <= 0)
    stop("control_genes must be a positive integer!")
  
  if(cooperation < 0 || cooperation > 1)
    stop('cooperation must be a positive real between 0 and 1!')
  
  progress <- match.arg(progress)
  #end input check
  
  dist.in <- (1:genes) ^ (-genes.in_degree)
  dist.in <- dist.in / sum(dist.in)
  
  dist.out <- (1:genes) ^ (-genes.out_degree)
  dist.out <- dist.out / sum(dist.out)
  
  total <- genes + mirnas + 3
  
  if(!requireNamespace('tcltk', quietly=T))
  {  
    warning('tcltk package needed for gui progress bar is not installed. Fallink back to text mode.')
    progress <- 'text'
  }
  
  pb <- switch(progress, none=NULL, text=txtProgressBar(0,total,0,style=3),
               gui=tcltk::tkProgressBar(title="Generating Regulatory Network...",
                                 label="Generating Regulatory Network...",min=0,
                                 max=total, initial=0))
  
  #generate transcriptional interactions
  ggnet <- .GGInteractions(n.genes = genes, target.distribution.in = dist.in, 
                          target.distribution.out=dist.out, 
                          target.clustering.coefficient = clustering_coefficient, 
                          motif.size.max=max_size, 
                          motif.regulators.max = max_regulators, 
                          motif.pool.size=pool_size, 
                          p.repression=repression, 
                          hub_genes = control_genes, progress=pb)
  
  mgnet <- gmnet <- mmnet <- mlabels <- NULL
  
  if(mirnas > 0)
  {
    mgnet <- .MGInteractions(ggnet, mirnas, mirnas.out_degree, progress=pb)
    gmnet <- .GMInteractions(ggnet, mgnet, genes.out_degree, repression)
    mmnet <- Matrix(0, mirnas,mirnas,sparse=T)
    mlabels <- paste('m', 1:mirnas, sep='')
  }
  
  gDegree <- apply(abs(ggnet), 1, sum)
  
  controlled <- .mySample(1:genes, size=control_genes, replace=F, gDegree)
  
  sgnet <- Matrix(0, nrow=control_genes, ncol=genes, sparse=T)
  
  #at the moment each signal influences only a gene
  for(i in 1:control_genes) {sgnet[i, controlled[i]] <- 1}  
  
  glabels <- paste('g', 1:genes, sep='')
  slabels <- paste('s', 1:control_genes, sep='')
  
  net <- rbind2(cbind2(ggnet, gmnet),cbind2(mgnet, mmnet))
  
  rownames(net) <- colnames(net) <- c(glabels, mlabels)
  rownames(sgnet) <- slabels
  colnames(sgnet) <- glabels
    
  regnet <- list()
  class(regnet) <- 'RegulatoryNetwork'
  #regnet$adjacency <- net
  regnet$interactions <- net
  regnet$sgInteractions <- sgnet
  regnet$signals <- control_genes
  regnet$genes <- genes
  regnet$mirnas <- mirnas
  regnet$genes.inDegree <- genes.in_degree
  regnet$genes.outDegree <- genes.out_degree
  regnet$mirnas.inDegree <- mirnas.in_degree
  regnet$mirnas.outDegree <- mirnas.out_degree
  regnet$CC <- clustering_coefficient
  regnet$maxRegulators <- max_regulators
  regnet$maxSize <- max_size
  regnet$poolSize <- pool_size
  regnet$repression <- repression
  
  if(!is.null(pb)) pb$up(total - 2)
  
  if(!is.null(pb))
  {
    if(!is.null(pb$label))
    {
      pb$label("Defining parameters...")
    }
  }
  params <- .regulatoryParameters(regnet)
  regnet$parameters <- params
  
  if(!is.null(pb))
  {
      if(!is.null(pb$label))
      {
          pb$label("Defining interactions...")
      }
  }
  rules <- .regulatoryInteractions(regnet, AND.p = cooperation)
  regnet$rules <- rules
  if(!is.null(pb)) pb$up(total - 1)
  
  if(!is.null(pb)) pb$up(total)
  
  if(!is.null(pb)) pb$kill()
  return(regnet)
}