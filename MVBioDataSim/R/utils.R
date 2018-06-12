.CC <- function(graph, directed=TRUE)
{
  nodes <- dim(graph)[1]
  edges <- which(graph!=0, arr.ind=T)
  
  C <- sapply(seq(1, nodes), FUN=function(n)
    {
      neighbors <- union(edges[edges[,1] == n, 2], edges[edges[,2] == n, 1])
      neighbors.count <- length(neighbors)
      if(neighbors.count < 2) return(0)
      edges.count <- sum(graph[neighbors,neighbors])
      if(!directed) edges.count <- 2 * edges.count
      return(edges.count / (neighbors.count * (neighbors.count - 1)))
    })
  return(list(mean(C), C))
}

.mySample <- function(x, size, replace=FALSE, prob=NULL)
{
    if(is.null(x)) return(NULL)
    if(!replace && size > length(x)) stop("Cannot take a sample larger than the population when 'replace = FALSE'")
    if(length(x) == 1) return(rep(x, size))
    return(sample(x, size, replace, prob))
}

.mySeq <- function(from=1, to=1)
{
    if(from <= to)
    {
        return(seq(from, to))
    }else{
        return(numeric(0))
    }
}

Mirnas <- function(regnet)
{
  return(grep('^m.*', rownames(regnet$interactions), value=T))
}

Genes <- function(regnet)
{
  return(grep('^g.*', rownames(regnet$interactions), value=T))
}

Signals <- function(regnet)
{
  return(rownames(regnet$sgInteractions))
}

#Vectorized
getIndex <- function(regnet, names)
{
  idx <- sapply(names, FUN=function(n){
#     nn <- switch(substr(n,1,1), g=rownames(regnet$interactions), 
#                     m=rownames(regnet$interactions), 
#                     s=rownames(regnet$sgInteractions))
    if(substr(n,1,1)=='s')
    {
      return(grep(paste('^',n,'$',sep=''), rownames(regnet$sgInteractions)))
    }
    else
    {
      return(grep(paste('^',n,'$',sep=''), rownames(regnet$interactions)))
    }
#     return(which(nn==n))
  })
  return(idx)
}

regulatingGenes <- function(regnet, node)
{
  return(grep('^g.*', names(which(regnet$interactions[,node]!=0)),value=T))
}

#miRNAs can only regulate genes
regulatingMirnas <- function(regnet, gene)
{
  return(grep('^m.*', names(which(regnet$interactions[,gene]!=0)),value=T))
}

controllingSignals <- function(regnet, gene)
{
  if(substr(gene,1,1)!='g') return(character(0))
  return(names(which(regnet$sgInteractions[,gene]!=0)))
}

# removeFromNetwork <- function(regnet, name)
# {
#     idx <- getIndex(regnet, name)
#     new_net <- regnet
#     new_net$adjacency <- regnet$adjacency[-idx, -idx]
#     return(new_net)
# }