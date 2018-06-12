.generateMotif <- function(nodes.max, regulators.max, p.repression)
{
    type <- .mySample(1:3, 1)
    return(switch(type, .motifStar(nodes.max[1], p.repression), 
                        .motifMultiple(nodes.max[2], regulators.max, p.repression),
                        .motifFFL(nodes.max[3], p.repression),
                        .moitfFBL(nodes.max[4], p.repression)))
}

.motifStar <- function(nodes.max, p.repression)
{
  motif.size <- .mySample(2:nodes.max, 1)
  motif <- matrix(0, nrow=motif.size, ncol=motif.size)
  regulator <- .mySample(1:motif.size, size=1)
  regulated <- .mySeq(1, motif.size)[-regulator]
      
  motif[regulator, regulated] <- .mySample(x=c(-1,1), size=motif.size - 1, 
                                           replace=T, prob=c(p.repression, 
                                                             1 - p.repression))
  return(list(motif=motif, 
              regulators=(1:motif.size) %in% regulator,
              regulated=(1:motif.size) %in% regulated, 
              size=motif.size,
              type="SIM"))
}

.motifMultiple <- function(nodes.max, regulators.max, p.repression)
{
  if(nodes.max < 3) {return(.motifStar(nodes.max, p.repression))}

  nodes <- .mySample(3:nodes.max, 1)
  regulators.max <- min(regulators.max, nodes - 1)
  regulators.size <- .mySample(2:regulators.max, 1)
  regulators <- .mySample(1:nodes, regulators.size, replace=F)
  
  toRegulate <- .mySeq(1, nodes)[-regulators]
  toRegulate.size <- length(toRegulate)
  
  motif <- matrix(0, nodes, nodes)
  for(r in regulators)
  {
      size <- trunc(runif(1, 1, toRegulate.size + 1))
      regulated <- .mySample(toRegulate, size=size, replace=F)
      motif[r, regulated] <- .mySample(x=c(-1,1), size=size, replace=T, 
                                       prob=c(p.repression, 1 - p.repression))
  }   
  
  #if(sum(diag(abs(motif)))) stop("uno\n")
  
  idx <- which((apply(abs(motif), 1, sum) + apply(abs(motif), 2, sum)) == 0)
  if(length(idx) > 0)
  {
      for(i in idx)
      {
          r <- .mySample(regulators, size=1)
          motif[r,i] <- .mySample(x=c(-1,1), size=1, prob=c(p.repression, 
                                                            1 - p.repression))
      }
  }
  
  #if(sum(diag(abs(motif)))) stop("due\n")
  
  return(list(motif=motif, 
              regulators=(1:nodes) %in% regulators, 
              regulated=(1:nodes) %in% toRegulate, 
              size=nodes,
              type="DOR"))
}

.motifFFL <- function(nodes.max, p.repression)
{
  if(nodes.max < 3) {return(.motifStar(nodes.max, p.repression))}
  
  n.nodes <- .mySample(.mySeq(3, nodes.max), 1, replace=F)
  
  motif <- matrix(0, nrow=n.nodes, ncol=n.nodes)
  perm <- .mySample(1:n.nodes, n.nodes, replace=F)

  edges <- rbind(c(perm[1], perm[2]),
                 cbind(perm[2], perm[3:n.nodes]))
  motif[edges] <- 1
  
  repressor <- rbinom(1, 1, p.repression)
  motif[cbind(perm[1], perm[3:n.nodes])] <- if(repressor) -1 else 1
  #motif[perm[n.nodes - 1], perm[n.nodes]] <- if(repressor) -1 else 1
  
  return(list(motif=motif, 
              regulators=(1:n.nodes) %in% perm[1],
              regulated=(1:n.nodes) %in% perm[-1],
              size=n.nodes,
              type="FFL"))
}

.moitfFBL <- function(nodes.max, p.repression)
{
  if(nodes.max < 3) {return(.motifStar(nodes.max, p.repression))}
  
  n.nodes <- .mySample(.mySeq(3, nodes.max), 1, replace=F)
  
  motif <- matrix(0, nrow=n.nodes, ncol=n.nodes)
  loop_participants <- .mySample(1:n.nodes, 2, replace=F)
  
  edges <- rbind(c(loop_participants[1], loop_participants[2]),
                 c(loop_participants[2], loop_participants[1]))
  
  motif[edges] <- if(rbinom(1,1,p.repression)) -1 else 1
  
  edges <- rbind(cbind(loop_participants[1], (1:n.nodes)[-loop_participants]),
                 cbind(loop_participants[2], (1:n.nodes)[-loop_participants]))
  
  motif[edges] <- .mySample(c(-1, 1), 2*(n.nodes-2), replace=T, 
                            prob=c(p.repression, 1-p.repression))    
  return(list(motif=motif, 
              regulators=(1:n.nodes) %in% loop_participants,
              regulated=(1:n.nodes) %in% (1:n.nodes)[-loop_participants],
              size=n.nodes,
              type="FBL"))
}