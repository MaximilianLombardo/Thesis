.GMInteractions <- function(ggnet, mgnet, gamma,
                            p.repression)
{
  ngenes <- nrow(ggnet)
  nmirnas <- nrow(mgnet)
  mirnas <- 1:nmirnas
  gmnet <- Matrix(0, ngenes, nmirnas, sparse=T)

  
  out.deg <- apply(abs(ggnet), 1, sum)
  regulators <- which(out.deg > 0)
  
  size.dist <- (1:(nmirnas))^(-gamma)
  size.dist <- size.dist / sum(size.dist)
  
  for(r in regulators)
  {
    size <- .mySample(1:nmirnas,size=1, prob=size.dist)
    mirna.scores <- rep(1, nmirnas) + apply(abs(mgnet), 1, 
                                            FUN=function(m){score <- sum(ggnet[r,] & m); return(m[r] * score + score)})
    regulated <- .mySample(1:nmirnas, size=size,replace=F,prob=mirna.scores)
    gmnet[r, regulated] <- .mySample(c(-1,1), size=size, replace=T,
                                   prob=c(p.repression, 1-p.repression))
  }
  return(gmnet)
}