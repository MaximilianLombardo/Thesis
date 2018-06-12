#n.genes = number of genes
#n.mirna = number of miRNAs
#gamma = parameter of the power-law distribution
#cf.cl = desired mean clustering coefficent of each node
#motif.size = max number of nodes per motif
#max.regulators = max number of regulators per motif
#motif.pool = size of the set of motifs to be generated per round
.GGInteractions <- function(n.genes = 200, target.distribution.in = NULL, 
                        target.distribution.out, 
                        target.clustering.coefficient = 0.1, 
                        motif.size.max=c(5,8,5,5), motif.regulators.max = 5, 
                        motif.pool.size=50, p.repression=0.4, hub_genes = 5, 
                        progress=NULL)
{
  #total number of interacting molecules
  N <- n.genes
  net <- Matrix(0, nrow=N, ncol=N,sparse=T)
  
  V <- 1:N
  H <- vector()
  
  motif.size.max.start <- motif.size.max
  motif.regulators.max.start <- motif.regulators.max
  total <- N
  
  while(length(V) + length(H) > 1)
  {
    if(length(V) < 2)
    {
        V <- union(V, H)
        H <- integer(0)
        motif.size.max <- motif.size.max.start
        motif.regulators.max <- motif.regulators.max.start
    }
    
    motif.size.max <- pmin(motif.size.max, length(V))
    motif.regulators.max <- min(motif.regulators.max, motif.size.max[2] - 1)
    
    net.degree.out <- apply(abs(net), 1, sum)
    net.degree.out.distribution <- tabulate(net.degree.out,nbins=total)
    net.degree.out.distribution <- if(sum(net.degree.out.distribution)>0) net.degree.out.distribution / sum(net.degree.out.distribution) else rep(0, n.genes)
    net.degree.in <- apply(abs(net), 2, sum)
    net.degree.in.distribution <- tabulate(net.degree.in,nbins=total) 
    net.degree.in.distribution <- if(sum(net.degree.in.distribution)>0) net.degree.in.distribution / sum(net.degree.in.distribution) else rep(0, n.genes)
    net.clustering.coefficient <- MVBioDataSim:::.CC(net,T)[[1]]
    
    current.difference.in <- if(!is.null(target.distribution.in)) target.distribution.in - net.degree.in.distribution else NULL
    current.difference.out <- target.distribution.out - net.degree.out.distribution   
    
    #GENERATION AND EVALUATION OF MODULES
    motif.pool <- lapply(1:motif.pool.size, function(mp){
        motif<- .generateMotif(nodes.max=motif.size.max, 
                               regulators.max=motif.regulators.max,
                               p.repression=p.repression)
        S <- .evaluateMotif(motif=motif, 
                            net=net, 
                            net.genes=V, 
                            net.degree.in=net.degree.in, 
                            net.degree.out=net.degree.out, 
                            net.clustering.coefficient=net.clustering.coefficient,
                            current.difference.in=current.difference.in, 
                            current.difference.out=current.difference.out, 
                            target.distribution.in=target.distribution.in, 
                            target.distribution.out=target.distribution.out, 
                            target.clustering.coefficient=target.clustering.coefficient)
        return(c(motif, S))
    })
    
    #SAMPLING OF MODULE
    scores <- t(sapply(motif.pool, function(m){return(c(sum(m$score),m$clustering.coefficient))}))
    
    scores[which(scores[,1] < 0)] <- 0
    p.scores <- p.cc <- rep(1/motif.pool.size, motif.pool.size)
    if(sum(scores[,1]) > 0) p.scores <- scores[,1] / sum(scores[,1])
    
    p.cc <- pmax(0, p.cc + scores[,2])
    if(sum(p.cc) > 0) p.cc <- p.cc / sum(p.cc)
    
    p <- p.scores * p.cc
    if(sum(p) <= 0) p <- NULL 
    
    idx <- .mySample(.mySeq(1, motif.pool.size), size=1,replace=F,prob=p)
    motif <- motif.pool[[idx]]
    
    #ADD MODULE TO CURRENT NETWORK
    mapping <- .mapMotif(motif, V)
    net[mapping, mapping] <- net[mapping, mapping] + motif$motif

    V <- setdiff(V, mapping[(1:motif$size)])
    H <- union(H, mapping[motif$regulators])
    
    if(!is.null(progress)) progress$up(total - length(V) - length(H))
  }
  return(net)
}