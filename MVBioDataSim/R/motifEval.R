.evaluateMotif<- function(motif, net, net.genes, net.degree.in, net.degree.out,
                          net.clustering.coefficient,
                          current.difference.in, current.difference.out, 
                          target.distribution.in, target.distribution.out, 
                          target.clustering.coefficient)
{
	motif.size <- motif$size
	motif.out.degree <- apply(abs(motif$motif), 1, sum)
	motif.clustering.coefficient <- .CC(motif$motif,T)[[1]]
  motif.genes <- (1:motif.size)
	    
    S.out <- matrix(0, nrow(net), motif.size)
    for(i in net.genes)
    {
        k.prime <- net.degree.out[i]
        for(j in motif.genes)
        {
          k.second <- k.prime + motif.out.degree[j]
          if(k.prime > 0) S.out[i,j] <- S.out[i,j] + sign(abs(current.difference.out[k.prime]) - abs(current.difference.out[k.prime] + 1)) * abs(current.difference.out[k.prime]) / target.distribution.out[k.prime]
          if(k.second > 0) S.out[i,j] <- S.out[i,j] + sign(abs(current.difference.out[k.second]) - abs(current.difference.out[k.second] - 1)) * abs(current.difference.out[k.second]) / target.distribution.out[k.second]
        }
    }
	
	if(!is.null(target.distribution.in))
	{
	    motif.in.degree <- apply(abs(motif$motif), 2, sum)
	    S.in <- matrix(0, nrow(net), motif.size)
	    for(i in net.genes)
	    {
        k.prime <- net.degree.in[i]
        for(j in motif.genes)
        {
            k.second <- k.prime + motif.in.degree[j]
            if(k.prime > 0) S.in[i,j] <- S.in[i,j] + sign(abs(current.difference.in[k.prime]) - abs(current.difference.in[k.prime] + 1)) * abs(current.difference.in[k.prime]) / target.distribution.in[k.prime]
            if(k.second > 0)S.in[i,j] <- S.in[i,j] + sign(abs(current.difference.in[k.second]) - abs(current.difference.in[k.second] - 1)) * abs(current.difference.in[k.second]) / target.distribution.in[k.second]
        }
	    }
        S.out <- (S.out + S.in) / 2
	}
    
    if(target.clustering.coefficient > net.clustering.coefficient)
    {
        CC <- motif.clustering.coefficient - net.clustering.coefficient
    }else{
        CC <- net.clustering.coefficient - motif.clustering.coefficient
    }
    
    return(list(score=S.out, clustering.coefficient=CC))	
}