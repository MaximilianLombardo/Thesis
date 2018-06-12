.mapMotif <- function(motif, genes)
{   
    motif$score[motif$score < 0] <- 0
    motif.genes <- (1:motif$size)
    mapping <- integer(length=motif$size)
    
    for(g in motif.genes)
    {
        s <- motif$score[genes, g]
        ss <- sum(s)
        if((ss <= 0) || is.infinite(s)) {p <- NULL}
        else{ p <- s/ss }
        
        if(length(genes) > 1) gene <- .mySample(1:length(genes), size=1, prob=p)
        else gene <- 1
        
        mapping[g] <- genes[gene]
        genes <- genes[-gene]
    }
    
    return(mapping)
}