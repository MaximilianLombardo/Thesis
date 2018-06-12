.MGInteractions <- function(ggnet, mirnas, lambda, progress=NULL)
{
    genes <- nrow(ggnet)
    total <- genes + mirnas
    
    if(!is.null(progress))
    {
        if(!is.null(progress$label))
        {
            progress$label("Adding miRNA mediation...")
        }
    }
    
    mgnet <- Matrix(0, nrow=mirnas, ncol=genes,sparse=T)
    
    genes.degree.out <- apply(abs(ggnet), 1, sum)
    genes.degree.in <- apply(abs(ggnet), 2, sum)
    flux <- genes.degree.out * genes.degree.in
    if(sum(flux)==0) flux=NULL
    flux[flux==0] <- 1/length(which(flux==0))
    
    for(i in 1:mirnas)
    {
        size <- ceiling(rexp(1, lambda))
        while(size > genes) {size <- ceiling(rexp(1, lambda))}
        s <- .mySample((1:genes), size=size, replace=F, prob=flux)
        mgnet[i, s] <- -1
        if(!is.null(progress)) progress$up(total - i)
    }
    return(mgnet)
}