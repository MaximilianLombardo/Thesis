.regulatoryInteractions <- function(regnet, AND.p = 0.5)
{
  rules <- list()
  ops <- c("AND", "OR")
  ngenes <- regnet$genes
  nsignals <- regnet$control_genes
  for(g in rownames(regnet$interactions))
  {
    idx <- regulatingGenes(regnet, g)
    #g.idx <- getIndex(regnet, g)
    rg <- sapply(idx, function(tmp){
            #y <- call('[', quote(y), tmp)
            #yh <- call('^', y, call('[', quote(h), tmp, g.idx))
            #th <- call('^', call('[', quote(theta), tmp, g.idx), 
            #           call('[', quote(h), tmp, g.idx))
            y <- substitute(y[tmp], list(tmp=tmp))
            
#             hN <- call('+', regnet$parameters$h[tmp, g], 
#                        call('[', quote(hBN), tmp, g))
#             tN <- call('+', regnet$parameters$theta[tmp, g], 
#                        call('[', quote(tBN), tmp, g))
            
            yh <- call('^', y, call('[', quote(hN), tmp, g))
            
            #th <-regnet$parameters$theta[tmp, g]^regnet$parameters$h[tmp,g]
            th <- call('^', call('[', quote(tN), tmp, g), 
                       call('[', quote(hN), tmp, g))
            res <- call('/', yh, call('+', yh, th))
            if(regnet$interactions[tmp, g] < 0)
              res <- call('-', 1, res)
            return(res)
        })
    while(length(rg) > 1)
    {
        idxG <- .mySample(1:length(rg), 2, replace=F)
        op <- .mySample(ops, 1, prob=c(AND.p, 1-AND.p))
        expr <- switch(op, AND=call('min', rg[[idxG[1]]], rg[[idxG[2]]]),
                       OR=call('min', 1, call('+', rg[[idxG[1]]], rg[[idxG[2]]])))
        rg[idxG] <- NULL
        rg <- c(rg, expr)
    }
    activation <- if(length(rg) > 0) rg[[1]] else quote(0)
    idx <- regulatingMirnas(regnet, g)
    rmir <- sapply(idx, function(tmp){
        #yh <- call('^', y, call('[', quote(h), tmp, g.idx))
#         th <- call('^', call('[', quote(theta), tmp, g.idx), 
#                    call('[', quote(h), tmp, g.idx))
#         call('*', call('[', quote(gamma), tmp-ngenes, g.idx), 
#              call('/', yh, call('+', yh, th)))
        y <- call('[', quote(y), tmp)
#         
#         hN <- call('+', regnet$parameters$h[tmp, g], 
#                    call('[', quote(hBN), tmp, g))
#         tN <- call('+', regnet$parameters$theta[tmp, g], 
#                    call('[', quote(tBN), tmp, g))
#         gN <- call('+', regnet$parameters$gamma[tmp, g], 
#                    call('[', quote(gBN), tmp, g))
          
        yh <- call('^', y, call('[', quote(hN), tmp, g))
        
        th <- call('^', call('[', quote(tN), tmp, g),
                   call('[', quote(hN), tmp, g))
        #th <- regnet$parameters$theta[tmp,g]^regnet$parameters$h[tmp, g]
        res <- call('*', call('[', quote(gN), tmp, g), 
                    call('/', yh, call('+', yh, th)))
        return(res)
    })
 
    cs <- controllingSignals(regnet, g)
    if(length(cs) > 0)
    {
      s <- call('+', activation, call(cs, quote(t)))
      activation <- call('min', 1, s)
    }
    repression <- if(length(rmir) > 0) rmir[[1]] else quote(0)
    if(length(rmir) > 1)
    {
        for(i in 2:length(rmir))
        {
            repression <- call('+', repression, rmir[[i]])
        }
        repression <- call('min', 1, repression)
    }
#     pN <- call('+', regnet$parameters$production[g], call('[', quote(pBN), g))
#     dN <- call('+', regnet$parameters$degradation[g], call('[', quote(dBN), g))
    f <- call('-', call('*', call('[', quote(pN), g), activation), 
              call('*', call('[', quote(y), g), 
                   call('+', call('[', quote(dN), g), repression)))
#     f <- call('-', call('*', call('[', quote(production), g.idx), 
#                               activation), 
#                     call('*', call('[', quote(y), g.idx), 
#                          call('+', call('[', quote(degradation), g.idx), repression)))
    rules[[g]] <- f
  }
  return(rules)
}