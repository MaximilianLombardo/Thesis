#TEST#1: generate 100 networks with 10-1000 nodes to test clustering coefficient scale invariance
#IT TAKES SOME TIME TO GENERATE THE DATASET
#TO REPRODUCE THE PLOTS THE LIBRARY GGPLOTS SHOULD BE INSTALLED, OTHERWISE FALLBACK PLOTS WILL BE PRODUCED
# nets <- 100
# dataset <- list()
# cat("Generating test#1\n")
# pb <- txtProgressBar(min=1,max=nets,initial=1,style=3)
# for(iter in 1:nets)
# {
#     setTxtProgressBar(pb, value=iter)
#     nodes <- sample(10:1000, 1)
#     
# 		regnet <- RegulatoryNetwork(genes=nodes, mirnas=0, control_genes=5, progress='none')
# 		
#     o <- apply(abs(regnet$interactions), 1, sum)
#     i <- apply(abs(regnet$interactions), 2, sum)
#     cc <- MVBioDataSim::.CC(abs(regnet$interactions))
#     
#     dataset[[iter]] <- list(adj = regnet$interactions, out_degree=o, in_degree=i, cc=cc[[2]])
# }
# close(pb)

has_ggplot = require(ggplot2)
load('test_clustering_coefficient.RData')
fits <- lapply(dataset, function(d){df <- aggregate(d$cc, by=list(k=d$out_degree+d$in_degree), FUN=mean);
                                      zero <- df[,1]==0 | df[,2]==0;
                                      return(lm(log(x)~log(k), data=df,subset=!zero))})
size <- sapply(dataset, function(d){return(length(d$out_degree))})
coeffs <- sapply(fits, function(f){return(f$coefficients[2])})
df=data.frame(size=size, coeff=coeffs)

if(has_ggplot)
{
  p=ggplot(df, aes(x=size, y=coeff)) + geom_point() + ylim(c(-2,0)) +
  labs(title='Clustering coefficient scale invariance', x='Network size', y='Clustering coefficient distribution parameter')
  print(p)
  #ggsave(filename = 'cluscoef.pdf', plot = p, width = 10, height = 10, units = 'in',dpi = 100)
}else{
  plot(df, ylim=c(-2,0), ylab='Clustering coefficient distribution parameter', xlab='Network size', main='Clustering coefficient scale invariance')
  fit=lm(coeff~size, df)
  abline(fit)
}