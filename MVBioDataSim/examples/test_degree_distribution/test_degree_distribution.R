#TEST#2: generate 50 networks with 1000 nodes to test degree distribution
#IT TAKES SOME TIME TO GENERATE THE DATASET
#TO REPRODUCE THE PLOTS THE LIBRARY GGPLOTS SHOULD BE INSTALLED, OTHERWISE FALLBACK PLOTS WILL BE PRODUCED
# nets <- 50
# nodes <- 1000
# dataset <- list()
# cat("Generating test#2\n")
# pb <- txtProgressBar(min=1,max=nets,initial=1,style=3)
# for(iter in 1:nets)
# {
#     setTxtProgressBar(pb, value=iter)
#     
# 		regnet <- RegulatoryNetwork(genes=nodes, mirnas=0, control_genes=5, progress='none')
# 		
#     i <- apply(abs(regnet$interactions), 2, sum)
#     o <- apply(abs(regnet$interactions), 1, sum)
#     
#     dataset[[iter]] <- list(adj = regnet$interactions, out_degree=o, in_degree=i, dist=tabulate(o,nodes)+tabulate(i,nodes))
# }
# close(pb)

has_ggplot = require(ggplot2)
load('test_degree_distribution.RData')
p <- unlist(lapply(dataset, '[[', 'dist'))
k <- rep(1:1000, 50)

df <- data.frame(k=k, p=p)
zero <- df[,1]==0|df[,2]==0

l <- lm(log(p)~log(k), data=df, subset=!zero)

summary(l)

if(has_ggplot)
{
  p=ggplot(df[!zero,], aes(x=k, y=p)) + scale_x_log10() + scale_y_log10() + geom_point() +
    geom_smooth(method='lm') + labs(title='Degree Distribution', x='log(k)', y='log(P(k))')
    print(p)
  #ggsave(filename = 'scalefree.pdf', plot = p, width = 10, height = 10, units = 'in',dpi = 100)
}else{
plot(log(k[!zero]), log(p[!zero]), ylab='log(P(k))',xlab='log(k)', main='Degree Distribution')
abline(l)
}
