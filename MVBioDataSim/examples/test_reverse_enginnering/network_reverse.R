#NOTE: do not run this code! it takes a lot of time and works only on linux.
#This is only to illustrate the procedure used to produce the the tables in the
#paper.
\dontrun{
library(minet)
library(igraph)

load('dataset')
load('../1000g100m.grn')

rownames(dataset$log) <- rownames(rete$adjacency)
gexpr <- dataset$log[1:1000, 1:75]

adj <- abs(rete$adjacency)
gene.adj <- abs(rete$adjacency[1:1000, 1:1000])

regnet <- graph.adjacency(adj)
gene.regnet <- graph.adjacency(gene.adj)

edgelist <- get.edgelist(regnet)
gene.edgelist <- get.edgelist(gene.regnet)

#PANDA network reconstruction
dir.create('panda_experiments')
write.table(x=dataset$log, file='panda_experiments/expression.txt',quote=F, 
            row.names=T, col.names=F, sep='\t')
write.table(x=gexpr, file='panda_experiments/gene_expression.txt',quote=F, 
            row.names=T, col.names=F, sep='\t')

ratios <- c(0.1, 0.25, 0.5, 0.75, 1.0, 1.1, 1.25, 1.5, 1.75, 2.0)
names <- rownames(adj)
gene.names <- rownames(gene.adj)

for(r in ratios)
{
  dir.name <- paste('panda_experiments/', as.character(r), '/', sep='')
  dir.create(dir.name)
  #random sample a portion of actual edges
  if(r <= 1.0)
  {
    edges <- edgelist[sample(nrow(edgelist),floor(nrow(edgelist)*r),replace=F),]
    gene.edges <- gene.edgelist[sample(nrow(gene.edgelist),floor(nrow(gene.edgelist)*r),replace=F),]
  }else{  #add random incorrect edges
    edges <- matrix(0, floor(nrow(edgelist)*r), 2)
    for(i in 1:floor(nrow(edgelist)*r))
    {
      is.edge <- 0
      while(!is.edge)
      {
        nodes <- sample(nrow(adj), 2, replace=T)
        is.edge <- adj[nodes[1], nodes[2]]
      }
      edges[i, ]<- names[nodes]
    }
    edges <-rbind(edges, edgelist)
    
    gene.edges <- matrix(0, floor(nrow(gene.edgelist)*r), 2)
    for(i in 1:floor(nrow(gene.edgelist)*r))
    {
      is.edge <- 0
      while(!is.edge)
      {
        nodes <- sample(nrow(gene.adj), 2, replace=T)
        is.edge <- gene.adj[nodes[1], nodes[2]]
      }
      gene.edges[i, ]<- gene.names[nodes]
    }
    gene.edges <-rbind(gene.edges, gene.edgelist)
  }
  
  write.table(cbind(edges, 1.0),file=paste(dir.name, 'edgelist.txt', 
                                              sep=''), 
              quote=F, row.names=F, col.names=F, sep='\t')
  write.table(cbind(gene.edges, 1.0),file=paste(dir.name, 
                                                   'gene_edgelist.txt', 
                                                   sep=''), 
              quote=F, row.names=F, col.names=F, sep='\t')
  
  system(paste(getwd(), '/panda/panda ', 
              '-e ', getwd(), '/panda_experiments/expression.txt ', 
              '-m ', getwd(), '/', dir.name, 'edgelist.txt ', 
              '-o ', getwd(), '/', dir.name, 'full', sep=''))
  
  system(paste(getwd(), '/panda/panda ', 
              '-e ', getwd(), '/panda_experiments/gene_expression.txt ', 
              '-m ', getwd(), '/', dir.name, 'gene_edgelist.txt ', 
              '-o ', getwd(), '/', dir.name, 'gene', sep=''))
}

# #ARACNE network reconstruction
gexpr.mim <- build.mim(gexpr, estimator='pearson')
colnames(gexpr.mim) <- rownames(gexpr.mim) <- colnames(rete$adjacency)[1:1000]
gene.reconet <- graph.adjacency(aracne(gexpr.mim) > 0.05)
gexpr.edgelist <- get.edgelist(gene.reconet)

#for each reconstructed edge, test if it is an actual edge or a path in the original network
paths <- apply(gexpr.edgelist, 1, FUN=function(x) {shortest.paths(regnet, v=x[1], to=x[2],mode='out')})
table(paths)
mean(paths[paths!=Inf])
}
