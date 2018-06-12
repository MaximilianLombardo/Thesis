#NOTE: to run this example the package Boruta must be installed
#install.packages('Boruta', dependencies=TRUE)
library(Boruta)

library(igraph)
load('feature_importance.RData')
load('../1000g100m.grn')

g=as.data.frame(t(Data$genes))
m=as.data.frame(t(Data$mirnas))
genes = paste('g', 1:1000, sep='')
mirnas = paste('m', 1:100, sep='')

s=as.character(Data$signallingGenes$x)

grn = graph.adjacency(abs(rete$adjacency), mode = 'directed')
pathLen = shortest.paths(grn, v = V(grn),to = V(grn))

gRelevance=Boruta(x = g,y = Data$subjLabel,doTrace = 2,maxRuns = 1000)
gRelMinPathLen = apply(pathLen[s, genes[gRelevance$finalDecision=='Confirmed']], 2, min)
gTab=table(gRelMinPathLen)

mRelevance=Boruta(x = m,y = Data$subjLabel,doTrace = 2,maxRuns = 1000)
mRelMinPathLen = apply(pathLen[s, mirnas[mRelevance$finalDecision=='Confirmed']], 2, min)
mTab=table(mRelMinPathLen)

totRelevance=Boruta(x=cbind(g,m), y = Data$subjLabel,doTrace = 2,maxRuns = 1000)
totRelMinPathLen = apply(pathLen[s, totRelevance$finalDecision=='Confirmed'], 2, min)
totTab=table(totRelMinPathLen)
