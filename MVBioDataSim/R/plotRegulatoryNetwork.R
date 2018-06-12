#' @title Plot regulatory networks
#' @description Plots a regulatory network
#' @param regnet A regulatory network
#' @param geneColor the color used to plot genes
#' @param mirnaColor the color used to plot mirnas
#' @param nodeSize the size of the nodes
#' @param arrowSize the size of the edges arrows
#' @param showLabel logical, should node labels be displayed?
#' 
plotRegulatoryNetwork <- function(regnet, geneColor="#b0b0b0", 
                                    mirnaColor="SkyBlue2", nodeSize=5, 
                                    arrowSize=0.2, showLabel=F)
{
    labels <- NA
    nodeColors <- c(rep(geneColor, regnet$genes), 
                    rep(mirnaColor, regnet$mirnas))
    
    if(showLabel) labels <- rownames(regnet$interactions)
    
    plot(graph.adjacency(abs(regnet$interactions)), vertex.size=nodeSize, 
         vertex.label=labels, edge.arrow.size=arrowSize, 
         vertex.color=nodeColors)
}