###QUAKE glioblastoma single cell data set

library(SIMLR)
library(Seurat)

#load quake data
quake <- read.table("C:/Thesis/data/QuakeSingleCell/rawData/GBM_raw_gene_counts.csv", header = TRUE, row.names = 1)


#Get variable genes to reduce dimensionality
quake <- CreateSeuratObject(raw.data = quake, min.cells = 3, min.genes = 200, 
                           project = "GBM_QUAKE")

mito.genes <- grep(pattern = "^MT-", x = rownames(x = quake@data), value = TRUE)
percent.mito <- Matrix::colSums(quake@raw.data[mito.genes, ])/Matrix::colSums(quake@raw.data)

quake <- AddMetaData(object = quake, metadata = percent.mito, col.name = "percent.mito")


quake <- NormalizeData(object = quake, normalization.method = "LogNormalize", 
                      scale.factor = 10000)

quake <- FindVariableGenes(object = quake, mean.function = ExpMean, dispersion.function = LogVMR, 
                          x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)

#Perform SIMLR on sampled genes
quake.norm <- read.table("C:/Thesis/data/QuakeSingleCell/rawData/GBM_normalized_gene_counts.csv", header = TRUE, row.names = 1)

quake.norm <- quake.norm[quake@var.genes,sample(1:ncol(quake.norm), 500)]

#Estimate number of clusters
est.num.clusters <- SIMLR_Estimate_Number_of_Clusters(quake.norm, NUMC = 2:10)

#Perform SIMLR
simlr.result <- SIMLR(quake.norm, c = 6, normalize = FALSE)
