#Author: Ang Andy Tu

#Function related to dealing with CDS objects and converting them to Seurat, and also running Monocle 3 on the data
#to get pseudotime trajectories.


#this takes a seurat object, and converts it to CDS object, then do a set of pre-processing to get cell cluters
cdsprocess <- function(seurat, res = 1e-4) {
  cds = importCDS(seurat)
  cds = updateCDS(cds)
  cds = estimateSizeFactors(cds)
  cds = estimateDispersions(cds)
  cds = preprocessCDS(cds, num_dim = 20, method = "PCA", pseudo_expr = 1.0, verbose = TRUE, fastpath = TRUE)
  cds <- reduceDimension(cds, reduction_method = 'UMAP', python_home = 'C:\\')
  cds = clusterCells(cds,method = 'louvain', louvain_iter = 5,res, verbose = TRUE, python_home ='C:\\')
  plot_cell_clusters(cds, color_by = 'Cluster')
  return(cds)
}

learnTrajectory <-  function(monocds, RGE_method = 'DDRTree'){
  monocds <- partitionCells(monocds)
  monocds <- learnGraph(monocds,  RGE_method = RGE_method)
  monocds
}

subprocess <- function(seurat, perp = 100, cutoff = 20) {
  seurat = FindVariableGenes(seurat, do.plot = FALSE)
  seurat = ScaleData(seurat, vars.to.regress = c('nUMI', 'percent.mito'), 
                     model.use = 'poisson', genes.use = seurat@var.genes)
  seurat = RunPCA(seurat, do.print = FALSE)
  seurat = RunTSNE(seurat, dims.use = 1:cutoff, perplexity = perp)
  seurat = SetAllIdent(seurat, 'expt')
  TSNEPlot(seurat)
  seurat = RunUMAP(seurat, dims.use = 1:cutoff)
  DimPlot(seurat, 'umap')
  seurat
}
