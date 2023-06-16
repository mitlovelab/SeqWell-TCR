#Author: Ang A. Tu
#Purpose: This is a common set of helper functions that are used throughout the paper to speed up analysis. Some custom changes 
#between datasets were made, those are specifically identified in each of the scripts.

# reccreate color hues of default ggplot colors
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

# plot clonal expansion on existing tsne map coordinates. assumes coordinates are 
# already in data_input
TR_expan_map <- function(data_input, TR_var = TRB_count, ft.size = 1,...){
  TR_var <- enquo(TR_var)
  
  TR_expan_plot <- ggplot(arrange(data_input,!!TR_var), aes(x = Tsne.x, y = Tsne.y)) + 
    geom_point(color = "grey87", size = ft.size,...) + geom_point(data = arrange(data_input,!!TR_var) %>% filter(!!TR_var>0), 
                                                              aes(color = log2(!!TR_var)), size = ft.size,...) +
    scale_color_viridis(direction = -1, na.value = "grey87")
  return(TR_expan_plot)
}

# add counts of TCRA and TCRB to seurat object. Also plots the expansion on the default plot output in rstudio.
Add_TRAB_count <- function(seurat){
  require(tidyverse)
  require(viridis)
  seurat@meta.data %>% dplyr::count(TRB_CDR3) %>% 
    filter(!is.na(TRB_CDR3)) %>% arrange(desc(n)) -> seurat_CDR3B_count
  colnames(seurat_CDR3B_count) <- c("TRB_CDR3","clones")
  seurat@meta.data %>% dplyr::count(TRA_CDR3) %>% 
    filter(!is.na(TRA_CDR3)) %>% arrange(desc(n)) -> seurat_CDR3A_count
  colnames(seurat_CDR3A_count) <- c("TRA_CDR3","clones")
  
  left_join(seurat@meta.data, seurat_CDR3B_count, by = "TRB_CDR3") %>% .$clones -> seurat@meta.data$TRB_count
  str(seurat@meta.data$TRB_count)
  
  ggplot(seurat@meta.data, aes(x = UMAP.x, y = UMAP.y, color = TRB_count)) + 
    geom_point() + scale_color_viridis(direction = -1, na.value = "grey87")
  
  left_join(seurat@meta.data, seurat_CDR3A_count, by = "TRA_CDR3") %>% .$clones -> seurat@meta.data$TRA_count
  str(seurat@meta.data$TRA_count)
  
  ggplot(seurat@meta.data, aes(x = UMAP.x, y = UMAP.y, color = TRA_count)) + 
    geom_point() + scale_color_viridis(direction = -1, na.value = "grey87")
  
  return(seurat)
}

# Simple wrapper that puts together a couple of the standard seurat processing steps. limit_vargenes is called if we
# want to limit the number of variable genes, usually to achieve better clustering without batch correction. Input is
# an already made seurat object.
Seurat_process <- function(seurat, perp = 30, cutoff = NULL, limit = FALSE,cf = 0.2, mf = 0.1,umap.run = TRUE,...) {
  seurat <- NormalizeData(object = seurat)
  if(limit!=FALSE){
    seurat <- limit_vargenes(seurat,limit,cf)
  }else{
    seurat <-  FindVariableGenes(seurat,x.low.cutoff = mf, y.cutoff = cf, do.plot = FALSE)
  }
  seurat <-  ScaleData(seurat, vars.to.regress = c('nUMI', 'percent.mito'), 
                       model.use = 'poisson', genes.use = seurat@var.genes)
  seurat <-  RunPCA(seurat, do.print = FALSE)
  if(is.null(cutoff)){
    cutoff <- find_pcCut(seurat)
  }
  print(paste("PCA cutoff used is",cutoff,"components"))
  seurat <-  RunTSNE(seurat, dims.use = 1:cutoff, perplexity = perp)
  TSNEPlot(seurat)

  if(umap.run){
    seurat <-  RunUMAP(seurat, dims.use = 1:cutoff)
    DimPlot(seurat, 'umap')
  }
  return(seurat)
}

# function to find the cutoff of the pc elbow plot, basically finding the specific number of ranked PCs that's furtherest
# away from the fitted line between max and mim PCs, thus consituting the "elbow" of the plot
find_pcCut <- function(seurat){
  #pull out the pca deviation from each PCA components
  PCA_dev <- seurat@dr$pca@sdev
  allCoor<-cbind(seq_along(PCA_dev),PCA_dev)
  line_Vec<-allCoor[length(PCA_dev),]-allCoor[1,]
  ## normalize the line vector
  lineVecN = line_Vec / sqrt(sum(line_Vec^2));
  ## find the distance from each point to the line: vector between all points and first point
  
  vecFromFirst<-allCoor-do.call("rbind", rep(list(allCoor[1,]), length(PCA_dev)))
  q<-do.call("rbind", rep(list(lineVecN), length(PCA_dev)))
  scalarProduct<-q[,1]
  for (i in 1:length(PCA_dev)) {
    scalarProduct[i]<- vecFromFirst[i,] %*% q[i,]
  }
  vecFromFirstParallel = scalarProduct * q
  vecToLine = vecFromFirst - vecFromFirstParallel
  distToLine = sqrt(rowSums(vecToLine^2))
  ##Point in elbow is point furthest from line connecting first and last point
  pcCut<-which.max(distToLine) 
  return(pcCut)
}

# function for limiting variable genes used, if desired
limit_vargenes <- function(seurat,n_var_genes = 1000, cf = 0.2){
  length_var <- length(seurat@var.genes)
  while (length_var>n_var_genes) {
    seurat <- FindVariableGenes(object = seurat, x.low.cutoff = 0.1, y.cutoff = cf, do.plot = FALSE)
    length_var <- length(seurat@var.genes)
    cf=cf+.05
  }
  return(seurat)
}

# Takes in meta data from seurat object, and returns a heatmap matrix that can be formatted to be plotted with pheatmap.
# clone cutoff determines which clones we want to exclude (e.g. plotting singletons don't make much sense)
TR_cluster_heatmap <- function(seurat_meta, TR_chain = TRB_CDR3, clone_cutoff = 10, cluster_var = res.1){
  TR_chain <- enquo(TR_chain)
  cluster_var <- enquo(cluster_var)
  
  TR_CDR3_count <- seurat_meta %>% count(!!TR_chain, sort = TRUE) %>% filter(!is.na(!!TR_chain))
  TR_list_hmp <- TR_CDR3_count %>% filter(n>=clone_cutoff) %>% select(!!TR_chain) %>% pull(!!TR_chain)
  
  Meta_data_hmp <- seurat_meta %>% filter(!!TR_chain %in% TR_list_hmp) %>% group_by(!!TR_chain) %>% 
    count(!!cluster_var) %>% mutate(normalized = n/sum(n)) %>% complete(!!cluster_var)
  
  Meta_data_hmp$n[is.na(Meta_data_hmp$n)] <- 0
  Meta_data_hmp$normalized[is.na(Meta_data_hmp$normalized)] <- 0
  return(Meta_data_hmp)
}

# the previous funciton returns the format in long, here, we'll put it in wide, so we can return heatmap-ready matrix.
TR_heatmap_spread <- function(Meta_heatmap, TR_chain = TRB_CDR3, cluster_var = res.1){
  TR_chain <- enquo(TR_chain)
  cluster_var <- enquo(cluster_var)
  Meta_heatmap_w <- Meta_heatmap %>% select(!!TR_chain, !!cluster_var, normalized) %>% spread(!!cluster_var,normalized)
  hmp_rownames <- Meta_heatmap_w %>% pull(!!TR_chain)
  
  Meta_heatmap_w <- Meta_heatmap_w  %>% ungroup() %>% select(-!!TR_chain) %>% as.matrix()
  rownames(Meta_heatmap_w) <- hmp_rownames
  
  return(Meta_heatmap_w)
}

# a function to make the heatmap plot. Directly plots it out in the current plot device, or save to png.
TR_heatmap_plot <- function(Meta_heatmap, x_axis ="cluster", y_axis = "TRB clones",name = "TRB_heatmap.png", color_levels = viridis(10),row_clust = TRUE){
  require(grid)
  png(filename = name ,width = 600, height = 450)
  setHook("grid.newpage", function() pushViewport(viewport(x=1,y=1,width=0.9, height=0.9, name="vp", just=c("right","top"))), action="prepend")
  pheatmap(Meta_heatmap, color = color_levels, border_color =NA, drop_levels = TRUE, fontsize = 8, cluster_rows = row_clust, cluster_cols = FALSE)
  setHook("grid.newpage", NULL, "replace")
  grid.text(x_axis, y=-0.07, x = 0.35, gp=gpar(fontsize=16))
  grid.text(y_axis, x=-0.07, rot=90, gp=gpar(fontsize=16))
  dev.off()
}

#The following is to count only when you need to institute a cutoff, mostly when you want
#to set a cutoff on what's considered expanded clones.
Count_TRB <- function(Meta_data, cutoff = 3, group_var = orig.ident){
  group_var <- enquo(group_var)
  count_tb <- Meta_data %>% group_by(!!group_var) %>% filter(!is.na(TRB_CDR3)) %>% dplyr::count(TRB_CDR3,TRBV,TRBJ, sort = TRUE)
  colnames(count_tb) <- c("sampleID", "TRB_CDR3","TRBV","TRBJ","TRB_count")
  return(count_tb %>% ungroup %>% filter(TRB_count>= cutoff))
}

Count_TRA <- function(Meta_data, cutoff = 3){
  count_tb <- Meta_data %>% group_by(orig.ident) %>% filter(!is.na(TRA_CDR3)) %>% dplyr::count(TRA_CDR3,TRAV,TRAJ, sort = TRUE)
  colnames(count_tb) <- c("sampleID", "TRA_CDR3","TRAV","TRAJ","TRA_count")
  return(count_tb %>% ungroup %>% filter(TRA_count>= cutoff))
}

Count_TRAB <- function(Meta_data,cutoff = 1){
  count_tb <- Meta_data %>% group_by(orig.ident) %>% filter(!is.na(TRA_CDR3),!is.na(TRB_CDR3)) %>% 
    dplyr::count(TRA_CDR3,TRAV,TRAJ,TRB_CDR3,TRBV,TRBJ, sort = TRUE)
  colnames(count_tb) <- c("sampleID", "TRA_CDR3","TRAV","TRAJ","TRB_CDR3","TRBV","TRBJ","TRAB_count")
  return(count_tb %>% ungroup %>% filter(TRAB_count>= cutoff))
}

# output TRAB tables with counts, using the helper functions above.
Make_TRAB_table <- function(Meta_data,cutoff = 1){
  TRB_temp <- Count_TRB(Meta_data,cutoff)
  TRA_temp <- Count_TRA(Meta_data,cutoff)
  TRAB_temp <- Count_TRAB(Meta_data,cutoff)
  
  output <- left_join(TRAB_temp,TRB_temp,by = c("sampleID","TRB_CDR3","TRBV","TRBJ")) %>%
    left_join(TRA_temp,by = c("sampleID","TRA_CDR3","TRAV","TRAJ"))
  
  return(output)
}

# A wrapper around the Findmarker fnction in Seurat. helps to generate DeG tables that can then be used to generate
# Volcano plots.
# Functions for generating volcano plots quickly and easily, along side the proper data wrangling

DE_test <- function(seurat, group_1, group_2 = NULL, method = "bimod",max_p = 10e-30, min_p = 0.001, 
                    foldchangethresh = 2,logfc.threshold = 0,max_cells = Inf, ...){
  require(tidyverse)
  require(Seurat)
  
  DE_output <- FindMarkers(seurat,logfc.threshold = logfc.threshold, ident.1 = group_1, ident.2 = group_2,
                           test.use = method, do.print = TRUE, max.cells.per.ident = max_cells, ...)
  DE_output <- as_tibble(DE_output, rownames = "genes")
  DE_output <- DE_output %>% mutate(Graph_adj_p = ifelse(p_val_adj > max_p, p_val_adj, max_p),
                                    significance = ifelse(p_val_adj< min_p, "p-significant","not significant"))
  DE_output <- DE_output %>% mutate(significance = ifelse((abs(avg_logFC)> log(foldchangethresh)) & significance == "p-significant", "significant",significance)) %>%
    mutate(significance = factor(significance, levels = c("not significant","p-significant","significant")))
  
  return(DE_output)
}

# Asuuming we ran the DE_test function above. Takes the results and make a formated volcano plot.
# modified so that if no list of genes are supplied, automatically look for top differnetially expressing genes in terms
# of increasing and decrreasing fold-change genes. also include the nudge parameter.
plot_volcano <- function(DE_results,gene_list_1=NULL, gene_list_2 = NULL, max_p = 10e-30, 
                         min_p = 0.001, right_nudge = 0.8, left_nudge = 0){
  require(ggrepel)
  require(ggplot2)

  volcano <- ggplot(DE_results, aes(avg_logFC,-log10(Graph_adj_p))) + 
    geom_point(aes(color = significance), size = 0.3) + geom_hline(yintercept = -log10(min_p), linetype = 2) +
  geom_hline(yintercept = -log10(max_p),linetype = 2) + scale_color_manual(values =c("grey87","magenta3","darkturquoise")) +
  labs(x = "ln(fold change)", y = "-log10(adjusterd p-value)")
  
  if(!is.null(gene_list_1)){
    right_labeled_genes <- DE_results %>% filter(genes %in% gene_list_1)
    volcano <- volcano + geom_text_repel(data = right_labeled_genes, aes(label = genes),
                                         segment.size = 0.2, segment.color = "grey50", size =2, nudge_y = -2,
                                         nudge_x = right_nudge - right_labeled_genes$avg_logFC, direction = "y", hjust = 1) 
  } else{
    right_labeled_genes <- DE_results %>% arrange(desc(p_val_adj),desc(avg_logFC)) %>% top_n(10,avg_logFC)
    volcano <- volcano + geom_text_repel(data = right_labeled_genes, aes(label = genes),
                                         segment.size = 0.2, segment.color = "grey50", size =2, nudge_y = -2,
                                         nudge_x = right_nudge - right_labeled_genes$avg_logFC, direction = "y", hjust = 1) 
  }
  if(!is.null(gene_list_2)){
    left_labeled_genes <- DE_results %>% filter(genes %in% gene_list_2)
    volcano <- volcano + geom_text_repel(data = left_labeled_genes, aes(label = genes),
                                         segment.size = 0.2, segment.color = "grey50", size =2, nudge_y = -2,
                                         nudge_x = left_nudge + left_labeled_genes$avg_logFC, direction = "y", hjust = 1)
  } else{
    left_labeled_genes <- DE_results %>% arrange(desc(p_val_adj),(avg_logFC)) %>% top_n(10,-avg_logFC)
    volcano <- volcano + geom_text_repel(data = left_labeled_genes, aes(label = genes),
                                         segment.size = 0.2, segment.color = "grey50", size =2, nudge_y = -2,
                                         nudge_x = left_nudge + left_labeled_genes$avg_logFC, direction = "y", hjust = 1)  
  }
  return(volcano)
}

# A function to save pheatmap plots
save_pheatmap <- function(x, filename, device = "png", unit = "in",width=3, height=3) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  if(device =="png"){
    png(filename, width = width, height = height, units = unit, res = 1200)
  } else if(device == "pdf"){
    pdf(filename, width = width, height = height, family = "ArialMT")
  }
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

# Use ggseqlogo to construct loco plots of TCR sequences of a given length
Sequence_logo <- function(Data, TR_column, TR_length, length, method_type = "bits"){
  require(ggseqlogo)
  TR_column <- enquo(TR_column)
  TR_length <- enquo(TR_length)
  Data %>% filter((!!TR_length) == length) %>% select(!!TR_column) %>% pull(!!TR_column)
  logo <- ggseqlogo(Data %>% filter(!!TR_length == length) %>% select(!!TR_column) %>% pull(!!TR_column), method = method_type)
  return(logo)
}

# Old fn to save pheatmap into pdf. saved for reference. 
save_pheatmap_pdf <- function(x, filename, width=8, height=2) {
  pdf(filename, width = width, height = height, fonts = "ArialMT")
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

# set of two functions that make the clonality graphs.
# First is a function to count TCR chain by some cluster identity
Count_clones <- function(Data, group = Clusters,TR_chain = TRB_CDR3){
  #This just gerenates a dataframe of sample ID, TRB CDR3, and frequency
  group <- enquo(group)
  TR_chain <- enquo(TR_chain)
  Data %>% group_by(!!group) %>% filter(!is.na(!!TR_chain)) %>%dplyr::count(TRB_CDR3,sort = TRUE) -> count_df
  #colnames(count_df) <- c("orig.ident", "TRB_CDR3","frequency")
  return(as_tibble(count_df))
}

# takes the count, and calculates the proprotion of all clones in each cluster identity that fall into a particular
# exapnsion interval, starts with singletons at the bottom, and go by 10 equal bins.
Clonal_proportion <- function(x,bin.size = 10,group = orig.ident){
  
  group <- enquo(group)
  #x is a  sample name
  count_df<-Count_clones()
  #then I think I want to pull out the count for each of the samples.then i'm not sure.
  count_df_x <- count_df %>% filter(!!group %in% x)
  total_rows <- count_df_x %>% nrow()
  #set intervals to group clones
  singletons_rows <- count_df_x %>% filter(frequency==1) %>% nrow()
  if(round(total_rows-singletons_rows,-1) == 0){
    seq.max = round(total_rows-singletons_rows)
  }
  else{
    seq.max = round(total_rows-singletons_rows, -1)
  }
  intervals <- c(seq(0,seq.max,bin.size),total_rows)
  #print(bin.size)
  text_labels <- intervals[-1]
  print(text_labels)
  if(bin.size > 1){
    text_labels[-1] <- paste(text_labels+1,lead(text_labels), sep = "-")
    text_labels[1] <- paste(1,text_labels[1],sep = "-")
  }
  #text_labels[1] <- paste("top",text_labels[1])
  cut(row(count_df_x[1]),intervals, labels = text_labels) -> interval_labels
  cbind(count_df_x,intervalL = interval_labels) -> count_df_x
  #We'll save our data to another data frame. Initiate the dataframe here
  output <- data.frame(interval = levels(interval_labels), placeholder = NA)
  names(output)[2] <- x
  #then we should have a data frame with columns of interval, sample name (x), and NA in each of the categories
  #how do we calculate the sum?
  clone_sum = sum(count_df_x$frequency)
  for(y in output$interval){
    count_df_x %>% filter(intervalL == y) %>% ungroup() %>% select(frequency) %>% sum() -> temp
    index <- match(y,output$interval)
    output[index,2]<- temp/clone_sum
  }
  return(output)
}

# helper functions to move txne and umap coordinates to meta data
TSNE_xy <- function(seurat) {
  seurat@meta.data$Tsne.x <- seurat@dr$tsne@cell.embeddings[,1]
  seurat@meta.data$Tsne.y <- seurat@dr$tsne@cell.embeddings[,2]
  return(seurat)
}

UMAP_xy <- function(seurat) {
  seurat@meta.data$UMAP.x <- seurat@dr$umap@cell.embeddings[,1]
  seurat@meta.data$UMAP.y <- seurat@dr$umap@cell.embeddings[,2]
  return(seurat)
}

# make the dot plot per count in ggplot. It's an explicit version of geom_dot, which behaves oddly.
ggplot_dot <- function(x, count, ...) {
  require(ggplot2)
  message("The count variable must be an integer")
  count = as.integer(count) # make sure these are counts
  n = sum(count) # total number of dots to be drawn
  x = rep(x, count) # make x coordinates for dots
  count = count[count > 0]  # drop zero cases 
  y = integer(0)  # initialize y coordinates for dots
  for (i in seq_along(count)) 
    y <- c(y, 1:(count[i]))  # compute y coordinates
  ggplot(data.frame(x = x, y = y), aes(x = x, y = y)) +
    geom_point(aes(fill = log2(x)),...)  # draw one dot per positive count
}

# A helper function to get mean gene expression by TCR. This helps to reuturn a gene expression matrix that is 
# clonotype by genes. cells_id and TRB_list should be paired, so TRB_list has repeats. TRB_order is a clean list of 
# distinct clonotypes, in the order they should be on the heatmap rows.
clonotype_gene_data <- function(seurat, genes, cells_id, TRB_list, TRB_order, min = -2, max = 2){
  require(seurat)
  require(tidyverse)
  Heatmap_data <- FetchData(seurat, vars.all = genes, cells.use = cells_id)
  Heatmap_data <- as_tibble(Heatmap_data) %>% mutate(cell_id = cells_id, TRB_CDR3 = TRB_list) %>%
    select(cell_id, TRB_CDR3, everything()) %>% mutate(TRB_CDR3 = factor(TRB_CDR3,levels = TRB_order))

  Heatmap_data_collapse <- Heatmap_data %>% select(-cell_id) %>% group_by(TRB_CDR3) %>% summarise_all(funs(mean)) %>%
    mutate_if(is.double,scale)

  Heatmap_data_collapse_matrix <- as.matrix(Heatmap_data_collapse[-1]) #get rid of the TRB from the matrix

  rownames(Heatmap_data_collapse_matrix) <- Heatmap_data_collapse$TRB_CDR3

  Heatmap_data_collapse_matrix[which(Heatmap_data_collapse_matrix<min)] <- min
  Heatmap_data_collapse_matrix[which(Heatmap_data_collapse_matrix>max)] <-  max

  return(Heatmap_data_collapse_matrix)
}

# helper funciton to sort the leaves of a dendrogram.
sort_hclust <- function(...,method = "min") as.hclust(dendsort(as.dendrogram(...),type = method))

# A function to plot specific clonotype onto a background of grey points on a umap. THis is helpful in visualizing 
# where specific clones are on a umap/tsne.
plot_TR <- function(Meta, clonotype, TR_chain = TRB_CDR3, color = "red",background = "grey80",
                    pt_size = 0.5, bk_size = 0.5, ...){
  TR_chain <- enquo(TR_chain)
  TR_plot <- ggplot(Meta, aes(x = Tsne.x, y = Tsne.y)) + geom_point(size = bk_size, stroke = 0, pch = 16, color = background) +
    geom_point(data = Meta %>% filter(!!TR_chain %in% clonotype), aes(x = Tsne.x, y = Tsne.y), color = color, ...)
  return(TR_plot)
}
