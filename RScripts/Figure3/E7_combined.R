#Author: Ang Andy Tu
#Purpose: Analysis of the E7 immunized mice in Figure 3. combine ex vivo and non ex vivo stim
#Note: The analysis starts with some house keeping to combine stim and non-stimmed seurat objects. A combined object is 
# already provided, so that can be directly loaded as well.

# Load libraries and functions --------------------------------------------

library(tidyverse)
library(viridis)
library(Seurat)
library(pheatmap)
require(reticulate)
library(RColorBrewer)
library(dendsort)
library(ggbeeswarm)
library(ggpubr)

use_condaenv("r-reticulate-monocle")
# load in the script that contains the common scripts for the functions that we need for figures and analysis
function_file<-choose.files(caption="Choose scripts that contain common functions")
source(function_file)

graphing_parameters <- choose.files(caption="Choose scripts that contain graphing parameters")
source(graphing_parameters)
# functions ---------------------------------------------------------------
#Some of the functions will be more custom to this particular dataset becasue of the conditions that it has to filter for

get_cell_id <- function(seurat_meta,stimulation,clonotype,TR_column = TRB_CDR3){
  TR_column <- enquo(TR_column)
  seurat_meta %>% filter(Stimulation == stimulation) %>% filter(!!TR_column %in% clonotype) %>% 
    select(TRB_CDR3,cell_id)
}

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

Test_res <- function(seurat, res_test = seq(from = 0.4, to = 1.0 ,by = 0.2), max_cells = 200,
                     dims_use = 1:5){
  
  markers_list <- vector(list,length(res_test))
  markers_list <- map(res_test,function(res_test) 
    FindClusters(SubsetData(seurat,max.cells.per.ident = max_Cells), dims.use=dims_use, resolution=res_test, print.output = 0))
  names(markers_list) <- str_c('res.',res_test)
  markers_all <- FindAllMarkers(seurat,test.use = "roc", do.print = TRUE)
  markers_use <- markers_all %>% filter(avg_diff>0&power>0.3) %>% 
    mutate(cluster = factor(cluster)) %>% arrange(cluster, desc(avg_diff))
  return(DoHeatmap(SubsetData(Combined_seurat,max.cells.per.ident = 100),
                   genes.use = markers_use$gene, slim.col.label = TRUE,remove.key = TRUE))
  
}
#Clonotype-based function

#This is modified version of Plot_TR that takes into account showing stimulation vs non stimulation
plot_TR <- function(Meta, clonotype, TR_chain = TRB_CDR3, color = "red",background = "grey80",
                    pt_size = 0.5, bk_size = 0.5, ...){
  TR_chain <- enquo(TR_chain)
  TR_plot <- ggplot(Meta, aes(x = Tsne.x, y = Tsne.y)) + geom_point(size = bk_size, stroke = 0, pch = 16, color = background) +
    geom_point(data = Meta %>% filter(!!TR_chain %in% clonotype), 
               aes(x = Tsne.x, y = Tsne.y, shape = Stimulation), size = pt_size, color = color,stroke = 0, ...)+
    scale_shape_manual(values=c(1, 16))
  return(TR_plot)
}

# Load in the seurat and meta dataset -------------------------------------
# this section deals with combining the stim and non-stim datasets. This was done because the analysis started out with 
# just the stim. Non stim data was added during review process. Skip and directly load the combined seurat objects if
# wanted.
Stim_meta <- read_rds("../E7s_meta_data.rds")
Nonstim_meta <- read_rds("../E7_meta_data.rds")

Stim_seurat <- read_rds("../E7s_Seurat.rds")
Nonstim_seurat <- read_rds("../E7_Seurat.rds")

#try to clean up some of the meta data columns to make things a bit easier
Stim_seurat@meta.data$orig.ident <- Stim_seurat@meta.data$orig.ident.2
Stim_seurat@meta.data$cluster.tet.s <- Stim_seurat@meta.data$res.3


#Also put the conditions in so it would be clear where the data are actually from later

Stim_seurat@meta.data$Stimulation <- TRUE
Nonstim_seurat@meta.data$Stimulation <- FALSE

Stim_seurat@meta.data$Sample <- str_c(Stim_seurat@meta.data$orig.ident,"stim",sep=".")
Nonstim_seurat@meta.data$Sample <- str_c(Nonstim_seurat@meta.data$orig.ident,"nonstim",sep = ".")

# this is necessary due to I think a bug in this version of seurat. but in any case, resetting the object works to help
# with the merge command later on

Stim_seurat <- SubsetData(Stim_seurat)
Nonstim_seurat <- SubsetData(Nonstim_seurat)

# Merge the data, and run standard analysis pipeline, make clustering heatmap ----------------------

Combined_seurat <- MergeSeurat(Stim_seurat,Nonstim_seurat,add.cell.id1 = "stim",add.cell.id2 = "nonstim")

View(Combined_seurat@meta.data)

# alternatively, we can also load the already combined seurat objects. This is larger object because all genes were
# scaled (as shown later in the script)
# Combined_seurat <- read_rds("../Combined_all_scaled.rds")


rm(list= c('Stim_seurat','Nonstim_seurat')) #remove files we don't need to reduce space usage

variance.cutoff <- 1 # default
mean.cutoff <- 0.1 # default

Combined_seurat <- Seurat_process(Combined_Seurat, mf = mean.cutoff, cf = variance.cutoff, umap.run = FALSE)
Combined_seurat <-  RunTSNE(Combined_seurat, dims.use = 1:5, perplexity = 100)

# testing different resolution values to adjust cluster calling

for(i in seq(from = 0.4, to = 1.0 ,by = 0.1)){
  Combined_seurat <- FindClusters(Combined_seurat,dims.use=1:5,resolution=i, print.output = 0)
}

Combined_seurat <- SetAllIdent(Combined_seurat, id = 'res.0.4')

# we found that res.0.4 is appropriate for this analysis

# find genes for heatmap generation.

markers_all_0.4 <- FindAllMarkers(SubsetData(Combined_seurat,max.cells.per.ident = 200),test.use = "roc", do.print = TRUE)
markers_use_0.4 <- markers_all_0.4 %>% filter(avg_diff>0) %>% 
  mutate(cluster = factor(cluster, levels = c(0,1,3,2,4,5))) %>% arrange(cluster, desc(avg_diff))


DoHeatmap(SubsetData(Combined_seurat,max.cells.per.ident = 200),
          genes.use = intersect(markers_use_0.4$gene,Combined_seurat@var.genes),
          slim.col.label = TRUE,remove.key = TRUE,
          disp.min = -2.5, disp.max = 2.5, cex.row = 6, group.cex = 8, group.spacing = 0.05,
          group.order = c(0,1,3,2,4,5))

# in the end, clustering using resolution of 0.4 was what gave us the best results.

# remove irrelevant clusters to clean up heatmap

cells_hi_gene_reads<- WhichCells(SetAllIdent(Combined_seurat, id = 'res.1'),ident = 14)

cleaned_cells <- setdiff(Combined_seurat@cell.names,cells_hi_gene_reads)

DoHeatmap(SubsetData(Combined_seurat,max.cells.per.ident = 200,cells.use = cleaned_cells),
          genes.use = intersect(markers_use_0.4$gene,Combined_seurat@var.genes),
          slim.col.label = TRUE,remove.key = FALSE,
          disp.min = -2.5, disp.max = 2.5, cex.row = 5, group.cex = 8, group.spacing = 0.1,
          group.order = c(0,1,3,2,4,5))

ggsave("ClusterHeatmap_6_groups.pdf", width = 5, height = 4.5, units = "in")

# re-run with lower cutoff to find more variable genes, and see if better clustering--------

# testing if we could use a more relaxed set of cutoff to find more variable genes, and improve the clustering results. 

variance.cutoff <- 0.6
mean.cutoff <- 0.05


Combined_seurat_0.6 <- Seurat_process(Combined_seurat, mf = mean.cutoff, cf = variance.cutoff, umap.run = FALSE)
Combined_seurat_0.6 <-  RunTSNE(Combined_seurat_0.6, dims.use = 1:5, perplexity = 100)

for(i in seq(from = 0.8, to = 0.4 ,by = -0.2)){
  Combined_seurat_0.6 <- FindClusters(Combined_seurat_0.6,dims.use=1:5,resolution=i, print.output = FALSE,
                                      force.recalc = TRUE)
}

markers_all_0.4_2 <- FindAllMarkers(SubsetData(Combined_seurat_0.6,max.cells.per.ident = 200),test.use = "roc", do.print = TRUE)
markers_use_0.4_2 <- markers_all_0.4_2 %>% filter(avg_diff>0) %>% 
  mutate(cluster = factor(cluster, levels = c(0,1,3,6,2,4,5))) %>% arrange(cluster, desc(avg_diff))


DoHeatmap(SubsetData(Combined_seurat_0.6,max.cells.per.ident = 200),
            genes.use = intersect(markers_use_0.4_2$gene,Combined_seurat_0.6@var.genes),
            slim.col.label = TRUE,remove.key = TRUE,
            disp.min = -2.5, disp.max = 2.5, cex.row = 6, group.cex = 8, group.spacing = 0.05,
            group.order = c(0,1,3,2,4,5,6))
#The results show that there are really not much more genes with the relaxed plotting data, and the plotting looks almost 
#exactly same with the less relaxed parameter. 

# Plot detailed tSNE with cluster and treatment ---------------------------

# add tsne coordinates to the seurat object and then pull out the meta data for easy plotting
Combined_seurat <- TSNE_xy(Combined_seurat)
Combined_seurat <- Add_TRAB_count(Combined_seurat)
Combined_seurat@meta.data$cluster <- Combined_seurat@ident # make a cluster category that just stores the cluster idents

Combined_seurat_Meta <- Combined_seurat@meta.data


ggplot(Combined_seurat_Meta, aes(x = Tsne.x, y = Tsne.y, color = res.0.4, shape = Stimulation)) + geom_point(stroke = 0, size = 1) +
  scale_shape_manual(values=c(1, 16)) + Axis_themes + No_axis_labels + labs(x = 'tSNE 1', y = 'tSNE 2', color = 'cluster')

# testing figure size to help with visibility

ggsave('Combined_tsne_6_groups.pdf', width = 4.5, height = 4, units = 'in', family = 'ArialMT')

# subsample to help with visibility
ggplot(Combined_seurat_Meta %>% sample_n(nrow(Combined_seurat_Meta)), aes(x = Tsne.x, y = Tsne.y, color = res.0.4)) + geom_point(stroke = 0, size = 0.5) +
  scale_shape_manual(values=c(1, 16)) + Axis_themes + No_axis_labels + labs(x = 'tSNE 1', y = 'tSNE 2', color = 'cluster')

ggsave('Combined_tsne_6_groups_small.pdf', width = 3.5, height = 3, units = 'in', family = 'ArialMT')
# Like here, I think people will wonder about stimulation vs non-stim in each of the cluster

ggplot(Combined_seurat_Meta %>% sample_n(nrow(Combined_seurat_Meta)), aes(x = Tsne.x, y = Tsne.y, color = orig.ident)) + geom_point(stroke = 0, size = 0.5) +
  scale_shape_manual(values=c(1, 16)) + Axis_themes + No_axis_labels + labs(x = 'tSNE 1', y = 'tSNE 2', color = 'animal')

ggsave('Combined_tsne_animals_small.pdf', width = 3.5, height = 3, units = 'in', family = 'ArialMT')

#The reviewers indicated that they thought the color was hard to see. Let's plot one of just the stimulation condition

ggplot(Combined_seurat_Meta %>% sample_n(nrow(Combined_seurat_Meta)), aes(x = Tsne.x, y = Tsne.y, color = Stimulation)) + geom_point(stroke = 0.1, size = 0.5) +
  scale_color_manual(values=c('grey70', 'grey20')) + Axis_themes + No_axis_labels + labs(x = 'tSNE 1', y = 'tSNE 2', color = 'Stimulation')

ggsave('Combined_tsne_animals_small_conditions.png', width = 3.5, height = 3, units = 'in', family = 'ArialMT')
ggsave('Combined_tsne_animals_small_conditions.pdf', width = 3.5, height = 3, units = 'in', family = 'ArialMT')

# Add TCR information, plot expansion ---------------------------------------------

TRB_breaks <- log2(2^(c(0,seq(10))))
TRB_breaks_label <- 2^(c(0,seq(10)))

#the part that might need to change is the 
TSNE_TRB_expan <- TR_expan_map(Combined_seurat_Meta, ft.size = 0.5, stroke = 0, shape = 16) + scale_color_viridis(direction = -1, name = "Clonal size \n(TCRB)", breaks = TRB_breaks, 
                                                                                 labels = TRB_breaks_label, limits = c(0,log2(1024))) + labs(x = "tSNE 1", y = "tSNE 2") + Axis_themes

TSNE_Stimulation <- ggplot(Combined_seurat_Meta %>% sample_n(14000), aes(x = Tsne.x, y = Tsne.y, color = Stimulation)) + geom_point(size = 0.2) + 
  labs(x = "tSNE 1", y = "tSNE 2") + Axis_themes

ggsave("TSNE_Stim.eps", plot = TSNE_Stimulation, device = "eps",width = 3, height = 2.5, units = "in", family = "ArialMT")
ggsave("TSNE_Stim.png", plot = TSNE_Stimulation, device = "png",width = 3, height = 2.5, units = "in")

ggsave("TSNE_TRB_expan.png", plot = TSNE_TRB_expan + No_legend + No_axis_labels, device = "png",width = 2.5, height = 1.7, units = "in")
ggsave("TSNE_TRB_expan.pdf", plot = TSNE_TRB_expan + No_legend + No_axis_labels, device = "pdf",width = 2, height = 2, units = "in", family = "ArialMT")

# Plot canonical genes ----------------------------------------------------
Naive_genes <- c("CCR7","TCF7","SELL","CD44")
Eff_Em_genes <- c("GZMB","PRF1","ID2","GZMA","KLRG1")
Activation_exhaustion <- c("PDCD1","LAG3","TIGIT","CTLA4")

gene_list <- list(Naive_genes,Eff_Em_genes,Activation_exhaustion)
names(gene_list) <- c("Naive","Eff_em","Act_exh")

Selected_T_cell_genes <- map(gene_list,function(gene_list) FeaturePlot(Combined_seurat, gene_list, pt.size = 0.2, do.return = TRUE))
names(Selected_T_cell_genes) <- c("Naive","Eff_em","Act_exh")

Selected_T_cell_genes <- map(Selected_T_cell_genes,function(Selected_T_cell_genes) map(Selected_T_cell_genes,`+`,TSNE_theme))

plot_grid(plotlist = Selected_T_cell_genes[["Naive"]], nrow = 1)

save_plot("Naive_genes.png", plot_grid(plotlist = Selected_T_cell_genes[["Naive"]],nrow = 1),ncol = 4, base_height = 1.2, base_aspect_ratio = 0.9, unit = "in")
save_plot("Naive_genes.pdf", plot_grid(plotlist = Selected_T_cell_genes[["Naive"]],nrow = 1),ncol = 4, base_height = 1.2, base_aspect_ratio = 0.9, unit = "in", family = "ArialMT")

save_plot("Eff_Em_genes.png", plot_grid(plotlist = Selected_T_cell_genes[["Eff_em"]],nrow = 1),ncol = 5, base_height = 1.2, base_aspect_ratio = 0.9, unit = "in")
save_plot("Eff_Em_genes.pdf", plot_grid(plotlist = Selected_T_cell_genes[["Eff_em"]],nrow = 1),ncol = 5, base_height = 1.2, base_aspect_ratio = 0.9, unit = "in", family = "ArialMT")

save_plot("Act_exh_genes.png", plot_grid(plotlist = Selected_T_cell_genes[["Act_exh"]],nrow = 1),ncol = 4, base_height = 1.2, base_aspect_ratio = 0.9, unit = "in")
save_plot("Act_exh_genes.pdf", plot_grid(plotlist = Selected_T_cell_genes[["Act_exh"]],nrow = 1),ncol = 4, base_height = 1.2, base_aspect_ratio = 0.9, unit = "in", family = "ArialMT")


test[[1]] + TSNE_theme + theme(legend.position = 'right', legend.text = element_text(size = 6),
                               legend.title = element_blank(), legend.key.size = unit(0.2,'in'),
                               legend.key.height = unit(0.8,'in'))

# see if we can get everything on the same scale

Selected_T_cell_genes_test <- map(Selected_T_cell_genes,function(Selected_T_cell_genes) map(Selected_T_cell_genes,`+`,theme(legend.position = 'right')))

Selected_T_cell_genes_test <- map(Selected_T_cell_genes_test,
                                  function(Selected_T_cell_genes) map(Selected_T_cell_genes,`+`,
                                                                      scale_color_gradient(low = "yellow", 
                                                                                           high = "red",space = "Lab", na.value = "red", guide = "colourbar",aesthetics = "colour", limits = c(0,6))))

# The best thing to do in this case I think is to scale each row separately, since the cytotoxic cells are just a lot more highly expressed. than anyone else.

scale_naive = c(0,4.5)
scale_Eff_em = c(0,6)
scale_Act_exh = c(0,4.5)
color_scale <- list(scale_naive, scale_Eff_em, scale_Act_exh)

Selected_T_cell_genes_test <- map(Selected_T_cell_genes,function(Selected_T_cell_genes) map(Selected_T_cell_genes,`+`,theme(legend.position = 'right', legend.text = element_text(size = 6),
                                                                                                                            legend.title = element_blank(), legend.key.size = unit(0.1,'in'),
                                                                                                                            legend.key.height = unit(0.1,'in'))))

test <- map2(Selected_T_cell_genes_test,color_scale,
             function(Selected_T_cell_genes_test,color_scale) map(Selected_T_cell_genes_test,`+`,
            scale_color_gradient(low = "yellow", high = "red",space = "Lab", na.value = "red", 
                                 guide = "colourbar",aesthetics = "colour", limits = color_scale)))


save_plot("Naive_genes_legend.pdf", plot_grid(plotlist = test[["Naive"]],nrow = 1),ncol = 4, base_height = 1.2, base_aspect_ratio = 0.9, unit = "in", family = "ArialMT")
save_plot("Eff_Em_genes_legend.pdf", plot_grid(plotlist = test[["Eff_em"]],nrow = 1),ncol = 5, base_height = 1.2, base_aspect_ratio = 0.9, unit = "in", family = "ArialMT")
save_plot("Act_exh_genes_legend.pdf", plot_grid(plotlist = test[["Act_exh"]],nrow = 1),ncol = 4, base_height = 1.2, base_aspect_ratio = 0.9, unit = "in", family = "ArialMT")

test <- map(test,function(test) map(test,`+`,theme(legend.position = 'none')))

save_plot("Naive_genes_scaled.pdf", plot_grid(plotlist = test[["Naive"]],nrow = 1),ncol = 4, base_height = 1.2, base_aspect_ratio = 0.9, unit = "in", family = "ArialMT")
save_plot("Eff_Em_genes_scaled.pdf", plot_grid(plotlist = test[["Eff_em"]],nrow = 1),ncol = 5, base_height = 1.2, base_aspect_ratio = 0.9, unit = "in", family = "ArialMT")
save_plot("Act_exh_genes_scaled.pdf", plot_grid(plotlist = test[["Act_exh"]],nrow = 1),ncol = 4, base_height = 1.2, base_aspect_ratio = 0.9, unit = "in", family = "ArialMT")

# Trying to do clonotype based DE -----------------------------------------

# The restuls here is not useful.In the end, the direct DE approach by clones was not super
# clean. looking for modulce and clusters of genes were more useful.

Combined_seurat@meta.data$cell_id <- Combined_seurat@cell.names
Combined_seurat@meta.data$cluster <- Combined_seurat@ident
Combined_seurat_Meta <- Combined_seurat@meta.data
Combined_seurat_Meta <- as_tibble(Combined_seurat_Meta)

# plot it in heatmap
TRB_cluster_hmp <- TR_cluster_heatmap(Combined_seurat_Meta,cluster_var = cluster)
TRB_cluster_hmp_w <- TR_heatmap_spread(TRB_cluster_hmp,cluster_var = cluster)

pheatmap(TRB_cluster_hmp_w)

# create a directory to throw all the DE results and plot into
dir.create("ClusterDifGenes")
Sys.chmod("ClusterDifGenes", mode = "0777", use_umask = TRUE)

markers_all <- FindAllMarkers(Combined_seurat,test.use = "roc", do.print = TRUE)
write_csv(markers_all,"ClusterDifGenes/Cleaned_Cluster_Enrichments.txt")

markers_use <- markers_all %>% filter(avg_diff>0&power>0.3) %>% select(gene) %>% unique()
write_csv(markers_use,"ClusterDifGenes/Cleaned_Cluster_Unique_markers_use.txt")

markers_use <- markers_all %>% filter(avg_diff>0&power>0.3) 

DoHeatmap(SubsetData(Combined_seurat,max.cells.per.ident = 100),
          genes.use = markers_use$gene, slim.col.label = TRUE,remove.key = TRUE)
ggsave("ClusterHeatmap.pdf", path = "ClusterDifGenes", width = 10, height = 8, units = "in")

clusters <- levels(Combined_seurat@ident)

# re-write our automatic plot generation to make it more clear

Cluster_markers <- map(clusters,function(clusters) markers_all$gene[markers_all$cluster == clusters]) %>%
  map(function(genes) genes[1:9]) %>% map(function(genes) genes[!is.na(genes)]) %>% compact() %>%
  map(function(genes) FeaturePlot(Combined_seurat,genes,do.return = TRUE))

plot_names <- map(clusters, function(clusters) paste0("ClusterDifGenes/MarkersClust", clusters, ".pdf"))

#this part needs to be fixed. After compact() above, the nubmer of clsuters wasn't varied.
pwalk(list(plot_names[2:15],Cluster_markers),
      function(plot_names,Cluster_markers) save_plot(plot_names,plot_grid(plotlist=Cluster_markers,ncol = 3, nrow = 3),base_height = 11))

#get the cells that are stim_cytotoxic

Combined_seurat_Meta %>% filter(Stimulation) %>% filter(cluster.tet.s == "Cytotoxic")

test <- Combined_seurat_Meta %>% count(TRB_CDR3) %>% filter(!is.na(TRB_CDR3)) %>% arrange(desc(n))

stim_top <- get_cell_id(Combined_seurat_Meta,stimulation = TRUE,clonotype = "CASSQDLGNYAEQFF",TR_column = TRB_CDR3)

nonstim_top <- get_cell_id(Combined_seurat_Meta,stimulation = FALSE,clonotype = "CASSQDLGNYAEQFF",TR_column = TRB_CDR3)

Combined_seurat <- SetIdent(Combined_seurat, cells.use = stim_top, ident.use = "stim_top1")
Combined_seurat <- SetIdent(Combined_seurat, cells.use = nonstim_top, ident.use = "nonstim_top1")

topclone_DE<- DE_test(Combined_seurat,group_1 = "stim_top1",group_2 = "nonstim_top1")

#the results show some differences, Some of it is not what I expected. For this clone, for example, there's a 
#large difference in IFNG and CCL3/4. IFNg might be just one clone, but the other ones are more likely to be 
#clear differences between the clones.

# Plot the dot plot for expansion -----------------------------------------

# create a temp objects so we don't modify the original object.
temp <- Combined_seurat_Meta %>% filter(TRB_count>=1) %>% distinct(TRB_CDR3,TRB_count)
temp <- temp %>% count(TRB_count)

ggplot_dot(x = temp$TRB_count, count = temp$n, pch = 21, color = "black", size = 0.5, stroke = 0.1) + 
  scale_x_log10(breaks = TRB_breaks_label, limits = c(2,1024)) + 
  scale_fill_viridis(direction = -1,breaks = TRB_breaks, labels = TRB_breaks_label,limits = c(0,log2(1024))) + 
  labs(x = "TRB Clonal size", y = "Count of unique clones", fill = "TRB Clonal size") + Axis_themes + ylim(0,70)

#the dot plot still looks better, but on a really small element, it's kind difficult, We'd just have to see in the
#final figure before we really make a decision.

#Just hasn't really come up with a good solution to this. I think it will need to be re-plotted
#I think the only way to do this is to plot faceted plots, and make the scale clear. Would also help to make the
#plot easier to read in a smaller space.
TRB_breaks_label_y <- c(1,2,3,TRB_breaks_label[3:11]+1)
log2(seq(1,1001,100))
ggplot(temp, aes(x = TRB_count, y = n, fill = log2(TRB_count))) + geom_col(width = 0.02,color = "grey87") +
  scale_x_log10(breaks = TRB_breaks_label, limits = c(0.8,1024)) + theme(panel.grid.major.y = element_line(colour = "grey80", size = 0.1)) +
  scale_y_continuous(breaks = TRB_breaks_label, labels = TRB_breaks_label) +coord_cartesian(ylim=c(0, 32)) +
  scale_fill_viridis(direction = -1, breaks = TRB_breaks,labels = TRB_breaks_label,limits = c(0,log2(1024))) + 
  labs(x = "TRB Clonal size", y = "Count of unique clones", fill = "TRB Clonal size") + Axis_themes

ggsave("dotplot.png",plot = last_plot(), device = "png",width = 4, height = 1.5, units = "in")
ggsave("TRB_expan_dotplot.eps",plot = last_plot(), device = "eps",width = 4, height = 1, units = "in", family = "ArialMT")

#I think the bar plot just doesn't really make sense at this point, Let's see if we can make something that's
#a little more readable by using the beeswarm plot.

temp3 <- Combined_seurat_Meta %>% filter(!is.na(TRB_count)) %>% select(TRB_CDR3,TRB_count,Stimulation) %>% distinct()
TRB_count_by_stim <- Combined_seurat_Meta %>% filter(!is.na(TRB_count)) %>% group_by(Stimulation) %>% count(TRB_CDR3)
  
TRB_count_by_stim <- TRB_count_by_stim %>% rename(count = n)
TRB_count_by_stim <- TRB_count_by_stim %>% mutate(count = as.numeric(count)) %>% ungroup()
text_data <- tibble(TRB_count = rep(c(1,2),2), Stimulation = c(T,T,F,F),value = c(466,33,558,32))

ggplot(TRB_count_by_stim, aes(x = Stimulation, y = log2(count), fill = log2(count))) + 
  geom_quasirandom(pch = 21, stroke = 0,color = 'grey70', varwidth = TRUE, size = 0.75,
                   method ="tukey") +
  scale_y_continuous(breaks = log2(TRB_breaks_label), labels = TRB_breaks_label)+
  scale_fill_viridis(direction = -1, breaks = TRB_breaks, labels = TRB_breaks_label, limits = c(0,log2(1024))) +
  Axis_themes + geom_text(data = text_data, aes(x = Stimulation, y = log2(TRB_count), label = value), 
                          size = 2,vjust = -0.5, angle = -90) +
coord_flip() + geom_hline(yintercept = log2(15), linetype = 'dashed') + No_legend


ggsave("dotplot_tukey.png",plot = last_plot(), device = "png",width = 3, height = 1.25, units = "in")
ggsave("dotplot_tukey.pdf",plot = last_plot(), device = "pdf",width = 3, height = 1.25, units = "in", family = 'ArialMT')  

# I think the beeswarm plot is good
temp2 <- Combined_seurat_Meta %>% select(TRB_CDR3, cluster,TRB_count) %>% filter(!is.na(TRB_CDR3)) #%>% distinct()

ggplot(temp2, aes(x = cluster, y =  log2(TRB_count), fill = cluster)) + geom_violin(scale = 'width') + theme_classic() +
  Axis_themes + labs(x = 'cluster',y = 'log2(TRB clonal size)', fill = 'cluster') + 
  scale_y_continuous(breaks = TRB_breaks, labels = TRB_breaks_label)

#looking at the plot, it's obvious that cluster 2 and 8 are the navie clsuter, but looking at cluster 6-4, it's no longer super
#clean from the other clusters. so this would not allow us to look at the clonotype specific to expansion level quickly



#I think it's just going to be hard to draw boundaries. One way is to not draw boundaries,
#and then try to look for some sort of enrichment of genes based on some order of cells. but
#that's really just praying for gene sthat separate these out. which I think might work to some 

# plot histogram of clonal expansion and count ----------------------------

ggplot(Combined_seurat_Meta, aes(x = (TRB_count))) + geom_freqpoly(binwidth = 0.1)

ggplot(Combined_seurat_Meta, aes(x = Stimulation, y = (TRB_count))) + geom_violin(scale='width')

ggplot(Combined_seurat_Meta, aes(y = (TRB_count))) + geom_freqpoly(binwidth = 0.1)

# Gene to clonotype mapping, initial pass -----------------------------------------------
#I just have to see if anything would show up to the same list of genes

heatmap_genes <- read_csv("genes_stim.csv")

#incrase TRB_count cutoff 10 to reduce number of clones to look at.
exp_cells <- Combined_seurat_Meta %>% select(TRB_CDR3, cluster,TRB_count) %>% filter(!is.na(TRB_CDR3))

# let's try to see with #15 cutoff if we still see the same results. My fear is that the Myc-activated cells would be
# several impacted.
exp_cells <- exp_cells %>% filter(TRB_count >=15) %>% arrange(desc(TRB_count)) %>% select(TRB_CDR3) %>% 
  distinct() %>% pull(TRB_CDR3)

# make a separate seurat object that has all the genes scaled, so that we don't need to worry about not finding the right
#scaled genes, since only variable genes are scaled by default and missing out as a result of that. 
Heatmap_seurat <- ScaleData(Combined_seurat, vars.to.regress = c('nUMI', 'percent.mito'), 
                            model.use = 'poisson', genes.use = NULL)

#actually we have to break the data into stim and non-stim. Otherwise it won't really make sense
#remember that later that naive cluster was changed, so this isn't quite the same anymore.

naive_single <- Combined_seurat_Meta %>% filter(Stimulation == TRUE,cluster == 8, TRB_count<=1) %>% select(TRB_CDR3)

# To make it more easily plotted, we'll randomly sample 10-15 of these cells, and add that to the end.

sampled_single <- naive_single %>% sample_n(15)

clone_list <- c(exp_cells,sampled_single$TRB_CDR3)

# now we can go ahead and make some of that heatmap with the same heatmap_genes...
clone_id <- get_cell_id(Combined_seurat_Meta,stimulation = TRUE,clone_list)
# Combined_seurat_Meta %>% filter(Stimulation == TRUE,TRB_CDR3 %in% clone_list) %>% select(cell_id,TRB_CDR3)

stim_hmp<- clonotype_gene_data(Heatmap_seurat,genes = heatmap_genes$value,cells_id = clone_id$cell_id,
                    TRB_list = clone_id$TRB_CDR3,TRB_order = clone_list)

pheatmap(stim_hmp, cluster_cols = TRUE, cluster_rows = TRUE)

#now let's run it again real quick, with the nonstim cells, and just see if we can get the same thing

naive_single_ns <- Combined_seurat_Meta %>% filter(Stimulation == FALSE,cluster == 8, TRB_count<=1) %>% select(TRB_CDR3)
sampled_single_ns <- naive_single_ns %>% sample_n(15)

clone_list_ns <- c(exp_cells,sampled_single_ns$TRB_CDR3)

clone_id_ns <- get_cell_id(Combined_seurat_Meta,stimulation = FALSE,clone_list_ns)

#plot heatmap. this section should be redone

stim_hmp_ns<- clonotype_gene_data(Heatmap_seurat,genes = heatmap_genes$value,cells_id = clone_id_ns$cell_id,
                               TRB_list = clone_id_ns$TRB_CDR3,TRB_order = clone_list_ns)

pheatmap(stim_hmp_ns, cluster_cols = TRUE, cluster_rows = TRUE)

saveRDS(Heatmap_seurat,file = "Combined_all_scaled.rds")
rm(list = c("Heatmap_seurat"))
# in essense, this does work. and is showing what we want to see.
#couple things that are a bit different. One is that the expansion itself is not a perfect
#correlation, though it's a strong one. So we need to move away from that I think. We can just show
#the difference in expansion later.

# two: when we make two heatmaps, we should try to get the z score to correlate acros them... That's a bit hard.
#I think one way is to stack them, normalize, then split them. that way they are scale evenly across.




# trying to see if genes have good overlap with some other sets of variable genes, gene sets, etc. --------


length(heatmap_genes$value)

setdiff(heatmap_genes$value,Combined_seurat@var.genes)

FeaturePlot(Combined_seurat, c("CCR5"))

length(Combined_seurat@var.genes)

#hmm I wonder if the lenght of the variable genes are just too stringent.

test <- Seurat_process(Combined_seurat, mf = mean.cutoff, cf = 0.6, umap.run = FALSE,perp = 100)

length(test@var.genes)

length(setdiff(heatmap_genes$value,immunegenes$name))

#immune genes won't show much. 
immunegenes <- read.csv("InnateDB_genes.csv")

setdiff(heatmap_genes$value,markers_use$gene)
# Testing pulling clonotype gene data into heatmap into heatmap------------------------------------

stim_hmp_test<- clonotype_gene_data(Heatmap_seurat,genes = markers_all$gene %>% unique(),cells_id = clone_id$cell_id,
                               TRB_list = clone_id$TRB_CDR3,TRB_order = clone_list)

pheatmap(stim_hmp_test, cluster_cols = TRUE, cluster_rows = TRUE)

hmp_cluster_cols <- hclust(dist(t(stim_hmp_test)))

hmp_cluster_cols <- sort_hclust(hmp_cluster_cols)
plot(hmp_cluster_cols, main = "sorted Dendrogram", xlab = "", sub = "")

pheatmap(stim_hmp_test, cluster_cols = hmp_cluster_cols, cluster_rows = TRUE)

hmp_cluster_rows <- hclust(dist((stim_hmp_test)))
hmp_cluster_rows <- sort_hclust(hmp_cluster_rows)

pheatmap(stim_hmp_test, cluster_cols = hmp_cluster_cols, cluster_rows = hmp_cluster_rows)

stim_hmp_test_ns<- clonotype_gene_data(Heatmap_seurat,genes = markers_all$gene %>% unique(),cells_id = clone_id_ns$cell_id,
                                  TRB_list = clone_id_ns$TRB_CDR3,TRB_order = clone_list_ns)

pheatmap(stim_hmp_test_ns, cluster_cols = hmp_cluster_cols, cluster_rows = TRUE, fontsize = 4)

#I'm not actually sure if looking at all cluster defining genes would help. and it's not really clear, and no 
#adjusted P-values here. what if we do a clean Seurat, and then 


Heatmap_seurat <- SetAllIdent(Heatmap_seurat,id = "cluster")

TSNEPlot(Heatmap_seurat)

group_order <-  c(4,6,7,13,3,5,9,8,2,0,1,10,11,12)

markers_use <- markers_all %>% filter(avg_diff>0&power>0.3) %>% 
  mutate(cluster = factor(cluster, levels = group_order)) %>% arrange(cluster, desc(avg_diff))

DoHeatmap(SubsetData(Heatmap_seurat,max.cells.per.ident = 100),
          genes.use = intersect(markers_use$gene,Heatmap_seurat@var.genes), 
          group.order = c(4,6,7,13,3,5,9,8,2,0,1,10,11,12,14), slim.col.label = TRUE,remove.key = TRUE,
          disp.min = -2, disp.max = 2, cex.row = 6, group.cex = 8, group.spacing = 0.05)

ggsave("ClusterHeatmap_nolegend.pdf", width = 3, height = 3, units = "in")
# messy heatmap not super helpful.

# Once we do this, we can see that there are a few clusters of differentiating and then we can I guess see if
#we do DE, if the results would be similar?

# DE_test based on 4 major groups of the heatmap to find all genes --------


Myc_DE <- DE_test(SubsetData(Heatmap_seurat,ident.remove = c(14,8)), group_1 = c(4,6,7,13))

genes_1_2 <- Myc_DE %>% filter(significance == "significant") %>% arrange(desc(avg_logFC)) %>% pull(genes)


nonstim_DE <- DE_test(SubsetData(Heatmap_seurat,ident.remove = c(14,8)), group_1 = c(3,5,9))

nonstim_genes <- nonstim_DE %>% filter(significance == "significant") %>% arrange(desc(avg_logFC)) %>% pull(genes)

stim_DE <- DE_test(SubsetData(Heatmap_seurat,ident.remove = c(14,8)), group_1 = c(0,1,2,10,11,12))

stim_genes <- stim_DE %>% filter(significance == "significant") %>% arrange(desc(avg_logFC)) %>% pull(genes)

naive_DE <- DE_test(SubsetData(Heatmap_seurat,ident.remove = c(14)), group_1 = c(8))

naive_genes <- naive_DE %>% filter(significance == "significant") %>% arrange(desc(avg_logFC)) %>% pull(genes)

gene_c <- unique(c(genes_1_2,nonstim_genes,stim_genes,naive_genes))

#let's work in the normalization, so that we're compareing evenly.


# Make clean clonotype-based heatmap-------------------------------
# I think first, let's find out how many clonotypes are shared across the Stimulation condition
# This creates the list of TRB that's shared between the two stimulation conditions

shared_TRB <- Combined_seurat_Meta %>% filter(!is.na(TRB_CDR3)) %>% select(TRB_CDR3, Stimulation) %>% distinct() %>% count(TRB_CDR3) %>%
  filter(n>1) %>% filter(!str_detect(TRB_CDR3,"\\*")) %>% pull(TRB_CDR3)

exp_cells <- intersect(exp_cells,shared_TRB)

#remember that later that naive cluster was changed, so this isn't quite the same anymore.
naive_single <- Combined_seurat_Meta %>% filter(Stimulation == TRUE,cluster == 8, TRB_count<=1) %>% select(TRB_CDR3)


sampled_single <- naive_single %>% sample_n(15)
sampled_single <- sampled_single %>% mutate(singleton_id = str_c("singleton",row_number()))


clone_list <- c(exp_cells,sampled_single$TRB_CDR3)

#Also remember that the naive clsuter has changed later. so the cluster also has to be re-run, if we need it to.
naive_single_ns <- Combined_seurat_Meta %>% filter(Stimulation == FALSE,cluster == 8, TRB_count<=1) %>% select(TRB_CDR3)
sampled_single_ns <- naive_single_ns %>% sample_n(15)

sampled_single_ns <- sampled_single_ns %>% mutate(singleton_id = str_c("singleton",row_number()))
clone_list_ns <- c(exp_cells,sampled_single_ns$TRB_CDR3)


clone_id <- get_cell_id(Combined_seurat_Meta,stimulation = TRUE,clone_list)
clone_id <- clone_id %>% left_join(sampled_single, by = "TRB_CDR3")
clone_id <- clone_id %>% mutate(TRB_CDR3 = if_else(!is.na(singleton_id),singleton_id,TRB_CDR3)) %>% 
  select(-singleton_id)

clone_id_ns <- get_cell_id(Combined_seurat_Meta,stimulation = FALSE,clone_list_ns)
clone_id_ns <- clone_id_ns %>% left_join(sampled_single_ns, by = "TRB_CDR3")
clone_id_ns <- clone_id_ns %>% mutate(TRB_CDR3 = if_else(!is.na(singleton_id),singleton_id,TRB_CDR3)) %>% 
  select(-singleton_id)

stim_hmp<- clonotype_gene_data(Heatmap_seurat,genes = gene_c, cells_id = clone_id$cell_id,
                                    TRB_list = clone_id$TRB_CDR3,
                                    TRB_order = c(exp_cells,str_c("singleton",c(1:15))),do_scale = FALSE)
stim_hmp_ns<- clonotype_gene_data(Heatmap_seurat,genes = gene_c, cells_id = clone_id_ns$cell_id,
                               TRB_list = clone_id_ns$TRB_CDR3,
                               TRB_order = c(exp_cells,str_c("singleton",c(1:15))),do_scale = FALSE)

combine <- rbind(stim_hmp,stim_hmp_ns)
combine <- scale(combine)

combine[which(combine<min)] <- min
combine[which(combine>max)] <- max

stim_hmp <- combine[1:81,]
stim_hmp_ns <- combine[82:162,]

#quick check to see how the data looks
pheatmap(stim_hmp, cluster_cols = TRUE, cluster_rows = TRUE, fontsize = 6)

hmp_cluster_cols <- hclust(dist(t(stim_hmp)), method = "ward.D2")
hmp_cluster_row <- hclust(dist((stim_hmp)), method = "ward.D2") #change method to ward.D2. Seems to work 
# much better, when it comes to making the clustering cleaner. That way we don't actually have to make a
# decision regarding which clonotype should go to which cluster.

# sorting the hclust dendrogram, so the results is a bit neater
hmp_cluster_cols <- sort_hclust(hmp_cluster_cols,method = 'min')
hmp_cluster_row <- sort_hclust(hmp_cluster_row)

plot(hmp_cluster_cols)

pheatmap(stim_hmp, cluster_cols = hmp_cluster_cols, cluster_rows = hmp_cluster_row, fontsize = 6)

pheatmap(stim_hmp_ns, cluster_cols = hmp_cluster_cols, cluster_rows = hmp_cluster_row, fontsize = 6)

#the row and col order of these two are indeed the same, which is nice. We can calculate a ratio easily

#the only difference would be we would have to be careful with 0's, and would probably have to deal in pseudocount
#or something or that nature, so that means we probably want to be a bit more careful about that.

#gene trees are cut into 4, since it seems there are another cluster of genes that change significantly between
#stim and unstim, though not really differentially between TCR.

gene_trees <- cutree(hmp_cluster_cols, k = 4)

gene_dgram_table <- data.frame(gene = gene_trees) 

clonotype_trees <- cutree(hmp_cluster_row, k = 3)
clonotype_dgram_table <- data.frame(TCR = clonotype_trees) 

heatmap_color$TCR <- brewer.pal(8,'Accent')[5:7]
names(heatmap_color$TCR) <- unique(clonotype_dgram_table$TCR)

heatmap_color$gene <- brewer.pal(8,'Accent')[1:4]
names(heatmap_color$gene) <- unique(gene_dgram_table$gene)


#Here plots the unfiltered, un 
pheatmap(stim_hmp, cluster_cols = hmp_cluster_cols, cluster_rows = hmp_cluster_row, 
         annotation_col = gene_dgram_table, annotation_row =clonotype_dgram_table, 
         annotation_color = heatmap_color,fontsize = 5,treeheight_row = 5, treeheight_col =5,
         show_colnames = FALSE, annotation_names_row = FALSE, annotation_names_col = FALSE,
         filename = "Stimulated_hmp.pdf", width = 4, height =4.7 ,legend = FALSE, annotation_legend = FALSE,
         cutree_rows = 3, cutree_cols = 4)

pheatmap(stim_hmp_ns, cluster_cols = hmp_cluster_cols, cluster_rows = hmp_cluster_row, 
         annotation_col = gene_dgram_table, annotation_row =clonotype_dgram_table,
         annotation_color = heatmap_color,fontsize = 5,treeheight_row = 5, treeheight_col = 5,
         show_colnames = FALSE, annotation_names_row = FALSE, annotation_names_col = FALSE,
         filename = "Non-stimulated_hmp.pdf", width = 4, height = 4.7,legend = FALSE, annotation_legend = FALSE,
         cutree_rows = 3, cutree_cols = 4)


# heatmap looks promising. Worth looking into more.
# finish the annotation of the heatmap, and see what we can get.
#then I think we need a way to get out the list of genes and TCR in a parallelized way, otherwise everything would take too long

Clonotype_by_tree <- vector(mode = 'list',length = length(unique(clonotype_dgram_table$TCR)))

Clonotype_by_tree <- map(unique(clonotype_dgram_table$TCR),function(test) Clonotype_by_tree[[test]] <- 
      rownames(clonotype_dgram_table)[which(clonotype_dgram_table==test)])

genes_by_tree <- vector(mode = 'list',length = length(unique(gene_dgram_table$gene)))

genes_by_tree <- map(unique(gene_dgram_table$gene),function(test) genes_by_tree[[test]] <- 
                           rownames(gene_dgram_table)[which(gene_dgram_table==test)])

# clonotype-gene mapping second pass --------------------------------------

#the issue with the current tree is that TCR group 1 is too large. should down sample the
#group so that there's more space on the figure.

table(clonotype_dgram_table)

#we want to sample the TCR in the 

Group1_sampled_TCR <- clonotype_dgram_table %>% rownames_to_column(var = 'CDR3') %>% 
  filter(TCR == 1) %>% sample_n(15) %>% pull(CDR3)

Group1_2_TCR <- clonotype_dgram_table %>% rownames_to_column(var = 'CDR3') %>% filter(!(TCR == 1)) %>%
  pull(CDR3)

#then we want to downsample the rows of the stim_hmp to the 

stim_hmp_short <- stim_hmp[c(Group1_sampled_TCR, Group1_2_TCR),]

stim_hmp_ns_short <- stim_hmp_ns[c(Group1_sampled_TCR, Group1_2_TCR),]

#There was a undone save, so we kinda have to re-input the the Group1_sampled_TCR, so that we can recreate
#the exact graph. 

Group1_sampled_TCR_2 <- c("CASSLSEYEQYF","CASSPWDSNERLFF","CASSLSSPSGNTLYF","CASSIRDIISNERLFF",
                          "CASSQGTGELYAEQFF","CATDRGNTEVFF","CASSLTEQDTQYF","CASSYWDNYAEQFF",
                          "CASSQDWGDSNERLFF","CASSSDDNYAEQFF","CASSRDWGDSAETLYF","CASSDAPGTGKTLYF",
                          "CASSRDPAWGDYAEQFF","CASSQDHGSAETLYF","CASSIWDSNQAPLF")

stim_hmp_short <- stim_hmp[c(Group1_sampled_TCR_2, Group1_2_TCR),]
stim_hmp_ns_short <- stim_hmp_ns[c(Group1_sampled_TCR_2, Group1_2_TCR),]

#worked, this was able to recap the previously sampled results.
hmp_cluster_cols_short <- hclust(dist(t(stim_hmp_short)), method = "ward.D2")
hmp_cluster_row_short <- hclust(dist((stim_hmp_short)), method = "ward.D2") #change method to ward.D2. Seems to work 
#much better, when it comes to making the clustering cleaner. That way we don't actually have to make a
#decision regarding which clonotype should go to which cluster.

#sorting the hclust dendrogram, so the results is a bit neater
hmp_cluster_cols_short <- sort_hclust(hmp_cluster_cols_short,method = 'min')
hmp_cluster_row_short <- sort_hclust(hmp_cluster_row_short)

clonotype_trees_short <- cutree(hmp_cluster_row_short, k = 3)
clonotype_dgram_table_short <- data.frame(TCR = clonotype_trees_short) 

pheatmap(stim_hmp_short, cluster_cols = hmp_cluster_cols, cluster_rows = hmp_cluster_row_short, 
         annotation_row =clonotype_dgram_table_short, annotation_col = gene_dgram_table, 
         annotation_color = heatmap_color, fontsize = 5,
         treeheight_row = 5, treeheight_col =5,show_colnames = FALSE, annotation_names_row = FALSE, 
         annotation_names_col = FALSE,legend = FALSE, annotation_legend = FALSE,cutree_rows = 3, 
         cutree_cols = 4, filename = "Stimulated_hmp_short.pdf", width = 3.5, height =3)

pheatmap(stim_hmp_ns_short, cluster_cols = hmp_cluster_cols, cluster_rows = hmp_cluster_row_short, 
         annotation_row =clonotype_dgram_table_short, annotation_col = gene_dgram_table, 
         annotation_color = heatmap_color, fontsize = 5,
         treeheight_row = 5, treeheight_col =5,show_colnames = FALSE, annotation_names_row = FALSE, 
         annotation_names_col = FALSE,legend = FALSE, annotation_legend = FALSE,cutree_rows = 3, 
         cutree_cols = 4, filename = "Stimulated_hmp_ns_short.pdf", width = 3.5, height =3)

#when we cut up the TCR, it ended up with a bunch of small groups, and 3 major groups, I think those are the
#ones we should pull do try to do DE.

#the short TCR is what we want to actually plot, so that the heatmap comes out more even.

genes_by_tree <- vector(mode = 'list',length = length(unique(gene_dgram_table$gene)))

genes_by_tree <- map(unique(gene_dgram_table$gene),function(test) genes_by_tree[[test]] <- 
                       rownames(gene_dgram_table)[which(gene_dgram_table==test)])



Group1_TCR <- clonotype_dgram_table %>% mutate(TRB = rownames(clonotype_dgram_table)) %>%
filter(clonotype_dgram_table,TCR == 1) %>% pull(TRB)

Group2_TCR <- clonotype_dgram_table %>% mutate(TRB = rownames(clonotype_dgram_table)) %>%
  filter(clonotype_dgram_table,TCR == 2) %>% pull(TRB)

Group3_TCR <- clonotype_dgram_table %>% mutate(TRB = rownames(clonotype_dgram_table)) %>%
  filter(clonotype_dgram_table,TCR == 3) %>% pull(TRB)

# Plot specific TRB on tSNE -----------------------------------------------

#We should pick specific clones from the heatmap to show on the inidividual 

Group2_clones <- c("CASSDSGGNTEVFF","CASGDWGGREQYF")
Group1_clones <- c("CASSLSEYEQYF","CASSQDWGDSNERLFF","CASSPWDSNERLFF")

TR_plot_list <- vector(mode = 'list', 5)

Group2_plots <- map(Group2_clones, function(Group2_clones) plot_TR(Combined_seurat_Meta,clonotype = Group2_clones, 
                                                                   pt_size =0.5, color = heatmap_color$TCR[[2]]) + TSNE_theme + labs(title = Group2_clones))
Group1_plots <- map(Group1_clones, function(Group1_clones) plot_TR(Combined_seurat_Meta,clonotype = Group1_clones, 
                                                                   pt_size =0.5, color = heatmap_color$TCR[[1]]) + TSNE_theme + labs(title = Group1_clones))


save_plot("TCR_clone.pdf", plot_grid(plotlist = c(Group1_plots, Group2_plots),nrow = 1),ncol = 5, base_height = 1.2, base_aspect_ratio = 0.9, unit = "in", family = "ArialMT")


# TCR based DE ------------------------------------------------------------
# unlike before, where we are comparing TCR directly across conditions. now we have groups of TCR, and we want to make
# comparison between the groups and conditions. This is a lot more powered.

Group1_s_cells <- get_cell_id(Combined_seurat_Meta,stimulation = TRUE,clonotype = Group1_TCR)
Group1_ns_cells <- get_cell_id(Combined_seurat_Meta,stimulation = FALSE,clonotype = Group1_TCR)

Group2_s_cells <- get_cell_id(Combined_seurat_Meta,stimulation = TRUE,clonotype = Group2_TCR)
Group2_ns_cells <- get_cell_id(Combined_seurat_Meta,stimulation = FALSE,clonotype = Group2_TCR)

# one of the things to note is that the number of cells are rather different. 

# DE_seurat is used for DE test, so we don't have to always switch between the various seurat objects

DE_seurat<- SetIdent(Heatmap_seurat, cells.use = Group1_s_cells$cell_id, ident.use = "Group1s")
DE_seurat<- SetIdent(DE_seurat, cells.use = Group1_ns_cells$cell_id, ident.use = "Group1ns")

DE_seurat<- SetIdent(DE_seurat, cells.use = Group2_s_cells$cell_id, ident.use = "Group2s")
DE_seurat<- SetIdent(DE_seurat, cells.use = Group2_ns_cells$cell_id, ident.use = "Group2ns")

# decided to not go for results in group 3. commented out.
#DE_seurat<- SetIdent(DE_seurat, cells.use = Group3_s_cells$cell_id, ident.use = "Group3s")
#DE_seurat<- SetIdent(DE_seurat, cells.use = Group3_ns_cells$cell_id, ident.use = "Group3ns")

TSNEPlot(DE_seurat, cells.use = c(Group1_s_cells$cell_id,Group1_ns_cells$cell_id,Group2_s_cells$cell_id,
                                  Group2_ns_cells$cell_id,Group3_s_cells$cell_id,Group3_ns_cells$cell_id))

# Group 1 is the cyto cells. group 2 is the myc cells.
Group1_DE <- DE_test(DE_seurat,group_1 = "Group1s",group_2 = "Group1ns", max_cells = 300)
Group2_DE <- DE_test(DE_seurat,group_1 = "Group2s",group_2 = "Group2ns", max_cells = 300)
Group1_2_ns_DE <- DE_test(DE_seurat, group_1 = "Group1ns",group_2 = "Group2ns", max_cells = 300)
Group1_2_s_DE <- DE_test(DE_seurat, group_1 = "Group1s",group_2 = "Group2s",max_cells = 300)

#Group3_DE <- DE_test(DE_seurat,group_1 = "Group3s",group_2 = "Group3ns")

DE_list <- list(Group1_DE,Group2_DE,Group1_2_s_DE, Group1_2_ns_DE)
names(DE_list) <- c("Group1","Group2","Group1_v_2s","Group1_v2ns")

DE_list <- map(DE_list, function(DE_list) DE_list %>% mutate(avg_logFC = avg_logFC * log10(exp(1))))
names(DE_list) <- c("Group1","Group2","Group1_v_2s","Group1_v2ns")
titles <- c("Group1","Group2","Group1_v_2s","Group1_v2ns")
DE_list_Volcano <- map2(DE_list,titles ,function(DE_list,titles) plot_volcano(DE_list, right_nudge = 1.5, left_nudge = -1) +
                         xlim(-1,1.5)+ theme(legend.position = 'none') + Axis_themes + labs(x = "log10(fold change)",
                                                                                            title=titles))


names(DE_list_Volcano) <- c("Group1","Group2","Group1_v_2s","Group1_v2ns")

plot_grid(plotlist = DE_list_Volcano, nrow = 2)

save_plot("Volcanoes.pdf",plot_grid(plotlist = DE_list_Volcano, nrow = 2),base_height = 5,family = "ArialMT")
# Repeat downsampled DE multiple times ------------------------------------


#let's test to see if we change the nuber of cells, we would get dramatically different results

Group1_DE_300 <- DE_test(DE_seurat,group_1 = "Group1s",group_2 = "Group1ns",max_cells = 300)

Group1_DE_300_volcano <- plot_volcano(Group1_DE_300,gene_list_1 = Group1_DE_300 %>% arrange(desc(p_val_adj),desc(avg_logFC)) %>% top_n(10,avg_logFC) %>% pull(genes),
                                  #gene_list_2 = Group1_DE_300 %>% arrange(desc(p_val_adj),(avg_logFC)) %>% top_n(10,-avg_logFC) %>% pull(genes))+ Axis_themes + labs(x = "log10(fold change)")

#the good news is that it doesn't seem to have changed the top listed genes very much, which is good news, considering that
#I haven't had to change a whole lot. and then everything is on the samescale without more input from me.
#let's do the same thing with group 2 cells as well.

Group2_DE_300 <- DE_test(DE_seurat,group_1 = "Group2s",group_2 = "Group2ns",max_cells = 300)

Group2_DE_300_volcano <- plot_volcano(Group2_DE_300,gene_list_1 = Group2_DE_300 %>% arrange(desc(p_val_adj),desc(avg_logFC)) %>% top_n(10,avg_logFC) %>% pull(genes),
                                      #gene_list_2 = Group2_DE_300 %>% arrange(desc(p_val_adj),(avg_logFC)) %>% top_n(10,-avg_logFC) %>% pull(genes))+ Axis_themes + labs(x = "log10(fold change)")

#would also seem like downsampling the cells in group2 also end up with the same sort of results as before. That makes it easier
#I think bimod just takes care of the cell imbalance. but a hard control is always better.

#let's extend downsampling test of the DE test to the non stim as well

Group1_3_ns_DE <- DE_test(DE_seurat, group_1 = "Group1ns",group_2 = "Group3ns", max_cells = 300)



Group1_3_s_DE <- DE_test(DE_seurat, group_1 = "Group1s",group_2 = "Group3s",max_cells = 300)



#also with stim, group 1 and 3, the large difference is more or less the same 

Group1_DE_300_list <- vector('list',10)
Rep <- c(1:10)
Group1_DE_300_list <- map(Rep,function(Rep) DE_test(DE_seurat,group_1 = "Group1s",group_2 = "Group1ns",max_cells = 300))

#seems like repeating the DE doesn't make that much of a difference. The results look extremely similar, or entirely
#the same. must be something to do with the way cells are sub-sampled, or the bi-mod test that

# plot differences in expression between group 1 and 2 TCR ------------------------------------------

#let's try a rough plot first. 

test <- full_join(Group1_DE,Group2_DE,by = "genes", suffix = c(".Grp1",".Grp2"))

test_filtered <- test %>% filter(!is.na(avg_logFC.Grp1) & !is.na(avg_logFC.Grp2))

test2 <- test
test2$avg_logFC.Grp1[is.na(test2$avg_logFC.Grp1)] <- 0

test_lm <- lm(avg_logFC.Grp2~avg_logFC.Grp1, data = test_filtered)
summary(test_lm)

ggplotRegression <- function (fit) {
  
  require(ggplot2)
  
  ggplot(fit$model, aes_string(x = names(fit$model)[2], y = names(fit$model)[1])) + 
    geom_point(size = 0.5, color = 'grey70') +
    stat_smooth(method = "lm", col = "red", size = 0.5) +
    labs(title = paste("Adj R2 = ",signif(summary(fit)$adj.r.squared, 5),
                       #"Intercept =",signif(fit$coef[[1]],5 ),
                       " Slope =",signif(fit$coef[[2]], 5),
                       " P =",signif(summary(fit)$coef[2,4], 5)))
}

Group1_2_FC <- ggplotRegression(test_lm)+ Axis_themes + geom_vline(xintercept = 0,linetype='dashed') + 
  geom_hline(yintercept=0,linetype='dashed') + labs(x = "Group1 log10(Fold change)",y = "Group2 log10(Fold change)")


#plot looks prety good, let's try to get an approximate of the final graph.
ggsave("Group1_2_stim_v_unstim.png", plot = Group1_2_FC, device = "png",width = 2.8, height = 3, units = "in")
ggsave("Group1_2_stim_v_unstim.pdf", plot = Group1_2_FC, device = "pdf",width = 3, height = 2.5, units = "in")



ggplot(test,aes(x = avg_logFC.Grp1,y = avg_logFC.Grp2)) + geom_point(size = 0.5, color = 'grey70') + 
  stat_smooth(method = "lm", col = "red") + Axis_themes + geom_vline(xintercept = 0,linetype='dashed') + geom_hline(yintercept=0,linetype='dashed')


#I think with Non-stim, we generally see the difference is in cytotoxic vs... something else. More naive? certainly we can say something about
#with stim, some of the difference is the same, but some of the others are very different. I guess the question is with the clonotype
#info, what are we highlighting here? that some cells seem to be in a different state of activation, based on their clonoality (and the quality
# of that clonality, such as expansion.)

# Looking at expansion level of the groups of TCR -------------------------

#I think in order to get a good tally of the total expansion of the TCR

TRB_count_by_stim <- Count_TRB(Combined_seurat_Meta, cutoff = 1, group_var = Stimulation)
TRB_count_by_stim <- TRB_count_by_stim %>% rename(Stimulation = sampleID)

TRB_count_by_stim <- TRB_count_by_stim %>% mutate(TRB_group = case_when(TRB_CDR3 %in% Group1_TCR ~ "Group1",
                                                                        TRB_CDR3 %in% Group2_TCR ~ "Group2"))

TRB_count_by_stim <- TRB_count_by_stim %>% mutate(Stimulation = case_when(Stimulation == TRUE ~ "Stimulated",
                                                                           Stimulation == FALSE ~ "Non-stimulated"))


my_comparisons <- list( c("Group1", "Group2"))

TRB_expan_group <- ggplot(TRB_count_by_stim %>% filter(!is.na(TRB_group)), aes(x = TRB_group, y = log2(TRB_count), fill= TRB_group)) + geom_violin(scale = 'width') + 
  facet_wrap(~Stimulation) + scale_fill_brewer(type = 'qual', palette = 'Accent') + 
  stat_compare_means(comparisons = my_comparisons, size = 2) + geom_boxplot(fill = 'white',width = 0.2) + 
  labs(x = '', y = 'Log2(TCRB clonal size)', fill = 'TCRB Group') + Axis_themes + 
  theme(strip.text = element_text(size = 6)) + theme(axis.text.x = element_text(angle = 90, hjust = 1))

ggsave("TRB_expan_between_groups.png", plot = TRB_expan_group, device = "png",width = 3.4, height = 2.6, units = "in")
ggsave("TRB_expan_between_groups.eps", plot = TRB_expan_group, device = "eps",width = 3.4, height = 2, units = "in",family = "ArialMT")


#can we get some test statistics on these populations? the TRB_count number is not really normaly distributed.
#we can log normalize it, but I think t test based on that may be kinda weird. So can we try a manwhittany test?

test <- TRB_count_by_stim %>% filter(Stimulation == 'Stimulated') %>% filter(TRB_group %in% c("Group1","Group3"))

wilcox.test(TRB_count~TRB_group,data = test)

#Group1 and 3 are significantly different but not between group 1 and 2, which is decent.

# module scores for cells--------------------------------------------------
# this si to calculate module scores in gene modules defined in Singer et al. the csv file is a direct port of one of 
# the supplement file from the paper. we'll get the genes and the cluster names from the supplement, and then use the
# AddModuleScore function to set the scoring for each groups of TCR that we have defined.

module_lists <- read_csv("../NIHMS813125_cluster.csv")

module_lists <- module_lists %>% rename(gene = x) %>% mutate(cluster = factor(as.character(cluster)))
module_lists <- module_lists %>% select(gene,cluster,everything())
module_lists <- module_lists %>% mutate(cluster = fct_relevel(cluster,"10",after = Inf))
# there are three thousand genes, and they are grouped into clusters as defined in the paper.

module_vector <- levels(module_lists$cluster)
cluster_names <- str_c("cluster",module_vector)

cluster_genes <- map(module_vector, function(module_vector) module_lists %>% filter(cluster == module_vector) %>% 
                       select(gene))

names(cluster_genes) <- cluster_names

spread(module_lists%>% select(gene,cluster),key = cluster, value = gene)

map2(cluster_genes,names(cluster_genes), function(x,y) write_csv(data.frame(y = x), 'test.csv', append= T,col_names = T))

# One thing that we'll have to look at is whether we can just first get an idea of how much overlap we have
# between the cluster and our dataset

All_genes <- rownames(DE_seurat@data)

cluster_genes <- map(cluster_genes, function(cluster_genes) intersect(cluster_genes$gene, All_genes))

#There is still a decent amount of overlap. can continue.
DE_seurat <- AddModuleScore(DE_seurat,genes.list = cluster_genes, ctrl.size = 10, enrich.name = "Cluster")

TRB_order <- c(exp_cells,str_c("singleton",c(1:30)))

DE_module_stim <- DE_seurat@meta.data %>% filter(cell_id %in% clone_id$cell_id) %>% rename(TRB_CDR3.old = TRB_CDR3) %>% left_join(clone_id, by = "cell_id") %>%
  mutate(TRB_CDR3 = factor(TRB_CDR3, levels = TRB_order)) %>% arrange(TRB_CDR3) %>% 
  select(cell_id, TRB_CDR3,starts_with("Cluster")) %>% as_tibble()

DE_module_nonstim <- DE_seurat@meta.data %>% filter(cell_id %in% clone_id_ns$cell_id) %>% rename(TRB_CDR3.old = TRB_CDR3) %>% left_join(clone_id_ns, by = "cell_id") %>%
  mutate(TRB_CDR3 = factor(TRB_CDR3, levels = TRB_order)) %>% arrange(TRB_CDR3) %>% 
  select(cell_id, TRB_CDR3,starts_with("Cluster")) %>% as_tibble()

#now that we have the DE modules in two separate dataframe, I think the best thing to do is collapse them first

DE_module_stim <- DE_module_stim %>% select(-cell_id,-cluster, -cluster.tet.s) %>% group_by(TRB_CDR3) %>% 
  summarise_all(funs(mean))

DE_module_stim <- DE_module_stim %>% column_to_rownames(var = "TRB_CDR3")

DE_module_stim <- as.matrix(DE_module_stim)

DE_module_nonstim <- DE_module_nonstim %>% select(-cell_id,-cluster, -cluster.tet.s) %>% group_by(TRB_CDR3) %>% 
  summarise_all(funs(mean))

DE_module_nonstim <- DE_module_nonstim %>% column_to_rownames(var = "TRB_CDR3")

DE_module_nonstim <- as.matrix(DE_module_nonstim)

#then we have to makesure we actually scale the module score across both conditions
#but first, let's switch them into modules
combine_module <- rbind(DE_module_stim,DE_module_nonstim)
combine_module <- scale(combine_module)

combine_module[which(combine_module<min)] <- min
combine_module[which(combine_module>max)] <- max

DE_module_stim <- combine_module[1:134,]
DE_module_nonstim  <- combine_module[135:268,]

DE_module_stim_short <- DE_module_stim[c(Group1_sampled_TCR_2, Group1_2_TCR),]
DE_module_nonstim_short <- DE_module_nonstim[c(Group1_sampled_TCR_2, Group1_2_TCR),]



hmp_module_cols <- hclust(dist(t(DE_module_stim_short)), method = "ward.D2")
hmp_module_cols <- sort_hclust(hmp_module_cols)

module_trees <- cutree(hmp_module_cols, k = 2)
module_dgram_table <- data.frame(module = module_trees)

heatmap_color$module <- brewer.pal(9,'Set1')[8:9]
names(heatmap_color$module) <- unique(module_dgram_table$module)

pheatmap(DE_module_stim_short, cluster_cols = hmp_module_cols, cluster_rows = hmp_cluster_row_short,
         annotation_col = module_dgram_table,annotation_row =clonotype_dgram_table_short,
         annotation_color = heatmap_color,fontsize = 6,treeheight_row = 0, treeheight_col = 10,
         show_colnames = TRUE, annotation_names_row = FALSE, annotation_names_col = FALSE, border_color = NA,
         cutree_rows = 3, cutree_cols = 2, width = 2.4, height = 5,legend = FALSE, annotation_legend = FALSE,
         filename = "stimulated_module_hmp_short.pdf")

pheatmap(DE_module_nonstim_short, cluster_cols = hmp_module_cols, cluster_rows = hmp_cluster_row_short,
         annotation_col = module_dgram_table,annotation_row =clonotype_dgram_table_short,
         annotation_color = heatmap_color,fontsize = 6,treeheight_row = 0, treeheight_col = 10,
         show_colnames = TRUE, annotation_names_row = FALSE, annotation_names_col = FALSE, border_color = NA,
         cutree_rows = 3, cutree_cols = 2, width = 2.4, height = 5,legend = FALSE, annotation_legend = FALSE,
         filename = "nonstimulated_module_hmp_short.pdf")

# MsigDb analysis and parsing ---------------------------------------------

# here we basically writes the genes in each cluster of genes, and the uploaded to MsigDB via their web interface.
# we then re-plot the results here.
write_csv(genes_by_tree[[1]] %>% as_tibble() %>% rename(Cluster1 = value), "Cluster1_genes.csv")
write_csv(genes_by_tree[[2]] %>% as_tibble() %>% rename(Cluster2 = value), "Cluster2_genes.csv")
write_csv(genes_by_tree[[3]] %>% as_tibble() %>% rename(Cluster3 = value), "Cluster3_genes.csv")
write_csv(genes_by_tree[[4]] %>% as_tibble() %>% rename(Cluster4 = value), "Cluster4_genes.csv")

#It's actually important to understand the gene set scores that are being enriched in each of the 3 gene clusters
#that we have identified

#Remember, cluster 1 is the Myc/cell cycle cluster. that's mapped against Hallmark (H) set
#cluster 2 is the Naive cluster, that's mapped against C2 and C7
#cluster 3 is the effector cluster, that's also mapped against C2 and C7
#updated the MsigDb results to the new separation of the gene lists. 

MsigDb <- read_csv("MsigDb/MsigDb_cleaned_2.csv")
MsigDb <- MsigDb %>% mutate(Cluster_subset = factor(Cluster_subset),Gene_Set_K = as.integer(Gene_Set_K),
                            log_qvalue = -log10(FDR_qvalue))

#filter the MsigDb. if we want unfiltered, can re-run lines above.
MsigDb <- MsigDb %>% filter(!is.na(Description_short)) %>% group_by(Cluster_subset) %>% arrange(Cluster_subset,desc(log_qvalue))

#add limits on the x axis so we're on the same scale.
#graph the results for cluster 1
MsigDb_1 <- ggplot(MsigDb %>% filter(Cluster_subset == "Cluster1") %>% slice(1:5),aes(x = fct_reorder(Description_short,log_qvalue), y = log_qvalue)) + 
  geom_bar(fill = heatmap_color[["gene"]][1], stat = "identity", width = 0.5) + coord_flip() + Axis_themes +
  labs(y = "-log10(FDR q-value)") + geom_hline(yintercept = -log10(0.01), linetype = 2) +
  scale_y_continuous(breaks =c(0,2,10,20,30), limits = c(0,25)) +theme(axis.title.y=element_blank())

ggsave("cluster1_Geneset.png",plot = MsigDb_1, device = "png",width = 3, height = 2, units = "in")
ggsave("cluster1_Geneset.pdf",plot = MsigDb_1, device = "pdf",width = 3, height = 2, units = "in",family = "ArialMT")
#graph the results for cluster 2
MsigDb_2 <- ggplot(MsigDb %>% filter(Cluster_subset == "Cluster2") %>% slice(1:5),aes(x = fct_reorder(Description_short,log_qvalue), y = log_qvalue)) + 
  geom_bar(fill = heatmap_color[["gene"]][2], stat = "identity", width = 0.5) + coord_flip() + Axis_themes +
  labs(y = "-log10(FDR q-value)") + geom_hline(yintercept = -log10(0.01), linetype = 2) +
  scale_y_continuous(breaks =c(0,2,10,20,30), limits = c(0,25)) +theme(axis.title.y=element_blank(),axis.text = element_text(size = 4))

ggsave("cluster2_Geneset.png",plot = MsigDb_2, device = "png",width = 4, height = 2, units = "in")
ggsave("cluster2_Geneset.pdf",plot = MsigDb_2, device = "pdf",width = 3, height = 2, units = "in",family = "ArialMT")
MsigDb_3 <- ggplot(MsigDb %>% filter(Cluster_subset == "Cluster3") %>% slice(1:5),aes(x = fct_reorder(Description_short,log_qvalue), y = log_qvalue)) + 
  geom_bar(fill = heatmap_color[["gene"]][3], stat = "identity", width = 0.5) + coord_flip() + Axis_themes +
  labs(y = "-log10(FDR q-value)") + geom_hline(yintercept = -log10(0.01), linetype = 2) +
  scale_y_continuous(breaks =c(0,2,10,20,30), limits = c(0,25)) +theme(axis.title.y=element_blank(),axis.text = element_text(size = 4))

ggsave("cluster3_Geneset.png",plot = MsigDb_3, device = "png",width = 4, height = 2, units = "in")
ggsave("cluster3_Geneset.pdf",plot = MsigDb_3, device = "pdf",width = 3, height = 2, units = "in",family = "ArialMT")

MsigDb_4 <- ggplot(MsigDb %>% filter(Cluster_subset == "Cluster4") %>% slice(1:5),aes(x = fct_reorder(Description_short,log_qvalue), y = log_qvalue)) + 
  geom_bar(fill = heatmap_color[["gene"]][4], stat = "identity", width = 0.5) + coord_flip() + Axis_themes +
  labs(y = "-log10(FDR q-value)") + geom_hline(yintercept = -log10(0.01), linetype = 2) +
  scale_y_continuous(breaks =c(0,2,10,20,30), limits = c(0,25)) +theme(axis.title.y=element_blank(),axis.text = element_text(size = 4))

ggsave("cluster4_Geneset.png",plot = MsigDb_4, device = "png",width = 4, height = 2, units = "in")
ggsave("cluster4_Geneset.pdf",plot = MsigDb_4, device = "pdf",width = 3, height = 2, units = "in",family = "ArialMT")

#then let's save just the bar portion, and then bring in the text later.
#I belive there's total of 2 inch in height and 7 in in width (so save more on the width for text sizes)


MsigDb_1 + theme(axis.text = element_blank(), axis.title = element_blank())
ggsave("cluster1_test.png",plot = last_plot(), device = "png",width = 2, height = 2, units = "in")
ggsave("cluster1_bars.pdf",plot = last_plot(), device = "pdf",width = 1.8, height = 1.2, units = "in",family = "ArialMT")

MsigDb_2 + theme(axis.text = element_blank(), axis.title = element_blank())
ggsave("cluster2_test.png",plot = last_plot(), device = "png",width = 2, height = 2, units = "in")
ggsave("cluster2_bars.pdf",plot = last_plot(), device = "pdf",width = 1.8, height = 1.2, units = "in",family = "ArialMT")

MsigDb_3 + theme(axis.text = element_blank(), axis.title = element_blank())
ggsave("cluster3_test.png",plot = last_plot(), device = "png",width = 2, height = 2, units = "in")
ggsave("cluster3_bars.pdf",plot = last_plot(), device = "pdf",width = 1.8, height = 1.2, units = "in",family = "ArialMT")

MsigDb_4 + theme(axis.text = element_blank(), axis.title = element_blank())
ggsave("cluster4_test.png",plot = last_plot(), device = "png",width = 2, height = 2, units = "in")
ggsave("cluster4_bars.pdf",plot = last_plot(), device = "pdf",width = 1.8, height = 1.2, units = "in",family = "ArialMT")


# Sifting GO term results -------------------------------------------------

MsigDb_GO <- read_csv("MsigDb/MsigDb_GO_cleaned.csv")
MsigDb_GO <- filter(MsigDb_GO, !is.na(Gene_Set_Name))
MsigDb_GO <- MsigDb_GO %>% mutate(Cluster_subset = factor(Cluster_subset),Gene_Set_K = as.integer(Gene_Set_K),
           log_qvalue = -log10(FDR_qvalue))

#I think just plotting the surface top ten go terms would be enough 

MsigDb_GO_2 <- ggplot(MsigDb_GO %>% filter(Cluster_subset == "Cluster2") %>% slice(1:10),aes(x = fct_reorder(Gene_Set_Name,log_qvalue), y = log_qvalue)) + 
  geom_bar(fill = heatmap_color[["gene"]][2], stat = "identity", width = 0.5) + coord_flip() + Axis_themes +
  labs(y = "-log10(FDR q-value)") + geom_hline(yintercept = -log10(0.01), linetype = 2) +
  scale_y_continuous(breaks =c(0,2,10,20), limits = c(0,15)) +theme(axis.title.y=element_blank(), axis.text = element_text(size = 6))
  
ggsave("cluster2_GO.png",plot = MsigDb_GO_2, device = "png",width = 6, height = 2, units = "in")
ggsave("cluster2_GO.pdf",plot = MsigDb_GO_2, device = "pdf",width = 6, height = 2, units = "in",family = "ArialMT")

MsigDb_GO_3 <- ggplot(MsigDb_GO %>% filter(Cluster_subset == "Cluster3") %>% slice(1:10),aes(x = fct_reorder(Gene_Set_Name,log_qvalue), y = log_qvalue)) + 
  geom_bar(fill = heatmap_color[["gene"]][3], stat = "identity", width = 0.5) + coord_flip() + Axis_themes +
  labs(y = "-log10(FDR q-value)") + geom_hline(yintercept = -log10(0.01), linetype = 2) +
  scale_y_continuous(breaks =c(0,2,10,20), limits = c(0,15)) +theme(axis.title.y=element_blank(), axis.text = element_text(size = 6))

ggsave("cluster3_GO.png",plot = MsigDb_GO_3, device = "png",width = 6, height = 2, units = "in")
ggsave("cluster3_GO.pdf",plot = MsigDb_GO_3, device = "pdf",width = 6, height = 2, units = "in",family = "ArialMT")

# Raw TCR Recovery QC --------------------------------------------------------------

#regraph the box plot about the recovery percentage to include all samples
#first read in the necessary csv files saved from the previous analysis

Filtered_coveraged_nonstim <- read_csv("filtered_coverage_unstim.csv")
Filtered_coveraged_stim <- read_csv("filtered_coverage_stim.csv")

#don't know if we'll need to do this, but might as well clean up the data, in case we need to graph paired
#relationship

Filtered_coveraged_nonstim <- Filtered_coveraged_nonstim %>% mutate(Stimulation = FALSE, 
                                                                    orig.ident = factor(orig.ident))
Filtered_coveraged_stim <- Filtered_coveraged_stim %>% mutate(Stimulation = TRUE, orig.ident = factor(orig.ident)) %>%
  mutate(orig.ident = fct_recode(orig.ident, "m3" = "m5","m4" = "Tol.Vac.4"))

#then let's combined them

Filtered_coverage <- rbind(Filtered_coveraged_stim,Filtered_coveraged_nonstim)
Filtered_coverage <- Filtered_coverage %>% mutate(TCR_Recovery = fct_relevel(TCR_Recovery, "TRA","TRB","TRAB"))
#then let's plot the boxplot

Filtered_coverage_plot <- ggplot(Filtered_coverage,aes(TCR_Recovery, Rate*100, color = TCR_Recovery, fill = TCR_Recovery)) + 
  geom_boxplot(fill = NA, width = 0.5) + geom_jitter(pch = 21, color = "black", width = 0.05, height = 0) + scale_fill_manual(values = cbPalette[c(8,3,2)]) +
  scale_color_manual(values = cbPalette[c(8,3,2)]) +
  labs(x = "",y = "Total TCR recovery (% of total cells)", color = "TCR mapping", fill = "TCR mapping") + 
  theme_classic() + Axis_themes + ylim(5,90)

ggsave("Filtered_coverage.pdf",plot = Filtered_coverage_plot, width = 3, height = 2.5, units = "in",family = 'ArialMT')

#Let's try it with the unfiltered data as well

Unfiltered_coveraged_nonstim <- read_csv("unfiltered_coverage_unstim.csv")
Unfiltered_coveraged_stim <- read_csv("unfiltered_coveraged_stim.csv")

#clean up the data as before

Unfiltered_coveraged_nonstim <- Unfiltered_coveraged_nonstim %>% mutate(Stimulation = FALSE, 
                                                                        Sample = factor(Sample))
Unfiltered_coveraged_stim <- Unfiltered_coveraged_stim %>% mutate(Stimulation = TRUE, Sample = factor(Sample)) %>%
  mutate(Sample = fct_recode(Sample, "m3" = "m5","m4" = "Tol.Vac.4"))

Unfiltered_coverage <- rbind(Unfiltered_coveraged_stim,Unfiltered_coveraged_nonstim)

Unfiltered_coverage <- Unfiltered_coverage%>% rename(orig.ident = Sample, TCR_Recovery = variable,
                                                     Rate = value)

Unfiltered_coverage <- Unfiltered_coverage %>% mutate(TCR_Recovery = fct_relevel(TCR_Recovery, "TRA","TRB","TRAB"))

Unfiltered_coverage_plot <- ggplot(Unfiltered_coverage,aes(TCR_Recovery, Rate*100, color = TCR_Recovery, fill = TCR_Recovery)) + 
  geom_boxplot(fill = NA, width = 0.5) + geom_jitter(pch = 21, color = "black", width = 0.05, height = 0) + scale_fill_manual(values = cbPalette[c(8,3,2)]) +
  scale_color_manual(values = cbPalette[c(8,3,2)]) +
  labs(x = "",y = "Total TCR recovery (% of total cells)", color = "TCR mapping", fill = "TCR mapping") + 
  theme_classic() + Axis_themes + ylim(5,90)

ggsave("Unfiltered_coverage.pdf",plot = Unfiltered_coverage_plot, width = 3, height = 2.5, units = "in",family = 'ArialMT')

# Plots for stim vs unstim per cluster ------------------------------------

#I think one of the question people will have is just how much of the cells in each of the large cluster
#consists of the stimulated condition and unstimulated condition.

ggplot(Combined_seurat_Meta, aes(x = res.0.4, fill = Stimulation)) + geom_bar(position = 'fill') + 
  geom_hline(yintercept = 6912/14424, linetype = 2) + labs(x = "Cluster", y = "Fraction of cells in cluster (0-1)") +
  Axis_themes

ggsave('stim_v_unstim_cluster.pdf',width = 3.5, height = 3, units = 'in',family = 'ArialMT')

#then I think plots just by the mouse number.

ggplot(Combined_seurat_Meta, aes(x = res.0.4, fill = orig.ident)) + geom_bar(position = 'fill') + 
  labs(x = "Cluster", y = "Fraction of cells in cluster (0-1)") +
  Axis_themes

ggsave('animals_per_cluster.pdf',width = 3.5, height = 3, units = 'in',family = 'ArialMT')


# Cumulative rank ---------------------------------------------------------
TRB_expan_count <- Combined_seurat_Meta %>% filter(!is.na(TRB_CDR3)) %>% group_by(orig.ident,Stimulation) %>%
  count(TRB_CDR3) %>% arrange(desc(n), .by_group = TRUE) %>% mutate(rank = row_number(desc(n))) %>%
  mutate(proportion = n/sum(n))

TRB_expan_count <- TRB_expan_count %>% mutate(Cumulative_proportion = cumsum(proportion))

#now plot the results
ggplot(TRB_expan_count %>% filter(rank<=20), 
       aes(rank, Cumulative_proportion, color = orig.ident, linetype = Stimulation)) + 
  geom_step(stat = "identity", size = 0.2) + geom_point(color = "black", size = 0.5, pch = 16)+theme_bw() + labs(x = "TRB Clonal Rank",y = "Cumulative Proportion") +
  xlim(1,20) + ylim(0,NA) + Axis_themes

ggsave("Total_clonal_cumulative.pdf",plot = last_plot(),width = 3.5, height = 2.5, units = "in", family = "ArialMT")

#then we want to produce stats that we can mention in the paper.

TRB_expan_count %>% filter(rank<=20) %>% arrange(desc(rank)) %>% slice(1) %>% pull(Cumulative_proportion) %>% summary()

TRB_expan_count %>% filter(rank<=20) %>% arrange(desc(rank)) %>% summarise(total_cells = sum(n)) %>% pull(total_cells) %>%
  sd()

# Public clone analysis ---------------------------------------------------

#First we have to find some way of qunatifying public clones across the 4 animals.
#then somehow, try to get a table of them. I think I may be able to just count, and then filter

TRB_by_animal <- Combined_seurat_Meta %>% count(TRB_CDR3, orig.ident)
TRB_by_animal <- spread(TRB_by_animal, orig.ident, n, fill = 0) %>% mutate(total =rowSums(.[2:5]))
TRB_by_animal <- TRB_by_animal %>% mutate(Animals_number = rowSums(.[2:5]!=0))
At_least_3 <- TRB_by_animal %>% filter(Animals_number>=3)
At_least_3 <- At_least_3 %>% filter(!is.na(TRB_CDR3))

write_csv(At_least_3,"E7_public_clones.csv")

m1_TRB <- Combined_seurat_Meta %>% filter(orig.ident == 'm1' & !is.na(TRB_CDR3)) %>% pull(TRB_CDR3)
m2_TRB <- Combined_seurat_Meta %>% filter(orig.ident == 'm2' & !is.na(TRB_CDR3)) %>% pull(TRB_CDR3)
m3_TRB <- Combined_seurat_Meta %>% filter(orig.ident == 'm3' & !is.na(TRB_CDR3)) %>% pull(TRB_CDR3)
m4_TRB <- Combined_seurat_Meta %>% filter(orig.ident == 'm4' & !is.na(TRB_CDR3)) %>% pull(TRB_CDR3)

library(VennDiagram)

animal_TRB_list <- list('Mouse 1' = m1_TRB,
                        'Mouse 2' = m2_TRB,
                        'Mouse 3' = m3_TRB,
                        'Mouse 4' = m4_TRB)

overlap <- calculate.overlap(x = animal_TRB_list)
venn.diagram(x = animal_TRB_list, filename = 'public_clones_venn.tiff',height = 3, width = 3,
             units = 'in', resolution = 1200, fill = gg_color_hue(4), fontfamily = 'sans')


# Pseudotime analysis -----------------------------------------------------

# An attempt to see if pseudotime analysis would be useful.

require(monocle)
library(grid)
library(VGAM)

monocle_file<-choose.files(caption="Choose scripts that common monocle functions")
source(monocle_file)

keep_cells <- FeaturePlot(Stim_seurat, c("CCR7"), do.identify = TRUE)


Stim_cds <- cdsprocess(SubsetData(Stim_seurat,cells.use = keep_cells, subset.raw = TRUE))

plot_cell_clusters(Stim_cds, color_by = 'TRB_count')

Stim_cds <- learnTrajectory(Stim_cds, 'DDRTree')

plot_cell_trajectory(Stim_cds, color_by = 'TRB_count')

Stim_cds <- orderCells(Stim_cds)

plot_cell_trajectory(Stim_cds, color_by = 'Pseudotime') + theme(legend.position = "right")

Stim_cds_DE <- differentialGeneTest(Stim_cds)

cells_genes <-  apply(Stim_seurat@data > 0, 1, sum) #here it is trying to count genes that are greater than 0
#the resulting array here has a name attribute, to get at the genes where at least 75 cells are expressing.
filtered_genes = names(cells_genes[cells_genes > 75])
length(filtered_genes)
sig = Stim_cds_DE[filtered_genes,]# select for significant genes
sig = sig[sig$qval < 1e-4,]#select for corrected pvalue of less than 0.0001
dim(sig)
sig  = sig[order(sig$qval),]
head(sig, 100)
sum(sig$use_for_ordering)

PT_heatmap <- plot_pseudotime_heatmap(Stim_cds[rownames(sig)[1:100]],
                                      show_rownames = TRUE, norm_method = "log",
                                      hmcols = magma(100),
                                      return_heatmap = TRUE)