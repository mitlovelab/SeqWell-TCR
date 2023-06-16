#Author: Ang Andy Tu
#Purpose: Analysis of the Peanut Allergy samples for Figure 4. Patient 77,74,71,110. Most of the focus of the figure is 
# just on Patient 77.
#Note: This is assuming we've done the regular incorporation of data into the Seurat object, then we only take in
#the samples corresponding to patient 77 cells for the monocle analysis.

# load libraries ----------------------------------------------------------
require(tidyverse)
require(Seurat)
require(monocle)
require(viridis)
library(pheatmap)
library(grid)
library(VGAM)
require(reticulate)

#change the name of the reticulate enivronment to the one matching the packages instruction
#provided on GitHub (list of package/versions for conda).
use_condaenv("r-reticulate-monocle")

# Custom functions --------------------------------------------------------

#load the common seurat and monocle functions

function_file<-choose.files(caption="Choose scripts that contain common functions")
source(function_file)

monocle_file<-choose.files(caption="Choose scripts that common monocle functions")
source(monocle_file)

graphing_parameters <- choose.files(caption="Choose scripts that contain graphing parameters")
source(graphing_parameters)

# Load data-----------
# first we want to load in the PNOIT data, that include the 4 patients (71,74,77,110). They are the ones
# we've included in this study.
PNOIT_seurat <- read_rds("../Hyporeactive_Seurat.rds") # seurat of all patients
P77_seurat <- read_rds("../P77_seurat.rds")# seurat of just Patient77

#then set up the tibble data from the PNOIT dataset
#get the Meta_dataset with barcode updated (this is also generated later in the script. No necessary to laod):
PNOIT_meta <- read_rds("../PNOIT_meta.rds")

#here we make the basic QC plot for the single-cell data. Just for general look.
RNASeq_QC <- VlnPlot(PNOIT_seurat, c("nGene", "nUMI", "reads", "percent.mito"), group.by = "orig.ident", 
        x.lab.rot = TRUE, nCol = 2, size.title.use = 10, point.size.use = 0, do.return = TRUE,return.plotlist = TRUE)

QC_titles <- list("Number of genes","Number of UMI","Number of Reads","Percent of Mitochondrial Genes")

RNASeq_QC <- map2(RNASeq_QC, QC_titles, function(test,QC_titles) test + labs(title = QC_titles, x = "") + 
                    Axis_themes + theme(plot.margin = unit(c(0, 0, 0, 0), "in")))
RNASeq_QC <- map(RNASeq_QC, function(test) test + geom_boxplot(fill = 'white',width = 0.2, outlier.shape = NA))

save_plot("RNASeq_QC.png", plot_grid(plotlist = RNASeq_QC), base_height = 4, unit = "in")
save_plot("RNASeq_QC.pdf", plot_grid(plotlist = RNASeq_QC), base_height = 4, unit = "in", family = "ArialMT")

# Notes: the construction of the trajectory depends on the UMAP construction, and different versions of the umap package
# in testing has shown to result in slightly different trajectory results. In the end, the differences are small: the
# general conclusions stays the same, but the shape of the trajectory could be slightly different.

# Process P77 seurats, filter ---------------------------------
P77_seurat <- subprocess(P77_seurat)

P77_seurat <- TSNE_xy(P77_seurat)
P77_seurat <- UMAP_xy(P77_seurat)
P77_seurat <- Add_TRAB_count(P77_seurat)

CDS_3_seurat <- cdsprocess(P77_seurat)

plot_cell_clusters(CDS_3_seurat, color_by = 'TRB_count')

CDS_3_seurat <- learnTrajectory(CDS_3_seurat, 'DDRTree')

plot_cell_trajectory(CDS_3_seurat, color_by = 'TRB_count')

CDS_3_seurat <- orderCells(CDS_3_seurat)

PT_trajectory_PT <- plot_cell_trajectory(CDS_3_seurat, color_by = 'Pseudotime') + theme(legend.position = "right")

# There are some cells that seem to have RBC contamination. remove them by HBB expression

HBB_cells <- FeaturePlot(P77_seurat, c("HBB"), reduction.use = "umap", do.identify = TRUE)

Final.cell.list <- P77_seurat@cell.names[!(P77_seurat@cell.names %in% HBB_cells)]

Seurat_final <- SubsetData(P77_seurat, cells.use = Final.cell.list, subset.raw = TRUE)

Seurat_final <- subprocess(Seurat_final, perp = 40)

# the clustering here is a bit arbitrary, but using 15 dimensions seems to work well

Seurat_final <- FindClusters(Seurat_final,dims.use=1:15,resolution=1, print.output = 0)
#get rid of the weird ordered factor problem

Seurat_final@ident <- factor(Seurat_final@ident, ordered = FALSE)
DimPlot(Seurat_final, 'umap')

FeaturePlot(Seurat_final,c("ident","IL4","IL5","IL13"), reduction.use = "umap")

Seurat_final <- TSNE_xy(Seurat_final)
Seurat_final <- UMAP_xy(Seurat_final)
Seurat_final <- Add_TRAB_count(Seurat_final)

Final_meta_data <- Seurat_final@meta.data


# Psuedotime analysis  -----------------------------------------------

##note: it's really important to re-load the conda-env command before running anything.

CDS_final <- cdsprocess(Seurat_final)

plot_cell_clusters(CDS_final, color_by = 'TRB_count')

CDS_final <- learnTrajectory(CDS_final, 'DDRTree')

plot_cell_trajectory(CDS_final, color_by = 'TRB_count')

CDS_final <- orderCells(CDS_final)

PT_trajectory_PT <- plot_cell_trajectory(CDS_final, color_by = 'Pseudotime', cell_size = 0.5) + theme(legend.position = "none") +
                    Axis_themes + scale_x_reverse()

# now that we can reverse the scaling for plotting, let's see if we can do the same trajectory mapped with TRB_data

PT_trajectory <- plot_cell_trajectory(CDS_final, color_by = 'TRB_count') + scale_x_reverse()
segment_trajectory <- PT_trajectory$layers[[1]]
PT_data <- PT_trajectory$data %>% as.tibble()
PT_data <- PT_data %>% dplyr::rename(Component1 = data_dim_1, Component2 = data_dim_2)
PT_data <- PT_data %>% arrange(!is.na(TRB_count),TRB_count) #this makes the TRB_count show NA first

TRB_count_PT <- ggplot(PT_data, aes(x = Component1, y = Component2)) + segment_trajectory + 
  geom_point(aes(x = Component1, y = Component2, color = TRB_count), size = 0.5) + 
  scale_color_viridis(direction = -1, name = "Pseudotime", na.value = "grey87") + Axis_themes + theme(legend.position = "none")+
  scale_x_reverse()

plot_grid(PT_trajectory_PT,TRB_count_PT,ncol = 1) #500x650

ggsave("PT_image_nolegend.eps",plot = last_plot(), device = "eps",width = 3, height = 3, units = "in", family = "ArialMT")
  
#remember to load VGAM. otherwise DE won't work.
CDS_final_DE <- differentialGeneTest(CDS_final)

#the first stuff is just counting how many cells are expressing each gene
cells_genes <-  apply(Seurat_final@data > 0, 1, sum) #here it is trying to count genes that are greater than 0
#the resulting array here has a name attribute, to get at the genes where at least 75 cells are expressing.
filtered_genes = names(cells_genes[cells_genes > 75])
length(filtered_genes)
sig = CDS_final_DE[filtered_genes,]# select for significant genes
sig = sig[sig$qval < 1e-4,]#select for corrected pvalue of less than 0.0001
dim(sig)
sig  = sig[order(sig$qval),]
head(sig, 50)
sum(sig$use_for_ordering)


#this is done with log normalization, then scaling the max and minimum to 3 and -3. default k = 6
PT_heatmap <- plot_pseudotime_heatmap(CDS_3_seurat[rownames(sig)[1:100]],
                                      show_rownames = TRUE, norm_method = "log",
                                      hmcols = magma(100),
                                      return_heatmap = TRUE)

save_pheatmap(PT_heatmap, "pseudotime_gene_heatmap.png",width = 4.5, height = 8)
save_pheatmap(PT_heatmap, "pseudotime_gene_heatmap.pdf",device = "pdf",width = 4.5, height = 8)


#Now we want to get the genes in each cluster
PT_heatmap_cluster_list <- cutree(PT_heatmap$tree_row, k =6)

Cluster_gene_list_k6<- tibble(genes = names(PT_heatmap_cluster_list), cluster = PT_heatmap_cluster_list)
write_csv(Cluster_gene_list_k6,"PT_heatmap_genes_k6.csv")

#let's try with k = 5

PT_heatmap_5 <- plot_pseudotime_heatmap(CDS_3_seurat[rownames(sig)[1:100]],
                                      show_rownames = TRUE, norm_method = "log",
                                      hmcols = magma(100),
                                      return_heatmap = TRUE, num_clusters = 5)

save_pheatmap(PT_heatmap_5, "pseudotime_gene_heatmap.png",width = 4.5, height = 8)
save_pheatmap(PT_heatmap_5, "pseudotime_gene_heatmap_k5.pdf",device = "pdf",width = 4.5, height = 8)

PT_heatmap_cluster_list_k5 <- cutree(PT_heatmap_5$tree_row, k =5)
Cluster_gene_list_k5<- tibble(genes = names(PT_heatmap_cluster_list_k5), cluster = PT_heatmap_cluster_list_k5)

write_csv(Cluster_gene_list_k5,"PT_heatmap_genes_k5.csv")

#to see if our gene enrichment is robust, we should expand the number of genes we include,
#just so we can see if we can get the same msig db enrichment

PT_heatmap_5_200genes <- plot_pseudotime_heatmap(CDS_3_seurat[rownames(sig)[1:200]],
                                        show_rownames = TRUE, norm_method = "log",
                                        hmcols = magma(100),
                                        return_heatmap = TRUE, num_clusters = 5)

save_pheatmap(PT_heatmap_5_200genes, "pseudotime_200gene_heatmap_k5.png",width = 4.5, height = 8)

PT_heatmap_cluster_list_k5_200genes <- cutree(PT_heatmap_5_200genes$tree_row, k =5)
Cluster_200gene_list_k5<- tibble(genes = names(PT_heatmap_cluster_list_k5_200genes), cluster = PT_heatmap_cluster_list_k5_200genes)

write_csv(Cluster_200gene_list_k5,"PT_heatmap_200genes_k5.csv")

#there are the num_cluster which will determine how many clusters the heatmap will try to draw. Need better
#thankfully, the labels are already at 6pt. so no need to change them for that reason.

# Incorporating TRB clonality ---------------------------------------------


#I think a better way to graph where the clones are falling on the pseudotime spectrum is to do box plot
#for each clones.

PT_data_TRB_box <- PT_data %>% filter(TRB_count >=5)
#sort the levels and whatnot so it graphs out correctly
PT_data_TRB_box <- PT_data_TRB_box %>% mutate(TRB_CDR3 = fct_drop(TRB_CDR3)) %>% 
                  mutate(TRB_CDR3 = fct_reorder(TRB_CDR3,Pseudotime, .desc = TRUE))
#axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)

TRB_CDR3_PT_box <- ggplot(PT_data_TRB_box, aes(x = TRB_CDR3, y = Pseudotime, fill = Pseudotime)) + geom_boxplot(outlier.shape = NA) + 
                    geom_jitter(width = 0.3, pch = 21,color = "black") + scale_fill_viridis(option = "plasma")  +
                    theme_minimal() + theme(axis.title.y=element_blank(),axis.title.x= element_blank()) +
                    coord_flip() + Axis_themes

ggsave("TRB_CDR3_PT_box.png",plot = TRB_CDR3_PT_box, device = "png",width = 5, height = 3.5, units = "in")
ggsave("TRB_CDR3_PT_box.eps",plot = TRB_CDR3_PT_box, device = "eps",width = 5, height = 3.5, units = "in")
##Then I think, let's try to graph a couple of the genes on Feature Plot, also get the code down for the styling
#of the TSNE/UMAP Plots

genes_to_plot <- c("SELL", "CCR7", "TCF7", "LEF1","CD69","TNF","IL2RA","GATA3","CD40LG")

genes_to_plot_rev <- c("SELL", "CCR7", "LEF1","CD69","CD40LG","TNFRSF4","GATA3","IL5","IL13")

plot_genes_in_pseudotime(CDS_final[c("FOS","NFKBIA","JUN","GATA3","IL2RA")], color_by = 'Pseudotime') + Axis_themes
ggsave("Selected_genes_PT.eps",plot = last_plot(), device = "eps",width = 3.5, height = 3.5, units = "in", family = "ArialMT")

plot_genes_in_pseudotime(CDS_final[c("IL5","IL9","IL13","IL17RB")], color_by = 'Pseudotime') + Axis_themes
ggsave("Th2_genes_PT.eps",plot = last_plot(), device = "eps",width = 3, height = 3.5, units = "in", family = "ArialMT")

#I think we can say that the pseudotime correspond with T cell activation, in that higher time should mean higher
#response. But I think it's hard to say that maybe some cells are on their way to full activation.

#But, to say that they actually have different mean, or averages between these groups, let's show some sort
#of ANOVA.

aovMod_Selected_CDR3 <- aov(Pseudotime ~ TRB_CDR3, data=PT_data_TRB_box)
summary(aovMod_Selected_CDR3)

#ANOVA actually shows that they are not significantly different, which makes sense since we did
#select for clones that are towards the end of psuedotime by selecting for expansion. 

#try Levene's Test
library(car)
leven_test <- leveneTest(Pseudotime ~ TRB_CDR3, data=PT_data_TRB_box)
#this is also not significant. I think we're just seeing a lot of TCR groups that are highly correlated
#which is somewhat fine, for this paper. I'm willing to bet we can further divide group, if we want to show
#a specific clone has different distribution (might just use T-test if we end up using paired test)

# Plot specific genes overlaid on UMAP/Tsne ---------------------------------
Umap_grid <- FeaturePlot(Seurat_final, genes_to_plot,reduction.use = "umap", do.return = TRUE, pt.size = 0.5)


Umap_grid <- lapply(Umap_grid, `+`, TSNE_theme)
#to plot these, we'd have to use plot_grid(plotlist = ...)
plot_grid(plotlist = Umap_grid) #saved as a square 288x288 (I think 320x320 would work)

#then we have to plot the expansion of the Tcell on the Umap as well.

UMAP_TRB_expan <- ggplot(arrange(Final_meta_data,TRB_count), aes(x = UMAP.x, y = UMAP.y)) + 
  geom_point(color = "grey87", size = 1) + geom_point(data = arrange(Final_meta_data,TRB_count) %>% filter(TRB_count>0), 
                                                      aes(x = UMAP.x, y = UMAP.y, color = TRB_count), size = 1) + 
  scale_color_viridis(direction = -1, name = "Clonal size\n(TCRB)") + 
  labs(title = "Patient 77", x = "UMAP 1", y = "UMAP 2") + Axis_themes

dev.off()

ggsave("UMAP_TRB_expan.eps",plot = UMAP_TRB_expan, device = "eps",width = 4, height = 3.5, units = "in", family = "ArialMT")
ggsave("UMAP_TRB_expan.png",plot = UMAP_TRB_expan, device = "png",width = 4, height = 3.5, units = "in")

#One of the reviewer comment is that we need to chagne it back to Tsne, since UMAP is confusing people.

#start with the grid of genes:
tSNE_grid <- FeaturePlot(Seurat_final, genes_to_plot_rev,reduction.use = "tsne", do.return = TRUE, pt.size = 0.3) 
tSNE_grid <- map(tSNE_grid, `+`, TSNE_theme)
tSNE_grid <- map(tSNE_grid, `+`, theme(legend.position = 'right',legend.title = element_blank(),
                                       legend.key.size = unit(0.2,'in'),
                                       legend.key.height = unit(0.8,'in')))

plot_grid(plotlist = tSNE_grid)

save_plot("tSNE_grid.png", plot_grid(plotlist = tSNE_grid), base_width = 3.2, base_height = 3.2, unit = "in")
save_plot("tSNE_grid.pdf", plot_grid(plotlist = tSNE_grid), base_width = 3.2, base_height = 3.2, unit = "in", family = "ArialMT")

#we have to re-plot this to get the plots to behave interms of constant scale, and for humna data, having the
#correct scale is much more important.
Naive_genes <- genes_to_plot_rev[1:3]
Activation_genes <- genes_to_plot_rev[4:6]
Th2_tsne_genes <- genes_to_plot_rev[7:9]

gene_list <- list(Naive_genes,Activation_genes,Th2_tsne_genes)
names(gene_list) <- c("Naive","Activation","Th2")

Selected_T_cell_genes <- map(gene_list,function(gene_list) FeaturePlot(Seurat_final, gene_list, pt.size = 0.2, do.return = TRUE))
Selected_T_cell_genes <- map(Selected_T_cell_genes,function(Selected_T_cell_genes) map(Selected_T_cell_genes,`+`,TSNE_theme))

Selected_T_cell_genes <- map(Selected_T_cell_genes,function(Selected_T_cell_genes) map(Selected_T_cell_genes,`+`,theme(legend.position = 'right', legend.text = element_text(size = 6),
                                                                                                                            legend.title = element_blank(), legend.key.size = unit(0.05,'in'),
                                                                                                                            legend.key.height = unit(0.1,'in'),
                                                                                                                       legend.margin=margin(t = 0, unit='cm'),
                                                                                                                       plot.margin = unit(c(0, 0, 0, 0), "cm"),
                                                                                                                       legend.background = element_blank())))
#let's try just not scaling them

save_plot("Naive_genes_scaled.pdf", plot_grid(plotlist = Selected_T_cell_genes[["Naive"]],nrow = 1),ncol = 3, base_height = 1.0, base_aspect_ratio = 1.1, unit = "in", family = "ArialMT")
save_plot("Act_genes_scaled.pdf", plot_grid(plotlist = Selected_T_cell_genes[["Activation"]],nrow = 1),ncol = 3, base_height = 1.0, base_aspect_ratio = 1.1, unit = "in", family = "ArialMT")
save_plot("Th2_genes_scaled.pdf", plot_grid(plotlist = Selected_T_cell_genes[["Th2"]],nrow = 1),ncol = 3, base_height = 1.0, base_aspect_ratio = 1.1, unit = "in", family = "ArialMT")


#below is if we're scaling them.
scale_naive = c(0,3.5)
scale_Eff_em = c(0,4)
scale_th2 = c(0,6.5)
color_scale <- list(scale_naive, scale_Eff_em, scale_th2)

Selected_T_cell_genes <- map2(Selected_T_cell_genes,color_scale,
             function(Selected_T_cell_genes,color_scale) map(Selected_T_cell_genes,`+`,
                                                                  scale_color_gradient(low = "yellow", high = "red",space = "Lab", na.value = "red", 
                                                                                       guide = "colourbar",aesthetics = "colour", limits = color_scale)))

save_plot("Naive_genes_legend.pdf", plot_grid(plotlist = Selected_T_cell_genes[["Naive"]],nrow = 1),ncol = 3, base_height = 1.2, base_aspect_ratio = 0.9, unit = "in", family = "ArialMT")
save_plot("Activation_genes_legend.pdf", plot_grid(plotlist = Selected_T_cell_genes[["Activation"]],nrow = 1),ncol = 3, base_height = 1.2, base_aspect_ratio = 0.9, unit = "in", family = "ArialMT")
save_plot("Th2_genes_legend.pdf", plot_grid(plotlist = Selected_T_cell_genes[["Th2"]],nrow = 1),ncol = 3, base_height = 1.2, base_aspect_ratio = 0.9, unit = "in", family = "ArialMT")

Selected_T_cell_genes_nolegend <- map(Selected_T_cell_genes,function(Selected_T_cell_genes) map(Selected_T_cell_genes,`+`,theme(legend.position = 'none')))

save_plot("Naive_genes_scaled.pdf", plot_grid(plotlist = Selected_T_cell_genes_nolegend[["Naive"]],nrow = 1),ncol = 3, base_height = 1.1, base_aspect_ratio = 1, unit = "in", family = "ArialMT")
save_plot("Eff_Em_genes_scaled.pdf", plot_grid(plotlist = Selected_T_cell_genes_nolegend[["Activation"]],nrow = 1),ncol = 3, base_height = 1.1, base_aspect_ratio = 1, unit = "in", family = "ArialMT")
save_plot("Act_exh_genes_scaled.pdf", plot_grid(plotlist = Selected_T_cell_genes_nolegend[["Th2"]],nrow = 1),ncol = 3, base_height = 1.1, base_aspect_ratio = 1, unit = "in", family = "ArialMT")

#TCR clone
tSNE_TRB_expan <- ggplot(arrange(Final_meta_data,TRB_count), aes(x = Tsne.x, y = Tsne.y)) + 
  geom_point(color = "grey87", size = 1) + geom_point(data = arrange(Final_meta_data,TRB_count) %>% filter(TRB_count>0), 
                                                      aes(x = Tsne.x, y = Tsne.y, color = TRB_count), size = 1) + 
  scale_color_viridis(direction = -1, name = "Clonal size\n(TCRB)") + 
  labs(title = "Patient 77", x = "tSNE 1", y = "tSNE 2") + Axis_themes

ggsave("tSNE_TRB_expan.eps",plot = tSNE_TRB_expan, device = "eps",width = 4, height = 3.5, units = "in", family = "ArialMT")
ggsave("tSNE_TRB_expan.png",plot = tSNE_TRB_expan, device = "png",width = 4, height = 3.5, units = "in")


#Let's just find the markers that are differentiating the expanded cells, just see if this kind of analysis would be
#helpful

FindMarkers(object = Seurat_final,ident.1 = c(4,2,1),test.use = "bimod") -> DE_genes_expanded

#the markers that are most associated with the expanded cells are the ones we already suspected. 


#I think just to test, try to find dot plots of these genes for these particular clones. 
#we can do the work with Seurat_final. First let's just get the cells for which we actually have TCRB for

Final_meta_data <- as.tibble(Final_meta_data)
TRB_cell_list <- Final_meta_data %>% filter(!is.na(TRB_CDR3)) %>% select(cell_id)

Seurat_final_TRB <- SubsetData(Seurat_final, cells.use = TRB_cell_list$cell_id, subset.raw = TRUE)
Seurat_final_TRB <- SetAllIdent(Seurat_final_TRB, id = "TRB_CDR3")

TRB_expanded_ident <- levels(PT_data_TRB_box$TRB_CDR3)
Seurat_final_TRB_expanded <- SubsetData(Seurat_final_TRB, ident.use = TRB_expanded_ident, subset.raw = TRUE)

Th2_genes <- c("IL5","IL9","IL13","IL17RB")
Plasma_10 <- plasma(10)

Seurat_final_TRB_expanded@ident <- factor(Seurat_final_TRB_expanded@ident, levels = TRB_expanded_ident)
TRB_TH2_genes <- DotPlot(Seurat_final_TRB_expanded, genes.plot = Th2_genes, cols.use = c(Plasma_10[1],Plasma_10[10]),plot.legend = TRUE, do.return = TRUE) + Axis_themes

TRB_TH2_genes <- TRB_TH2_genes + theme_minimal() + Axis_themes + theme(legend.position = "left",
                                                                       axis.title.x= element_blank(),
                                                                       axis.title.y= element_blank())

ggsave("TRB_TH2_genes.png",plot = TRB_TH2_genes, device = "png",width = 4, height = 3.5, units = "in")
ggsave("TRB_TH2_genes.eps",plot = TRB_TH2_genes, device = "eps",width = 4, height = 3.5, units = "in")
#I think the other thing is to explain that the most expanded TCR cells are not always the most active ones. 

TCR_activation_genes_PT<- plot_genes_in_pseudotime(CDS_final[c("FOS","NFKBIA","JUN","GATA3","CD40LG", "IL5","IL9","IL13","IL17RB")], color_by = 'Pseudotime', cell_size = 0.2,
                         panel_order = c("FOS","JUN","NFKBIA","GATA3","CD40LG", "IL5","IL9","IL13","IL17RB")) + Axis_themes +theme(legend.position = "right",strip.text.x = element_text(size = 8))
TCR_activation_genes_PT

# MsigDb results ----------------------------------------------------------

# activation signature from the analysis above was taken to msigdb online portal, and queried against available 
# signatures. the releveant results are then loaded it here for visualization

MsigDb <- read_tsv("k5_200genes_MsigDB.txt") # 
MsigDb <- MsigDb[1:40,]
MsigDb <- MsigDb %>% mutate(Cluster_subset = factor(Cluster_name),Genes_Gene_Set_K = as.integer(Genes_Gene_Set_K))

MsigDb <- MsigDb %>% mutate(log_qvalue = -log10(FDR_q_value))

gg_color <- gg_color_hue(5)

MsigDb <- MsigDb %>% mutate(Cluster_subset = fct_recode(Cluster_subset,"Early" = "Early_PT","Late" = "Late_PT"))
MsigDb <- MsigDb %>% mutate(Short_Description = str_remove_all(Short_Description, "Up-regulated in "))
MsigDb <- MsigDb %>% mutate(Short_Description = firstup(Short_Description))

ggplot(MsigDb %>% filter(!is.na(Short_Description)),aes(x = fct_reorder(Systematic_name,log_qvalue), y = log_qvalue, fill = Cluster_subset)) + 
  geom_bar(stat = "identity") + coord_flip() + geom_text(aes(label = Short_Description), y = 2.1, hjust = 0, size = 2) + Axis_themes +
  labs(y = "-log10(FDR q-value)", x = "Gene set systematic ID", fill = "Psuedotemporal Enrichment") + geom_hline(yintercept = -log10(0.01), linetype = 2) +
  scale_y_continuous(breaks =c(0,2,10,20,30)) + theme(legend.position = "bottom")

ggsave("PT_MsigDb.png",plot = last_plot(), device = "png",width = 3, height = 3.5, units = "in")
ggsave("PT_MsigDb.eps",plot = last_plot(), device = "eps",width = 3, height = 3.5, units = "in", family = "ArialMT")

firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}
MsigDb %>% mutate(Short_Description = firstup(Short_Description))

# Linear Regression -------------------------------------------------------

ggplot(PT_data, aes(x = TRB_count, y = Pseudotime)) + geom_point()
PT_data_cleaned <- PT_data %>% filter(!is.na(TRB_count))


linearMod <- lm(Pseudotime ~ TRB_count, data=PT_data %>% filter(!is.na(TRB_count)))
summary(linearMod)
#linear model is not a good fit

#I think anova is best suited if we want to test whether the means of the bins
#are significantly different.
aovMod <- aov(TRB_count ~ Pseudotime, data=PT_data_cleaned)
summary(aovMod)

#Spearman and Kendall

Spear <- cor.test(PT_data_cleaned$TRB_count, PT_data_cleaned$Pseudotime, method = "spearman", alternative = "two.sided")
Kendall <- cor.test(PT_data_cleaned$TRB_count, PT_data_cleaned$Pseudotime, method = "kendall")

# Investigating what separates the patients ------------------------------

#Let's see why the the patients are not really overlapping on our umap/tsne representations.

# PNOIT_seurat contains all the cells from the hyporeactive patients, which is relevant ones to this study

PNOIT_seurat <- SetAllIdent(PNOIT_seurat, id = 'orig.ident')

levels(PNOIT_seurat@ident)

DE_P77 <- DE_test(PNOIT_seurat, group_1 = "Patient77", logfc.threshold = 0.1)
DE_P71 <- DE_test(PNOIT_seurat, group_1 = "Patient71", logfc.threshold = 0.1)
DE_P74 <- DE_test(PNOIT_seurat, group_1 = "Patient74", logfc.threshold = 0.1)
DE_P110 <- DE_test(PNOIT_seurat, group_1 = "Patient110", logfc.threshold = 0.1)

#it does seem like mostof the differentiating genes are ones that are metabolic, except
#for some th2 genes, which we know the P77 just has more TH2 phenotypic cells, 

patient_genes <- c("LHCGR","NACA2","MTRNR2L8",
                   "SNORD89","KIAA0319L","GBP5","RPS4Y1","XIST")

patient_genes_plots <- FeaturePlot(PNOIT_seurat, patient_genes,pt.size = 0.5,do.return = TRUE)

patient_genes_plots <- map(patient_genes_plots,`+`, TSNE_theme)

patient_genes_vln <- VlnPlot(PNOIT_seurat,patient_genes, point.size.use = 0, do.return = TRUE, 
                             group.by = 'orig.ident', size.x.use = 6, size.y.use = 6, size.title.use = 8,
                             x.lab.rot = TRUE, remove.legend = TRUE, return.plotlist = TRUE)

patient_genes_vln <- map(patient_genes_vln,`+`, Axis_themes)
patient_genes_vln <- map(patient_genes_vln,`+`,geom_boxplot(fill = 'white',width = 0.1,outlier.shape = NA))
patient_genes_vln <- map(patient_genes_vln,`+`,theme(axis.title.x = element_blank()))
plot_grid(plotlist = patient_genes_plots)
plot_grid(plotlist = patient_genes_vln, ncol = 4)
save_plot("patient_genes_plots.png", plot_grid(plotlist = patient_genes_plots), base_width = 3.2, base_height = 3.2, unit = "in")
save_plot("patient_genes_plots.eps", plot_grid(plotlist = patient_genes_plots), base_width = 3.2, base_height = 3.2, unit = "in", family = "ArialMT")

save_plot("patient_genes_vln.png", plot_grid(plotlist = patient_genes_vln, ncol = 4), base_width = 8, base_height = 3.2, unit = "in")
save_plot("patient_genes_vln.pdf", plot_grid(plotlist = patient_genes_vln, ncol = 4), base_width = 8, base_height = 3.2, unit = "in", family = "ArialMT")

# Investigate modules from literature ---------------------------------------

#we need to map some known modules, and see how they work.
module_lists <- read_csv("modules_CD4.csv")
module_lists <- as.list(module_lists[,1:14])

module_lists <- map(module_lists,function(list) list[!is.na(list)])

All_genes <- rownames(Seurat_final@data)

module_lists_fiiltered <- map(module_lists, function(module_lists) intersect(module_lists, All_genes))

Seurat_final <- AddModuleScore(Seurat_final,genes.list = module_lists_fiiltered, ctrl.size = 50, enrich.name = "Cluster")

#eh mapping the modules don't look that great to begin with. Module 7, which is naive, maps okay.
#Th2 map only has a slight increase in mapping, the signals are just not that strong

Plot_titles <- names(module_lists_fiiltered)
module_grid_1 <- FeaturePlot(Seurat_final, str_c("Cluster",1:7),reduction.use = "tsne", do.return = TRUE, pt.size = 0.3) 
module_grid_1 <- map2(module_grid_1, Plot_titles[1:7], function(module_grid_1,Plot_titles) module_grid_1 + TSNE_theme + labs(title = Plot_titles) + theme(plot.title = element_text(size = 6)))

module_grid_2 <- FeaturePlot(Seurat_final, str_c("Cluster",8:14),reduction.use = "tsne", do.return = TRUE, pt.size = 0.3) 
module_grid_2 <- map2(module_grid_2, Plot_titles[8:14], function(module_grid_2,Plot_titles) module_grid_2 + TSNE_theme + labs(title = Plot_titles)+ theme(plot.title = element_text(size = 6)))



save_plot("module_grid_1.pdf", plot_grid(plotlist = module_grid_1), base_width = 3.2, base_height = 3.2, unit = "in", family = "ArialMT")
save_plot("module_grid_2.pdf", plot_grid(plotlist = module_grid_2), base_width = 3.2, base_height = 3.2, unit = "in", family = "ArialMT")

# plot the shorten modules with wei-et-al. signatures. This is taken from Figure 3 of the 2009 paper in Immunity.

module_short <- read_csv("modules_CD4_short.csv")
module_short <- as.list(module_short[,1:5])
module_short <- map(module_short,function(list) list[!is.na(list)])

module_short_fiiltered <- map(module_short, function(module_short) intersect(module_short, All_genes))

Seurat_final <- AddModuleScore(Seurat_final,genes.list = module_short_fiiltered, ctrl.size = 50, enrich.name = "Short")
short_titles <- names(module_short_fiiltered)

module_grid_short <- FeaturePlot(Seurat_final, str_c("Short",1:5),reduction.use = "tsne", do.return = TRUE, pt.size = 0.3) 
module_grid_short <- map2(module_grid_short, short_titles, function(module_grid_short,short_titles) module_grid_short + TSNE_theme + 
                            labs(title = short_titles) + theme(plot.title = element_text(size = 6)))

save_plot("Wei_09_grid.pdf", plot_grid(plotlist = module_grid_short), base_width = 3.2, base_height = (3.2/3)*2, unit = "in", family = "ArialMT")

module_grid_short <- map(module_grid_short,`+`,theme(legend.position = 'right', legend.text = element_text(size = 6),
                                                                                                                       legend.title = element_blank(), legend.key.size = unit(0.05,'in'),
                                                                                                                       legend.key.height = unit(0.1,'in'),
                                                                                                                       legend.margin=margin(t = 0, unit='cm'),
                                                                                                                       plot.margin = unit(c(0, 0, 0, 0), "cm"),
                                                                                                                       legend.background = element_blank()))

save_plot("Wei_09_grid_legend.pdf", plot_grid(plotlist = module_grid_short), base_width = 3.2, base_height = (3.2/3)*2, unit = "in", family = "ArialMT")

# Plot specific clones on tsne --------------------------------------------
TRB_expanded_ident

plot_TR(Final_meta_data, clonotype = TRB_expanded_ident[6])