#Author: Ang Andy Tu
#Purpose: Analysis of the E7 immunized mice in Figure 3.
#Note: To self: an older version of the code is stored off github. That includes intermeidate test steps that 
#ultimately was not needed to generate the figures.

# Load libraries and functions --------------------------------------------

library(tidyverse)
library(viridis)
library(Seurat)
library(pheatmap)
require(reticulate)
use_condaenv("r-reticulate-duncan")

## Custom Functions --------------------------------------------------------
function_file<-choose.files(caption="Choose scripts that contain common functions")
source(function_file)

graphing_parameters <- choose.files(caption="Choose scripts that contain graphing parameters")
source(graphing_parameters)

# custom load files -------------------------------------------------------

Tol.Vac<-choose.files(caption="Tol.vac TCR summary file")

Tol.Vac_summary <- read_tsv(Tol.Vac)

TOl.Vac.4_BC <- Meta_data %>% filter(orig.ident == "Tol.Vac.4") %>% select(BC)

Tol.Vac.4_summary <- Tol.Vac_summary %>% filter(BC %in% TOl.Vac.4_BC$BC)

#then let's read in the summary from M1 and M2.

file <- choose.files(caption="M1 + M2 TCR summary file")

M1_tet.s <- read_tsv(file[1])
M2_tet.s <- read_tsv(file[2])

M1_M2_Tol.Vac.4_summary <- bind_rows(M1_tet.s,M2_tet.s, Tol.Vac.4_summary)
#when M5 data comes in, we can just add that to this file
TCR_df <- M1_M2_Tol.Vac.4_summary

file <- choose.files(caption="M5 TCR summary file")

M5_tet.s <- read_tsv(file)

tet.s_summary <- bind_rows(M1_M2_Tol.Vac.4_summary, M5_tet.s)

TCR_df <- tet.s_summary 

#the code seemed to stop at TRBV stage. Something werid happen at the beginning, where variables were filtered weird

#it seems that the M2 library are of lower quality. Well let's keep it, maybe when M5 comes in, we'll throw out M2, if it
#seems like it's also a lot of things that we end up clearing out during the filtering step.

#I think the strategy now is still to just focus on the three mice that we have data for, then add in the data as we go.

#I think there was a bit of issue with summary, since I think starting with Tibble would change the the way that
#the data actaully is summarized. 

#I think with the new data coming, let's just process the whole tsne, with the same dimension reduction that we 
#were looking for, then hopefully add in data as they get processed.

#Ran through all the mice, as it turns out, TOl.Vac.4 (earliest one) was actually the worst recoery and 
#Now, I think we want to check to see if the clonality would make sense
#though important to remember, we might want to somehow account for sample to sample differences

# Initial analysis of full set --------------------------------------------


umi_seurat <- TSNE_xy(umi_seurat)

umi_seurat <- Add_TRAB_count(umi_seurat)

Meta_data <- as.tibble(umi_seurat@meta.data)

#spot check to see the relationhsip between expansion, and tsne
TR_expan_map(Meta_data, ft.size = 0.5)

TRB_expan_by_animal<- Meta_data %>% group_by(orig.ident) %>% count(TRB_CDR3) %>% rename(TRB_count_by_animal = n)

Meta_data <- left_join(Meta_data,TRB_expan_by_animal, by = c("orig.ident","TRB_CDR3"))

ggplot(arrange(Meta_data,TRB_count), aes(x = Tsne.x, y = Tsne.y)) + 
  geom_point(color = "grey87", size = 1) + 
  geom_point(data = arrange(Meta_data,TRB_count_by_animal) %>% filter(TRB_count>0), 
                                                      aes(x = Tsne.x, y = Tsne.y, color = log2(TRB_count_by_animal)), size = 1) +
  scale_color_viridis(direction = -1, na.value = "grey87") + 
  facet_wrap(~orig.ident)

#The results seems to show that there may still be some differences, but some of it is definitely driven by
#animals.

# Clean Seurat, and limit variable genes ----------------------------------


#first let's clean out the B cells, and cells of low quality, that came out as a separate cluster.

E7_cleaned_seurat <- SubsetData(umi_seurat, ident.remove = c("13","15"), subset.raw = TRUE)

E7_cleaned_seurat <- NormalizeData(object = E7_cleaned_seurat)
#Then we find out the length of the variable genes we found here, and set the number of variable genes we want

length_var <- length(E7_cleaned_seurat@var.genes)
n_var_genes <- 1000
cf <- 0.2

#in FindVariableGenes, y axis is the log(variance/mean), while x axis is mean expression. 
#let's just write a function for this, since we'll likely need to repeat this many times

limit_vargenes(seurat = E7_cleaned_seurat)

E7_lim_1k <- limit_vargenes(seurat = E7_cleaned_seurat)

E7_lim_1k <- ScaleData(object = E7_lim_1k, vars.to.regress = c("percent.mito","nUMI"), genes.use = E7_lim_1k@var.genes, model.use = "poisson")

E7_lim_1k <- RunPCA(object = E7_lim_1k, pc.genes = E7_lim_1k@var.genes, do.print = FALSE, pcs.compute = 40, maxit = 500, weight.by.var = FALSE)

#now that we have the PCA ran on the data, we can try to find a cutoff, then look for 
s <- E7_lim_1k@dr$pca@sdev


#This part is doing some geometry calculation. we should write this into a function.
allCoor<-cbind(1:length(s),s)
lineVec<-allCoor[length(s),]-allCoor[1,]
## normalize the line vector
lineVecN = lineVec / sqrt(sum(lineVec^2));
## find the distance from each point to the line:
## vector between all points and first point
vecFromFirst<-allCoor-do.call("rbind", rep(list(allCoor[1,]), length(s)))
q<-do.call("rbind", rep(list(lineVecN), length(s)))
scalarProduct<-q[,1]
for (i in 1:length(s)) {
  scalarProduct[i]<- vecFromFirst[i,] %*% q[i,]
}
vecFromFirstParallel = scalarProduct * q
vecToLine = vecFromFirst - vecFromFirstParallel
distToLine = sqrt(rowSums(vecToLine^2))
##Point in elbow is point furthest from line connecting first and last point
pcCut<-which.max(distToLine)

E7_lim_1k <- RunTSNE(object = E7_lim_1k, dims.use = 1:pcCut, do.fast = TRUE)

E7_lim_1k <- RunUMAP(object = E7_lim_1k, dims.use = 1:pcCut)

#I think it's keeping some of the structure that we are seeing before, but some new 
E7_lim_1k <- TSNE_xy(E7_lim_1k)
E7_lim_1k <- UMAP_xy(E7_lim_1k)
E7_lim_1k <- Add_TRAB_count(seurat = E7_lim_1k)

# Start analysis with 1k variable genes -----------------------------------


Meta_data_1k <- E7_lim_1k@meta.data

TR_expan_map(Meta_data_1k) + facet_wrap(~orig.ident)

TRB_expan_by_animal<- Meta_data_1k %>% group_by(orig.ident) %>% count(TRB_CDR3) %>% rename(TRB_count_by_animal = n) %>% filter(!is.na(TRB_CDR3))

Meta_data_1k <- left_join(Meta_data_1k,TRB_expan_by_animal, by = c("orig.ident","TRB_CDR3"))

TR_expan_map(Meta_data_1k,TRB_count_by_animal) + facet_wrap(~orig.ident)




#this is good to know that the clonal differences that we see are still there. I think we need to run the rest
#of the standard seurat analysis, and see what we get.

# re-run basic Seurat output ----------------------------------------------


direct<-dirname(file)
base<-basename(file)
setwd(direct)
if (file.exists("Plots")){
} else {
  dir.create("Plots")
}

pdf("Plots/tSNEClusterPlot_Samples.pdf")
TSNEPlot(E7_lim_1k, group.by = "orig.ident")
dev.off()

E7_lim_1k<-FindClusters(E7_lim_1k,dims.use=1:pcCut,resolution=res, print.output = 0,temp.file.location=paste(direct,"/",sep=""))

pdf("Plots/tSNEClusterPlot_Clusters.pdf")
TSNEPlot(E7_lim_1k)
dev.off()

pdf("Plots/tSNEClusterPlot_ClustersLabeled.pdf")
TSNEPlot(E7_lim_1k,pt.size = 1.5,do.label=TRUE)
dev.off()

markers_all <- FindAllMarkers(E7_lim_1k,test.use = "roc", do.print = TRUE)
write_csv(markers_all,"Cluster_Enrichments.txt")

markers_use <- markers_all %>% filter(avg_diff>0&power>0.3) %>% select(gene) %>% unique()

write_csv(markers_use,"Cluster_Unique.txt")

clustMark<-as.numeric(levels(E7_lim_1k@ident))

nClust<-max(clustMark)

#get an idea of the number of cells
table(E7_lim_1k@ident)

pdf("Plots/ClusterHeatmap.pdf")
DoHeatmap(SubsetData(E7_lim_1k,max.cells.per.ident = 100),
          genes.use = markers_use$gene, slim.col.label = TRUE,remove.key = TRUE)
dev.off()

# Plotting DE genes for clusters ------------------------------------------


if (file.exists("Plots/ClusterDifGenes")){
} else {
  dir.create("Plots/ClusterDifGenes")
}

function(cluster_in,seurat,markers_df){
  cluster_markers <- markers_df %>% filter(cluster == cluster_in) %>% .$gene
  if(length(cluster_markers)==0){return()
    } else{
    pdf(paste0("Plots/ClusterDifGenes/MarkersClust", cluster_in, ".pdf"))
    FeaturePlot(seurat,cluster_markers[1:min(9,length(cluster_markers))],pt.size = 1, ncol = 3)
    def.off()
    }
}


#I think ideally, we would really like to re-write this portion to be a bit more efficient. 
for (x in 0:nClust) {
  m<-markers_all$gene[clustMark == x]
  if (length(m)> 0) {
    if (length(m)>8){
      pdf(paste0("Plots/ClusterDifGenes/MarkersClust", x, ".pdf"))
      FeaturePlot(E7_lim_1k,m[1:9],pt.size = 1)
      dev.off()
    } else {
      add<-rep("ACTB",(9-length(m)))
      pdf(paste0("Plots/ClusterDifGenes/MarkersClust", x, ".pdf"))
      FeaturePlot(E7_lim_1k,c(m,add),pt.size = 1)
      dev.off()
    }
  }
  
}

#okay, so now we need to think of a few things. One is to plot just the degree of expansion of the clones, and I think
#the next thing is looking at shared clones, then lastly looking at the gene differences associated with specific clones

# Make general tsnePlots, first look at expansion at overall level --------------------------------------------------


TRB_expan_count <- Meta_data_1k %>% filter(!is.na(TRB_CDR3)) %>% count(TRB_CDR3) %>% 
  mutate(rank = row_number(desc(n))) %>% arrange(rank) %>% mutate(proportion = n/sum(n))

TRB_expan_count <- TRB_expan_count %>% mutate(Cumulative_proportion = cumsum(proportion))

ggplot(TRB_expan_count %>% filter(rank<=20), aes(rank, Cumulative_proportion)) + geom_step(stat = "identity") +
  geom_point(color = cbPalette[3], size = 0.5)+theme_bw() + labs(x = "TRB Clonal Rank",y = "Cumulative Proportion") +
  xlim(1,20) + ylim(0,NA) + Axis_themes

ggsave("Total_clonal_cumulative.png",plot = last_plot(), device = "png",width = 3, height = 3, units = "in")

TSNE_TRB_expan <- TR_expan_map(Meta_data_1k,ft.size = 0.5) + scale_color_viridis(direction = -1, name = "Clonal size \n(TCRB)", breaks = my_breaks, 
                                                 labels = my_breaks_label, limits = c(0,log2(512))) + labs(x = "tSNE 1", y = "tSNE 2") + Axis_themes


ggsave("TSNE_TRB_expan.png",plot = TSNE_TRB_expan + theme(legend.position = "none"), device = "png",width = 2.8, height = 3.2, units = "in")
ggsave("TSNE_TRB_expan.eps",plot = TSNE_TRB_expan + theme(legend.position = "none"), device = "eps",width = 2.8, height = 3.2, units = "in", family = "ArialMT")

my_breaks <- log2(2^(c(0,seq(9))))
my_breaks_label <- 2^(c(0,seq(9)))

selected_genes_plot <- FeaturePlot(E7_lim_1k, c("GZMB","MYC","NOLC1","NOP16","CCR7","TCF7"), pt.size = 0.5, nCol = 3, do.return = TRUE)
selected_genes_plot <- lapply(selected_genes_plot, `+`, TSNE_theme)

dev.off()

plot_grid(plotlist = selected_genes_plot)
ggsave("selected_genes_2.png",plot = last_plot(), device = "png",width = 4.5, height = 3.5, units = "in")
ggsave("selected_genes_2.eps",plot = last_plot(), device = "eps",width = 4, height = 3.2, units = "in", family = "ArialMT")

#now let's plot the actual TSNE as well
#Here actually it might be easier to just plot from the datapoints, instead of calling Seurat functions
E7_1k_TSNE <- list(TSNEPlot(E7_lim_1k, pt.size = 0.5, do.return = TRUE, group.by = "res.2") + theme_classic() + Axis_themes + labs(x = "tSNE 1", y = "tSNE 2"))

#E7_lim_1k@meta.data$orig.ident.2 <- Meta_data_1k$orig.ident.2

E7_1k_TSNE <- list(E7_1k_TSNE[[1]], TSNEPlot(E7_lim_1k, pt.size = 0.5, group.by = "orig.ident.2",do.return = TRUE) + theme_classic() + Axis_themes + labs(x = "tSNE 1", y = "tSNE 2"))

plot_grid(plotlist = E7_1k_TSNE)

ggsave("E7_1k_TSNE.png",plot = last_plot(), device = "png",width = 8, height = 3.5, units = "in")
ggsave("E7_1k_TSNE.eps",plot = last_plot(), device = "eps",width = 8, height = 3.2, units = "in", family = "ArialMT")

TSNE_TRB_expan_by_animal <- TR_expan_map(Meta_data_1k,TR_var = TRB_count_by_animal,ft.size = 0.5) + scale_color_viridis(direction = -1, name = "Clonal size \n by animal \n(TCRB)", breaks = my_breaks, 
                                                                                 labels = my_breaks_label, limits = c(0,log2(512))) +
  labs(x = "tSNE 1", y = "tSNE 2") + Axis_themes

TSNE_TRB_expan_by_animal <- TSNE_TRB_expan_by_animal + facet_wrap(~orig.ident) + TSNE_theme

ggsave("TSNE_TRB_expan_by_animal.png",plot = TSNE_TRB_expan_by_animal, device = "png",width = 4, height = 4, units = "in")
ggsave("TSNE_TRB_expan_by_animal_legend.eps",plot = TSNE_TRB_expan_by_animal, device = "eps",width = 3, height = 3, units = "in", family = "ArialMT")

# Re-process TSNE for individual animals ----------------------------------

#After showing the process data to Shalek and co, they suggested that trying to do the TSNEing separately, and then
#do 
#save the large seurat file, then remove it from the environment
saveRDS(umi_seurat,"E7_umi_seurat.RDS")
rm(umi_seurat)

#then I think take the m1, m2, m5 mouse out separately,
Meta_data_1k <- as.tibble(Meta_data_1k)

m1_cells <- Meta_data_1k %>% filter(str_detect(cell_id,"^m1")) %>% select(cell_id)

m2_cells <- Meta_data_1k %>% filter(str_detect(cell_id,"^m2")) %>% select(cell_id)

m5_cells <- Meta_data_1k %>% filter(str_detect(cell_id,"^m5")) %>% select(cell_id)

m1_lim_1k <- SubsetData(E7_lim_1k, cells.use = m1_cells$cell_id, subset.raw = TRUE)
m1_lim_1k <- Seurat_process(m1_lim_1k)#no real amount of tinkering really changed the structure of the data.
m1_lim_1k <- TSNE_xy(m1_lim_1k)
m1_lim_1k <- UMAP_xy(m1_lim_1k)
m1_lim_1k <- Add_TRAB_count(seurat = m1_lim_1k)
TR_expan_map(m1_lim_1k@meta.data,ft.size = 1) + Axis_themes
ggsave("m1_TRB_expan.png", device = "png",width = 4, height = 3.2, units = "in")
length(m1_lim_1k@var.genes)
temp <- FeaturePlot(m1_lim_1k, c("GZMB","LRIG1","TNFSF8"), do.return = TRUE, pt.size = 0.5)
temp <- lapply(temp, `+`, TSNE_theme)
plot_grid(plotlist = temp)
ggsave("m1_selected_genes.png",plot = last_plot(), device = "png",width = 3.5, height = 3.5, units = "in")


m2_lim_1k <- SubsetData(E7_lim_1k, cells.use = m2_cells$cell_id, subset.raw = TRUE)
m2_lim_1k <- Seurat_process(m2_lim_1k)
m2_lim_1k <- TSNE_xy(m2_lim_1k)
m2_lim_1k <- Add_TRAB_count(seurat = m2_lim_1k)
length(m2_lim_1k@var.genes)
TR_expan_map(m2_lim_1k@meta.data,ft.size = 1) + Axis_themes
ggsave("m2_TRB_expan.png", device = "png",width = 4, height = 3.2, units = "in")
temp <- FeaturePlot(m2_lim_1k, c("GZMB","LRIG1","TNFSF8"), do.return = TRUE, pt.size = 0.5)
temp <- lapply(temp, `+`, TSNE_theme)
plot_grid(plotlist = temp)
ggsave("m2_selected_genes.png",plot = last_plot(), device = "png",width = 3.5, height = 3.5, units = "in")

m5_lim_1k <- SubsetData(E7_lim_1k, cells.use = m5_cells$cell_id, subset.raw = TRUE)
m5_lim_1k <- Seurat_process(m5_lim_1k)
length(m5_lim_1k@var.genes)
m5_lim_1k <- TSNE_xy(m5_lim_1k)
m5_lim_1k <- Add_TRAB_count(seurat = m5_lim_1k)
TR_expan_map(m5_lim_1k@meta.data,ft.size = 1) + Axis_themes
ggsave("m3_5_TRB_expan.png", device = "png",width = 4, height = 3.2, units = "in")
temp <- FeaturePlot(m5_lim_1k, c("GZMB","LRIG1","TNFSF8"), do.return = TRUE, pt.size = 0.5)
temp <- lapply(temp, `+`, TSNE_theme)
plot_grid(plotlist = temp)
ggsave("m3_5_selected_genes.png",plot = last_plot(), device = "png",width = 3.5, height = 3.5, units = "in")

#I think overall, the structure and the data is definitely there, but I'm not sure if it's a great idea
#to do all separate tSNE, like Alex had suggested.


# Calculate TRB trends in cluster 4 (of the whole TSNE) -------------------
#In the end, it wasn't super easy to just show 4 tsne, since it wasn't exactly just a replica, with different
#number of cells/animals. I think the overall data still makes sense... and it just needs to show some possibilty of
#being useful. Need to work up some figures, and then ask alex/chris what they think.

Total_cluster <- TR_cluster_heatmap(Meta_data_1k,cluster_var = res.2)
Total_cluster <- TR_heatmap_spread(Total_cluster,cluster_var = res.2)
TR_heatmap_plot(Total_cluster)
#The results are actually quite interesting, in that the "intermediate" cluster that we found (cluster 4)
#actually does have special clonotypes, and those clonotypes span across animals, and have similar A/B matching.
#that would actually be the most important thing to show, to emphasize that we are getting shared clones across animals
#with shared transcriptome behaviors

#but first, we need to determine what to do about cluster 10. The best thing to do is to find the markers, then
#see what gene sets tend to be enriched. I think there are a lot of themed proteins, so it's worth looking at.

markers_c10<- FindMarkers(E7_lim_1k, ident.1 = "10",test.use = "bimod")
markers_c10$gene <- rownames(markers_c10)
markers_c10 <- markers_c10 %>% arrange(desc(avg_logFC),p_val_adj) %>% as.tibble()
write_csv(markers_c10,path = "markers_cluster_10.csv")

#took the genelist, and submitted a large subset of it to msigdb, and some of the list doesn't make too much sense,
#but in general, there seems to be some indication that it could be some sort of CD8 memory cell set, and that
#it's some sort of not particularly stimulated cell sets, at least from gene list that are included on there.

#right now, I would list that separately from others, but don't get rid of them.
 
Meta_data_1k <- Meta_data_1k %>% mutate(res.3 = fct_collapse(res.2, Cytotoxic = c("0","1","2","3","5","6","7","8","11"), Myc_activated = c("4"), Naive = "9", Memory = "10")) %>% ungroup()
Meta_data_1k <- Meta_data_1k %>% mutate(res.3 = fct_recode(res.3,"HSP_high" = "Memory"))
Collapsed_cluster <- TR_cluster_heatmap(Meta_data_1k %>% filter(res.3 %in% c("Cytotoxic","Myc_activated")), cluster_var = res.3, clone_cutoff = 10)
Collapsed_cluster <- TR_heatmap_spread(Collapsed_cluster,cluster_var = res.3)

Collapsed_cluster <- Collapsed_cluster[-46,] #get rid of non-functional clone

Collapsed_cluster <- Collapsed_cluster[order(Collapsed_cluster[,1]),]

TR_heatmap_plot(Collapsed_cluster, name = "collapsed_cluster.png",color_levels = magma(20), row_clust = FALSE)

ggsave("TR_cyto_myc_heatmap.png", device = "png",width = 3.5, height = 3.2, units = "in")
ggsave("collapsed_identiy_TSNE.eps", device = "eps",width = 4, height = 3.2, units = "in", family = "ArialMT")

###Chris wants to see a cluster, that looks at where each of the clone is detected.
TRB_dot_grid <- Meta_data_1k %>% filter(TRB_CDR3 %in% clone_list) %>% count(TRB_CDR3,orig.ident.2)

TRB_dot_grid <- TRB_dot_grid %>% mutate(n = ifelse(is.na(n),NA,1))

TRB_dot_grid <- TRB_dot_grid %>% mutate(TRB_CDR3 = factor(TRB_CDR3, levels = rev(clone_list)))

dot_grid <- ggplot(TRB_dot_grid) + geom_point(aes(x = orig.ident.2,y = TRB_CDR3,color = orig.ident.2),size = 1.5) +
  theme_minimal() + Axis_themes + labs(x = "mouse")+ scale_x_discrete(position = "top") +theme(legend.position = "none",
                                                                           axis.title.y = element_blank(),
                                                                           axis.text.y = element_blank())

ggsave("dot_grid.png",device = "png",width = 0.5,height = 3.5, units = "in")
ggsave("dot_grid.eps",device = "eps",width = 0.5,height = 5.25, units = "in", family = "ArialMT")


Collapsed_cluster <- TR_cluster_heatmap(Meta_data_1k %>% filter(res.3 %in% c("Cytotoxic","Myc_activated")), cluster_var = res.3, clone_cutoff = 10)

Myc_TRB <- rownames(Collapsed_cluster)[1:5]

Meta_1k_myc <- Meta_data_1k %>% filter(TRB_CDR3 %in% Myc_TRB) %>% mutate(orig.ident= fct_recode(orig.ident, "m3" = "m5", "m4" = "Tol.Vac.4"))

ggplot(Meta_1k_myc %>% filter(res.3 %in% c("Cytotoxic","Myc_activated"))) + geom_bar(aes(x = orig.ident, fill = orig.ident)) + facet_grid(res.3~TRB_CDR3) + theme_bw()

write_csv(Make_TRAB_table(Meta_1k_myc),"Myc_activated_all_mice.csv")

#after we've made res.3 in our table, we can similarly filter for cytotoxic cells, and look for common clones

Meta_1k_cyt <- Meta_data_1k %>% filter(res.3 == "Cytotoxic") %>% mutate(orig.ident= fct_recode(orig.ident, "m3" = "m5", "m4" = "Tol.Vac.4"))

write_csv(Make_TRAB_table(Meta_1k_cyt),"cytotoxic_all_mice.csv")
#What we see is that there are definitely similar sharing of CDR3 alpha and beta chain in the cytotoxic part as well.
#the question now is how to filter it in a way that allows for graphing.

collapsed_identiy_TSNE<- ggplot(Meta_data_1k) + geom_point(aes(x=Tsne.x, y =Tsne.y, color = res.3), size = 0.5) + Axis_themes +
  labs(x = "tSNE 1", y = "tSNE 2", color = "ClusterID") +scale_color_brewer(type = "qual", palette = 2)
ggsave("collapsed_identiy_TSNE.png", plot = collapsed_identiy_TSNE + theme(legend.position = "none"), device = "png",width = 2.8, height = 3.2, units = "in")
ggsave("collapsed_identiy_TSNE.eps", plot = collapsed_identiy_TSNE + theme(legend.position = "none"), device = "eps",width = 2.8, height = 3.2, units = "in", family = "ArialMT")

Collapsed_cluster_low <- TR_cluster_heatmap(Meta_data_1k, cluster_var = res.3,clone_cutoff = 5)
Collapsed_cluster_low <- TR_heatmap_spread(Collapsed_cluster_low,cluster_var = res.3)
Collapsed_cluster_low <- Collapsed_cluster_low[order(Collapsed_cluster_low[,1]),]
TR_heatmap_plot(Collapsed_cluster_low,name = "collapsed_cluster_low.png",color_levels = magma(20), row_clust = FALSE)
dev.off()

# Meta_data_collap_hmp <- Meta_data_collap_hmp %>% group_by(TRB_CDR3, res.1) %>% summarise(n = sum(n), normalized = sum(normalized))
# 
# Meta_data_collap_hmp_w <- Meta_data_collap_hmp %>% select(TRB_CDR3, res.1, normalized) %>% spread(res.1,normalized)%>% arrange(desc(Cytotoxic))
# collap_hmp_rownames <- Meta_data_collap_hmp_w$TRB_CDR3
# 
# Meta_data_collap_hmp_w  <- Meta_data_collap_hmp_w  %>% ungroup() %>% select(-TRB_CDR3) %>% as.matrix()
# rownames(Meta_data_collap_hmp_w ) <- collap_hmp_rownames

ggplot(Meta_data_1k) +geom_point(aes(Tsne.x,Tsne.y),color = "grey87") + geom_point(data = Meta_1k_myc,aes(Tsne.x,Tsne.y,color = TRB_CDR3))

#integrate re3 back into seurat object 
E7_lim_1k@meta.data$res.3 <- Meta_data_1k$res.3
E7_lim_1k@meta.data$res.2 <- E7_lim_1k@ident

VlnPlot(SubsetData(E7_lim_1k,cells.use = Meta_1k_myc$cell_id), features.plot = "MYC", group.by = "TRB_CDR3")
#I think we see that it is somewhat true that the cells that are in the myc group are higher in myc, no matter which cluster
#they full into on the overall TSNE, but it's also that one of the clonotype is particularly low. Overall, this is 
#telling us that even within these clonotypes, there are differences.
#and of course, just because it's low on MYC, doesn't mean it doesn't have MYC-activated genes.

#let's not dig into it. This is stuff that comes from the review. Ultimately, I think we want to do both myc single-genes
#and modules of myc-related genes.



# TR cluster heatmap by animal --------------------------------------------

#I think maybe generate the cluster heatmap with each mouse, and then we can recombind the matrices, then adjust
#the missing clones, then maybe try to get an idea of what the clonotype distribution might look like, if we break
#it out by animal.
#Several variable we should look at.
View(Meta_data_1k)
View(Collapsed_cluster)
#TR_cluster_heatmap(Meta_data_1k %>% filter(res.3 %in% c("Cytotoxic","Myc_activated")), cluster_var = res.3, clone_cutoff = 10)
#couple of the function we need are example as followed.
#let's get a list of total TRB that we should car about

Collapsed_cluster_list <- rownames(Collapsed_cluster)

#we'd have to first change the Meta_data_1k object to recode the animal identity to the ones that we want to use.

Meta_data_1k <- Meta_data_1k %>% mutate(orig.ident.2= fct_recode(orig.ident, "m3" = "m5", "m4" = "Tol.Vac.4"))

m1_cyt_myc_data <- Meta_data_1k %>% filter(res.3 %in% c("Cytotoxic","Myc_activated")) %>% filter(orig.ident.2 == "m1") %>%
  filter(TRB_CDR3 %in% Collapsed_cluster_list) %>% mutate(res.3 = fct_drop(res.3))

m1_cyt_myc_cluster <- TR_cluster_heatmap(m1_cyt_myc_data, cluster_var = res.3, clone_cutoff = 1)

#m1_cyt_myc_cluster <- TR_heatmap_spread(m1_cyt_myc_cluster,cluster_var = res.3)

View(m1_cyt_myc_cluster)

#TR_heatmap_plot(m1_cyt_myc_cluster)

#can't really easily join, because they are not simple data fucking frame. 

Collapsed_cluster <- TR_cluster_heatmap(Meta_data_1k %>% filter(res.3 %in% c("Cytotoxic","Myc_activated")) %>% mutate(res.3 = fct_drop(res.3)),
                   cluster_var = res.3, clone_cutoff = 10)

collapsed_spread <- Collapsed_cluster %>% select(-n) %>% spread(res.3,normalized)
m1_cyt_myc_spread <- m1_cyt_myc_cluster %>% select(-n) %>% spread(res.3,normalized)

#I think the strategy of spreading, then joining is a good one. let's create a correctly spreaded matrix for each, then join
#I also think it's probably not a great idea to depend on the suffix portion of _join. It's just not extendable.

m2_cyt_myc_data <- Meta_data_1k %>% filter(res.3 %in% c("Cytotoxic","Myc_activated")) %>% filter(orig.ident.2 == "m2") %>%
  filter(TRB_CDR3 %in% Collapsed_cluster_list) %>% mutate(res.3 = fct_drop(res.3))

m2_cyt_myc_cluster <- TR_cluster_heatmap(m2_cyt_myc_data, cluster_var = res.3, clone_cutoff = 1)
m2_cyt_myc_spread <- m2_cyt_myc_cluster %>% select(-n) %>% spread(res.3,normalized)

m3_cyt_myc_data <- Meta_data_1k %>% filter(res.3 %in% c("Cytotoxic","Myc_activated")) %>% filter(orig.ident.2 == "m3") %>%
  filter(TRB_CDR3 %in% Collapsed_cluster_list) %>% mutate(res.3 = fct_drop(res.3))

m3_cyt_myc_cluster <- TR_cluster_heatmap(m3_cyt_myc_data, cluster_var = res.3, clone_cutoff = 1)
m3_cyt_myc_spread <- m3_cyt_myc_cluster %>% select(-n) %>% spread(res.3,normalized)

m4_cyt_myc_data <- Meta_data_1k %>% filter(res.3 %in% c("Cytotoxic","Myc_activated")) %>% filter(orig.ident.2 == "m4") %>%
  filter(TRB_CDR3 %in% Collapsed_cluster_list) %>% mutate(res.3 = fct_drop(res.3))

m4_cyt_myc_cluster <- TR_cluster_heatmap(m4_cyt_myc_data, cluster_var = res.3, clone_cutoff = 1)
m4_cyt_myc_spread <- m4_cyt_myc_cluster %>% select(-n) %>% spread(res.3,normalized)

test_join_1_2 <- full_join(m1_cyt_myc_spread,m2_cyt_myc_spread,by = "TRB_CDR3", suffix = c(".m1",".m2"))

test_join_3_4 <- full_join(m3_cyt_myc_spread,m4_cyt_myc_spread,by = "TRB_CDR3", suffix = c(".m3",".m4"))

m_join <- full_join(test_join_1_2,test_join_3_4,by = "TRB_CDR3")

total_join <- left_join(collapsed_spread,m_join, by = "TRB_CDR3")

rownames(total_join) <- total_join$TRB_CDR3

total_join <- total_join %>% arrange(Cytotoxic.m1, Cytotoxic.m2,Cytotoxic.m3,Cytotoxic.m4)
rownames(total_join) <- total_join$TRB_CDR3
total_join_matrix <- total_join %>% ungroup() %>% select(-1) %>% as.matrix()


TR_heatmap_plot(total_join_matrix,row_clust = FALSE,color_levels = magma(20))


#The plot I think has useful data, but is kind of sparce. I think in order to fit into an interesting
#narrative, we would have to do filtering.
#The Narrative: This is to show that we have similar clones across animals that are driven
#to similar T cell programming. so to do that, first we have to filter 

#one way to find filter(length(unique(orig.ident))>1)

#let's just take total_join, and count the number of NA's across the row

total_join$Animal_no <- total_join %>% is.na() %>% `!` %>% rowSums

total_join <- total_join %>% mutate(Animal_no = (Animal_no-4)/2)

total_join_3_animals <- total_join %>% filter(Animal_no >= 3)

total_join_3_animals <- total_join_3_animals %>% left_join(TRB_expan_count %>% select(TRB_CDR3,n),by = "TRB_CDR3")

total_join_3_animals <- total_join_3_animals %>% mutate(name = paste0(TRB_CDR3," (",n,")"))



rownames(total_join_3_animals) <- total_join_3_animals$name
total_join_3_animals_matrix <- total_join_3_animals %>% ungroup() %>% select(2,4,6,8,10,3,5,7,9,11) %>% as.matrix()

#make sample annotations

sample_col <- data.frame(sample = rep(c("total","m1","m2","m3","m4"),c(2,2,2,2,2)))
rownames(sample_col) <- colnames(total_join_3_animals_matrix)
sample_col$cluster <- rep(c("cytotoxic","Myc activated"),c(5))

TR_heatmap_plot(total_join_3_animals_matrix,row_clust = FALSE,color_levels = magma(20), na_col = "grey40",
                annotation_col= sample_col)
mat_colors <- list(sample = c(gg_color_hue(8)[8],gg_color_hue(4)), cluster = Dark2_4[1:2])
names(mat_colors$sample) <- c(unique(as.character(sample_col$sample)))
names(mat_colors$cluster) <- c(unique(sample_col$cluster))



TRB_3_animal<- pheatmap(total_join_3_animals_matrix,cluster_rows = FALSE,cluster_cols = FALSE,color = magma(20), na_col = "grey40",
                annotation_col= sample_col,border_color =NA, annotation_colors = mat_colors[1:7],
         show_colnames = FALSE, width=8, height=2, fontsize = 6)

#need to make an adjustment, to include the number of time the TRB is detected in the clone name



save_pheatmap_pdf <- function(x, filename, width=8, height=2) {
  pdf(filename, width = width, height = height, fonts = "ArialMT")
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

save_pheatmap_pdf(TRB_3_animal,"3_animal_heatmap_2.pdf", height = 2.2)

#now let's make the accompanying logo plot
total_join_3_animals %>% filter(str_detect(TRB_CDR3,"^CASSQD")) %>% select(TRB_CDR3) -> filtered_logo_input

#let's also include the number of 
filtered_logo_input <- filtered_logo_input %>% left_join(TRB_expan_count %>% select(TRB_CDR3,n),by = "TRB_CDR3")

"CASSQDIDNYAEQFF"->filtered_logo_input$TRB_CDR3[4]
#remember that leucine (I) here is just a place hoder for gap.


logo_TRB_freq <- rep(filtered_logo_input$TRB_CDR3,filtered_logo_input$n)
logo_3_animal <- ggseqlogo(filtered_logo_input$TRB_CDR3, method = "probability")
filtered_logo_input <- filtered_logo_input %>% mutate(proportion = n/sum(n))
logo_3_animal <- ggseqlogo(logo_TRB_freq, method = "probability")

logo_3_animal + Axis_themes

ggsave("logo_3_animal.png",plot = logo_3_animal + Axis_themes, device = "png",width = 3, height = 1.8, units = "in", family = "ArialMT")
ggsave("logo_3_animal.pdf",plot = logo_3_animal + Axis_themes, device = "pdf",width = 3, height = 1.8, units = "in", family = "ArialMT")


#Next thing to do would be to see if the alpha chain would have similar congruency. To do that we need to maybe first
#build the relational matrix, which should give us the right matrix to use for these calculation

relation_3_animal <- Meta_data_1k %>% filter(TRB_CDR3 %in% total_join_3_animals$TRB_CDR3) %>% 
  filter(!str_detect(TRA_CDR3,"\\*"), !str_detect(TRB_CDR3,"\\*")) %>% filter(!is.na(TRB_CDR3),!is.na(TRA_CDR3))

TRB_r3_names <- relation_3_animal %>% count(TRB_CDR3) %>% mutate(name = paste0(TRB_CDR3," (",n,")")) %>% select(-n)
TRA_r3_names <- relation_3_animal %>% count(TRA_CDR3) %>% mutate(name = paste0(TRA_CDR3," (",n,")")) %>% select(-n)

relation_3_animal_pair <- relation_3_animal %>% count(TRB_CDR3,TRA_CDR3,orig.ident.2)
relation_3_animal_pair <- relation_3_animal_pair %>% inner_join(TRB_r3_names,by = "TRB_CDR3") %>% inner_join(TRA_r3_names, by = "TRA_CDR3") %>% select(name.x,name.y,orig.ident.2,n)
relation_3_animal_pair <- relation_3_animal_pair %>% rename(TRB_CDR3 = name.x, TRA_CDR3 = name.y)

relation_3_animal_pair_filtered <- relation_3_animal_pair %>% filter(str_detect(TRB_CDR3, "^CASSQD"))
relation_3_animal_pair_filtered <- relation_3_animal_pair_filtered %>% mutate(logn = log2(n+1))
relation_3_animal_pair_filtered <- relation_3_animal_pair_filtered %>% filter(n>1)


library(circlize)

dev.off()
circos.clear()
par() -> default_par
par(mar = default_par$mar - 2)
small_gap <-  2
big_gap <-  20
gg_colors <- gg_color_hue(4)

#grid.col <- c(brewer.pal(8,"Set1"), brewer.pal(7,"Set2"), brewer.pal(3,"Accent")[-2])

col_link <- relation_3_animal_pair_filtered %>% mutate(colors = recode(orig.ident.2, "m1" = gg_colors[1],
                                                    "m2" = gg_colors[2],
                                                    "m3" = gg_colors[3],
                                                    "m4" = gg_colors[4])) %>% .$colors
border_link <- relation_3_animal_pair_filtered %>% mutate(border = NA) %>% mutate(border = ifelse(TRB_CDR3 %in% c("CASSQDLGNYAEQFF (364)"),"black",NA)) %>% .$border


circos.par(gap.after = c(rep(small_gap, length(unique(relation_3_animal_pair_filtered$TRB_CDR3))-1), big_gap, 
                         rep(small_gap, length(unique(relation_3_animal_pair_filtered$TRA_CDR3))-1), big_gap),
           start.degree = -10)


chordDiagram(relation_3_animal_pair_filtered[c(-3,-4)], annotationTrack = "grid", 
             preAllocateTracks = 1,
             link.sort = TRUE,
             link.decreasing = FALSE,
             grid.col = NULL,#grid.col,
             col = col_link,
             link.lwd = 2,
             link.border = border_link)
circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  sector.name = get.cell.meta.data("sector.index")
  
  # print labels & text size (cex)
  circos.text(mean(xlim), ylim[1] + .2, sector.name,
              facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5), cex=0.6, family = "sans")
  
  # print axis
  #circos.axis(h = "top", labels.cex = 0.5, major.tick.percentage = 0.2,
  #sector.index = sector.name, track.index = 2)
}, bg.border = NA)

circos.clear()
abline(h = -0.05, lty = 5, col = "#00000080")

text(1.34,0,"TCR\u03B1")
text(1.34,-0.09,"TCR\u03B2")

#1225x900 in size
dev.off()



#This is saved in a rather, interesting 
# Volcano, cumulative expansion -------------------------------------------
#First thing to do is to generate some DE results, then plot them.
#First let's do the myc vs cytotoxic
E7_lim_1k <- SetAllIdent(E7_lim_1k, id = "res.3")

Cyto_v_Myc <- DE_test(E7_lim_1k, group_1 = "Cytotoxic",group_2 = "Myc_activated")

Cyto_v_Myc <- Cyto_v_Myc %>% mutate(avg_logFC = avg_logFC * log10(exp(1))) # change from natural to base-10 log

Cyto_v_Myc_volcano <- plot_volcano(Cyto_v_Myc,gene_list_1 = c("CCL5","CCR2","GZMA","GZMB","CXCR6"),
                                   gene_list_2 = c("TNFSF8","TNFRSF4","LRIG1","MYC"))

Cyto_v_Myc_volcano <- Cyto_v_Myc_volcano + Axis_themes + labs(x = "log10(fold change)")

ggsave("Cyto_v_Myc_volcano.eps",plot = Cyto_v_Myc_volcano, device = "eps",width = 4, height = 3.2, units = "in", family = "ArialMT")
ggsave("Cyto_v_Myc_volcano.png",plot = Cyto_v_Myc_volcano, device = "png",width = 4, height = 3.2, units = "in")

Mem_v_all <- DE_test(E7_lim_1k, group_1 = "Memory",group_2 = c("Myc_activated","Cytotoxic"))

Heat_shock_genes <- Mem_v_all %>% filter(significance == "significant") %>% filter(str_detect(genes,"^HSP")) %>% select(genes)

Mem_v_all <- Mem_v_all %>% mutate(avg_logFC = avg_logFC * log10(exp(1))) # change from natural to base-10 log

Mem_v_all_volcano <- plot_volcano(Mem_v_all, gene_list_1 = Heat_shock_genes$genes)

Mem_v_all_volcano <- Mem_v_all_volcano + Axis_themes + labs(x = "log10(fold change)")

ggsave("Mem_v_all_volcano.eps",plot = Mem_v_all_volcano, device = "eps",width = 4, height = 3.2, units = "in", family = "ArialMT")
ggsave("Mem_v_all_volcano.png",plot = Mem_v_all_volcano, device = "png",width = 4, height = 3.2, units = "in")

write_csv(Cyto_v_Myc,"Cyto_v_Myc.csv")
write_csv(Mem_v_all,"high_HSP_v_all.csv")

write_csv(Cyto_v_Myc %>% filter(significance == "significant",avg_logFC>0) %>% select(genes),"Cyto_genes.csv")

write_csv(Cyto_v_Myc %>% filter(avg_logFC>0) %>% arrange(p_val_adj) %>% select(genes,p_val_adj),"Cyto_genes_2.csv")

write_csv(Cyto_v_Myc %>% filter(significance == "significant",avg_logFC<0) %>% select(genes),"Myc_genes.csv")

write_csv(Cyto_v_Myc %>% filter(avg_logFC<0) %>% arrange(p_val_adj) %>% select(genes,p_val_adj),"Myc_genes_2.csv")


#filter and write genes into specific gene lists for easy retrieval

#now make cumulative proportion plots

TRB_expan_by_animal <- TRB_expan_by_animal %>% group_by(orig.ident) %>% arrange(desc(TRB_count_by_animal), .by_group = TRUE) %>%
  mutate(rank = row_number(desc(TRB_count_by_animal)))

TRB_expan_by_animal <- TRB_expan_by_animal %>% mutate(proportion = TRB_count_by_animal/sum(TRB_count_by_animal))

TRB_expan_by_animal <- TRB_expan_by_animal %>% mutate(Cumulative_proportion = cumsum(proportion))

TRB_expan_by_animal <- TRB_expan_by_animal %>% ungroup() %>% mutate(orig.ident= fct_recode(orig.ident, "m3" = "m5", "m4" = "Tol.Vac.4"))

#now plot the results
ggplot(TRB_expan_by_animal %>% filter(rank<=20), aes(rank, Cumulative_proportion, color = orig.ident)) + geom_step(stat = "identity") +
  geom_point(color = "black", size = 0.5)+theme_bw() + labs(x = "TRB Clonal Rank",y = "Cumulative Proportion") +
  xlim(1,20) + ylim(0,NA) + Axis_themes #+ facet_wrap(~orig.ident)

ggsave("Total_clonal_cumulative.png",plot = last_plot(), device = "png",width = 3.5, height = 2.8, units = "in")
ggsave("Total_clonal_cumulative.eps",plot = last_plot(), device = "eps",width = 3.5, height = 2.8, units = "in", family = "ArialMT")

#Make calculation from the cumulative_prop of each of the animals

std_mean <- TRB_expan_by_animal %>% group_by(rank) %>% summarise(average_cmp = mean(Cumulative_proportion), std_dev = sd(Cumulative_proportion),
                                                                 Min = min(Cumulative_proportion, na.rm = TRUE),
                                                                 Max = max(Cumulative_proportion, na.rm = TRUE)) %>% ungroup()

write_csv(x = std_mean,path = "cumulative_rank_stat.csv")

#another thing we should calculate is how many cells are scored in the top 20, just to give the analysis more credence

Scoring_cells <- TRB_expan_by_animal %>% group_by(orig.ident) %>% 
  filter(rank<=20) %>% summarise(cells = sum(TRB_count_by_animal))
Scoring_cells_stat <- Scoring_cells %>% ungroup() %>% summarise(mean = mean(cells), std_dev = sd(cells), min = min(cells), max = max(cells))

# Circus Plots ------------------------------------------------------------

library(RColorBrewer)
library(circlize)
par() -> default_par
par(mar = default_par$mar - 2)
small_gap <-  2
big_gap <-  20
gg_colors <- gg_color_hue(4)
circos.clear()

TRB_names <- Meta_1k_myc %>% filter(!is.na(TRB_CDR3),!is.na(TRA_CDR3)) %>% count(TRB_CDR3) %>% mutate(name = paste0(TRB_CDR3," (",n,")")) %>% select(-n)
TRA_names <- Meta_1k_myc %>% filter(!is.na(TRB_CDR3),!is.na(TRA_CDR3)) %>% count(TRA_CDR3) %>% mutate(name = paste0(TRA_CDR3," (",n,")")) %>% select(-n)


relation <- Meta_1k_myc %>% count(TRB_CDR3,TRA_CDR3,orig.ident) %>% filter(!is.na(TRB_CDR3)) %>% filter(!is.na(TRA_CDR3)) #%>% arrange(desc(n))
relation <- relation %>% filter(!str_detect(TRA_CDR3,"\\*"))
relation <- relation %>% left_join(TRB_names,by = "TRB_CDR3") %>% left_join(TRA_names, by = "TRA_CDR3") %>% select(name.x,name.y,orig.ident,n)
relation <- relation %>% rename(TRB_CDR3 = name.x, TRA_CDR3 = name.y)


grid.col <- c(brewer.pal(8,"Set1"), brewer.pal(7,"Set2"), brewer.pal(3,"Accent")[-2])
col_link <- relation %>% mutate(colors = recode(orig.ident, "m1" = gg_colors[1],
                                                "m2" = gg_colors[2],
                                                "m3" = gg_colors[3],
                                                "m4" = gg_colors[4])) %>% .$colors
border_link <- relation %>% mutate(border = NA) %>% mutate(border = ifelse(TRB_CDR3 %in% c("CASSDSGGNTEVFF (15)","CASSQSSYEQYF (7)"),"black",NA)) %>% .$border


circos.par(gap.after = c(rep(small_gap, length(unique(relation$TRB_CDR3))-1), big_gap, 
                         rep(small_gap, length(unique(relation$TRA_CDR3))-1), big_gap),
           start.degree = -16)


chordDiagram(relation[-3], annotationTrack = "grid", 
             preAllocateTracks = 1,
             link.sort = TRUE,
             link.decreasing = FALSE,
             grid.col = grid.col,
             col = col_link,
             link.lwd = 2,
             link.border = border_link)
circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  sector.name = get.cell.meta.data("sector.index")

  # print labels & text size (cex)
  circos.text(mean(xlim), ylim[1] + .2, sector.name,
              facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5), cex=1, family = "sans")

  # print axis
  #circos.axis(h = "top", labels.cex = 0.5, major.tick.percentage = 0.2,
              #sector.index = sector.name, track.index = 2)
}, bg.border = NA)

circos.clear()
abline(h = -0.05, lty = 5, col = "#00000080")

text(1.34,0,"TCR\u03B1")
text(1.34,-0.09,"TCR\u03B2")

#1225x900 in size
dev.off()


# Circos for cytotoxic cells ----------------------------------------------
require(RColorBrewer)
require(circlize)

cyt_circos<- Meta_1k_cyt %>% filter(TRB_count >= 10) %>% filter(!str_detect(TRA_CDR3,"\\*"), !str_detect(TRB_CDR3,"\\*"))
cyt_circos <- cyt_circos %>% filter(!is.na(TRB_CDR3),!is.na(TRA_CDR3))


TRB_cyt_names <- cyt_circos %>% count(TRB_CDR3) %>% mutate(name = paste0(TRB_CDR3," (",n,")")) %>% filter(n>1) %>% select(-n)
TRA_cyt_names <- cyt_circos %>% count(TRA_CDR3) %>% mutate(name = paste0(TRA_CDR3," (",n,")")) %>% filter(n>1) %>% select(-n)

relation_cyt <- cyt_circos %>% count(TRB_CDR3,TRA_CDR3,orig.ident)
relation_cyt <- relation_cyt %>% inner_join(TRB_cyt_names,by = "TRB_CDR3") %>% inner_join(TRA_cyt_names, by = "TRA_CDR3") %>% select(name.x,name.y,orig.ident,n)
relation_cyt <- relation_cyt %>% rename(TRB_CDR3 = name.x, TRA_CDR3 = name.y)

relation_cyt <- relation_cyt %>% filter(n>1)#further filter the links we want to show
relation_cyt <- relation_cyt %>% mutate(n = log2(n+1))

#even after all this, the relationship is still too complicated. Will just filter down to the number
#of chains that shows shared sharing across animals.

relation_cyt <- relation_cyt %>% group_by(TRA_CDR3) %>% filter(length(unique(orig.ident))>1) %>% ungroup()
relation_cyt <- relation_cyt %>% group_by(TRB_CDR3) %>% filter(length(unique(orig.ident))>1) %>% ungroup()

dev.off()
circos.clear()
par() -> default_par
par(mar = default_par$mar - 2)
small_gap <-  2
big_gap <-  20
gg_colors <- gg_color_hue(4)

#grid.col <- c(brewer.pal(8,"Set1"), brewer.pal(7,"Set2"), brewer.pal(3,"Accent")[-2])

col_link <- relation_cyt %>% mutate(colors = recode(orig.ident, "m1" = gg_colors[1],
                                                "m2" = gg_colors[2],
                                                "m3" = gg_colors[3],
                                                "m4" = gg_colors[4])) %>% .$colors
#border_link <- relation %>% mutate(border = NA) %>% mutate(border = ifelse(TRB_CDR3 %in% c("CASSDSGGNTEVFF (15)","CASSQSSYEQYF (7)"),"black",NA)) %>% .$border


circos.par(gap.after = c(rep(small_gap, length(unique(relation_cyt$TRB_CDR3))-1), big_gap, 
                         rep(small_gap, length(unique(relation_cyt$TRA_CDR3))-1), big_gap),
           start.degree = -10)


chordDiagram(relation_cyt[-3], annotationTrack = "grid", 
             preAllocateTracks = 1,
             link.sort = TRUE,
             link.decreasing = FALSE,
             grid.col = NULL,#grid.col,
             col = col_link,
             link.lwd = 2)#,
             #link.border = border_link)
circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  sector.name = get.cell.meta.data("sector.index")
  
  # print labels & text size (cex)
  circos.text(mean(xlim), ylim[1] + .2, sector.name,
              facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5), cex=0.6, family = "sans")
  
  # print axis
  #circos.axis(h = "top", labels.cex = 0.5, major.tick.percentage = 0.2,
  #sector.index = sector.name, track.index = 2)
}, bg.border = NA)

circos.clear()
abline(h = -0.05, lty = 5, col = "#00000080")

text(1.34,0,"TCR\u03B1")
text(1.34,-0.09,"TCR\u03B2")

#1225x900 in size
dev.off()

# MsigDb parsing ----------------------------------------------------------
##Trying to read in the MsigDb results
MsigDb <- read_csv("Overlap_genesets.csv")
MsigDb <- MsigDb[1:40,]
MsigDb <- MsigDb %>% mutate(Cluster_subset = factor(Cluster_subset),Genes_Gene_Set_K = as.integer(Genes_Gene_Set_K))
#log the q value to -log_10
MsigDb <- MsigDb %>% mutate(log_qvalue = -log10(FDR_qvalue))
MsigDb <- MsigDb %>% mutate(Description_short = ifelse(is.na(Description_short),"",Description_short))

Dark2_4 <- brewer.pal(4,"Dark2")
ggplot(MsigDb %>% filter(Cluster_subset == "Cytotoxic"),aes(x = fct_reorder(Systematic_name,log_qvalue), y = log_qvalue)) + 
  geom_bar(fill = Dark2_4[1], stat = "identity") + coord_flip() + geom_text(aes(label = Description_short), y = 2.1, hjust = 0, size = 2) + Axis_themes +
  labs(y = "-log10(FDR q-value)", x = "Gene set systematic ID") + geom_hline(yintercept = -log10(0.01), linetype = 2) +
  scale_y_continuous(breaks =c(0,2,10,20,30))


ggsave("Cyto_Geneset.eps",plot = last_plot(), device = "eps",width = 5, height = 3, units = "in", family = "ArialMT")
#ehh not sure if the gene-set names actually make any sense. Need better way to shorten the description
ggsave("Cyto_Geneset.png",plot = last_plot(), device = "png",width = 5, height = 3, units = "in")

ggplot(MsigDb %>% filter(Cluster_subset == "Myc_activated"),aes(x = fct_reorder(Systematic_name,log_qvalue), y = log_qvalue)) + 
  geom_bar(fill = Dark2_4[2], stat = "identity") + coord_flip() + geom_text(aes(label = Description_short), y = 2.1, hjust = 0, size = 2) + Axis_themes +
  labs(y = "-log10(FDR q-value)", x = "Gene set systematic ID") + geom_hline(yintercept = -log10(0.01), linetype = 2) +
  scale_y_continuous(breaks =c(0,2,10,20,30),limits = c(0,30))

ggsave("Myc_Geneset.png",plot = last_plot(), device = "png",width = 5, height = 3, units = "in")
ggsave("Myc_Geneset.eps",plot = last_plot(), device = "eps",width = 5, height = 3, units = "in", family = "ArialMT")

# more analysis... --------------------------------------------------------
# one thing that we should start to write is calculating sharing of alpha/beta chain across samples/animals.
#that's more or less done for now, at least we've picked out the interesting features, and the interesting versions
#Now we have two alpha and beta that are apparently super commonly shared. We should search for those
#and just tabulate how frequent those are.

Cyt_Special_clone <- Meta_data_1k %>% filter(TRB_CDR3 == "CASSQDLGNYAEQFF") %>% filter(!is.na(TRA_CDR3)) %>% 
  count(TRB_CDR3,TRA_CDR3,TRA.2_CDR3)

#double TRA was observed 110 times, and single A of either chain was found 39 and 188 times.

#we also want to know
Cyt_Special_clone_3 <- Meta_data_1k %>% filter(TRB_CDR3 == "CASSQDLGNYAEQFF") %>% filter(!is.na(TRA_CDR3)) %>% 
  count(TRB_CDR3,TRA_CDR3,TRA.2_CDR3,orig.ident.2, TRAV,TRAJ,TRBV,TRBJ,TRAV.2,TRAJ.2)

write_csv(Cyt_Special_clone_3,path = "Special_clone_2.csv")
#assume that we start with 

# Logo Plots --------------------------------------------------------------
require(ggseqlogo)

#let's start with the cytotoxic cells
View(Meta_1k_cyt)

Cyt_TRB_logo <- Meta_1k_cyt %>% filter(!is.na(TRB_CDR3),!str_detect(TRB_CDR3,"\\*")) %>% mutate(TRB_length = str_length(TRB_CDR3)) %>%
  select(TRB_CDR3,TRB_length) %>% arrange(desc(TRB_length))

Cyt_TRB_logo_collapsed <- Cyt_TRB_logo %>% distinct(TRB_CDR3,TRB_length)

Cyt_TRB_logo_collapsed <- Cyt_TRB_logo_collapsed %>% mutate(TRB_CDR3_trimmed = str_sub(TRB_CDR3,4,-4))

#Need to write function to do this more quickly.

Sequence_logo <- function(Data, TR_column, TR_length, length, method_type = "bits"){
  TR_column <- enquo(TR_column)
  TR_length <- enquo(TR_length)
  Data %>% filter((!!TR_length) == length) %>% select(!!TR_column) %>% pull(!!TR_column)
  logo <- ggseqlogo(Data %>% filter(!!TR_length == length) %>% select(!!TR_column) %>% pull(!!TR_column), method = method_type)
  return(logo)
}


Sequence_logo(Cyt_TRB_logo_collapsed,TR_column = TRB_CDR3_trimmed,TR_length = TRB_length, length = 17)
#I think it might be a good observation to have, just in case someone ask, but nmer analysis would be more generalizable
Sequence_logo(Cyt_TRB_logo_collapsed,TR_column = TRB_CDR3_trimmed,TR_length = TRB_length, length = 16)
Sequence_logo(Cyt_TRB_logo_collapsed,TR_column = TRB_CDR3_trimmed,TR_length = TRB_length, length = 15)
Sequence_logo(Cyt_TRB_logo_collapsed,TR_column = TRB_CDR3_trimmed,TR_length = TRB_length, length = 14)

#it seems liek overall legnth of 17 has the highest agreement for conserved regions...

Cyt_TRB_logo_distr <- Meta_1k_cyt %>% filter(!is.na(TRB_CDR3),!str_detect(TRB_CDR3,"\\*")) %>% mutate(TRB_length = str_length(TRB_CDR3))

ggplot(Cyt_TRB_logo_distr, aes(x = TRB_length)) + geom_bar() #shows that the mean of the length by cells is actually legnth
#15 and 14.

ggplot(Cyt_TRB_logo_collapsed, aes(x = TRB_length)) + geom_bar() #not even by collapsed number of unique CDR3 sequences




# Recovery percentage -----------------------------------------------------
# One last figure that we at least want in our back-pocket is the recovery rate for each of the 
# the annoying thing is that we just have to re-run the old code for the general TCR seq-well code.

#recreate umi_seuart object for the sake of easy running of the previously written code

umi_seurat <- E7_cleaned_seurat

#then run of the rest of the TCR_SeqWell code, but then we still need an overall percentage.


Coverage_summary_plot <- ggplot(Coverage_summary_sample_long,aes(variable, value*100, color = variable, fill = variable)) + 
  geom_boxplot(fill = NA) + geom_jitter(pch = 21, color = "black", width = 0.05, height = 0) + scale_fill_manual(values = cbPalette[c(8,3,2)]) +
  scale_color_manual(values = cbPalette[c(8,3,2)]) +
  labs(x = "",y = "Total TCR recovery (% of total cells)", color = "TCR mapping", fill = "TCR mapping") + 
  theme_classic() + Axis_themes + ylim(5,80)

ggsave("Coverage_summary.png",plot = last_plot(), device = "png", width = 3.2, height = 2.8, units = "in")
ggsave("Coverage_summary.eps",plot = last_plot(), device = "eps", width = 3.2, height = 2.8, units = "in", family = "ArialMT")

#to be save, we'll save a separate meta data file tibble, so that we don't interfere with anything that's already made

coverage_meta_data <- umi_seurat@meta.data %>% as.tibble()

TRA_tempt <- coverage_meta_data %>% filter(as.integer(as.character(TRAC_bin)) > 0) %>% group_by(orig.ident) %>% count(TRA_Recovery) %>% rename(TRA = n)
TRB_tempt <- coverage_meta_data %>% filter(as.integer(as.character(TRBC2_bin)) > 0) %>% group_by(orig.ident) %>% count(TRB_Recovery) %>% rename(TRB = n)

TRA_tempt <- TRA_tempt %>% mutate(TRA_rate = TRA/sum(TRA))
TRB_tempt <- TRB_tempt %>% mutate(TRB_rate = TRB/sum(TRB))

TRA_tempt <- TRA_tempt %>% mutate(TRA_Recovery = factor(TRA_Recovery))
TRB_tempt <- TRB_tempt %>% mutate(TRB_Recovery = factor(TRB_Recovery))

TRA_tempt <- TRA_tempt %>% rename(TCR_Recovery = TRA_Recovery, Rate = TRA_rate)
TRB_tempt <- TRB_tempt %>% rename(TCR_Recovery = TRB_Recovery, Rate = TRB_rate)

TRAB_tempt <- coverage_meta_data %>% filter(as.integer(as.character(TRBC2_bin)) > 0,as.integer(as.character(TRAC_bin)) > 0) %>% group_by(orig.ident) %>% 
  count(TRA_Recovery,TRB_Recovery) %>% rename(TRAB = n)

TRAB_tempt <- TRAB_tempt %>% mutate(TRAB_recovery = ifelse(TRA_Recovery == "TRA Recovered" & TRB_Recovery == "TRB Recovered", "TRAB Recovered", "No Recovery")) %>%
  select(-TRA_Recovery, -TRB_Recovery)

TRAB_tempt <- TRAB_tempt %>% mutate(TRAB_recovery = factor(TRAB_recovery))
TRAB_tempt <- TRAB_tempt %>% group_by(orig.ident, TRAB_recovery) %>% summarise(TRAB = sum(TRAB)) %>% ungroup() %>%
  group_by(orig.ident) %>% mutate(TRAB_rate = TRAB/sum(TRAB)) %>% ungroup()

TRAB_tempt <- TRAB_tempt %>% rename(TCR_Recovery = TRAB_recovery, Rate = TRAB_rate)

Coverage_summary_normalized <- bind_rows(TRA_tempt %>% select(-TRA) %>% filter(!(TCR_Recovery == "No Recovery")), 
                                         TRB_tempt %>% select(-TRB) %>% filter(!(TCR_Recovery == "No Recovery")),
                                         TRAB_tempt %>% select(-TRAB) %>% filter(!(TCR_Recovery == "No Recovery")))

Coverage_summary_normalized <- Coverage_summary_normalized %>% ungroup() %>% mutate(TCR_Recovery = factor(TCR_Recovery))

Coverage_summary_normalized <- Coverage_summary_normalized %>% mutate(TCR_Recovery = fct_recode(TCR_Recovery, "TRA" = "TRA Recovered",
                                                                                                "TRB" = "TRB Recovered",
                                                                                                "TRAB" = "TRAB Recovered"))

Coverage_summary_normalized <- Coverage_summary_normalized %>% mutate(TCR_Recovery = fct_relevel(TCR_Recovery,"TRA","TRB"))

Coverage_summary_Normalized_plot <- ggplot(Coverage_summary_normalized,aes(TCR_Recovery, Rate*100, color = TCR_Recovery, fill = TCR_Recovery)) + 
  geom_boxplot(fill = NA) + geom_jitter(pch = 21, color = "black", width = 0.05, height = 0) + scale_fill_manual(values = cbPalette[c(8,3,2)]) +
  scale_color_manual(values = cbPalette[c(8,3,2)]) +
  labs(x = "",y = "Total TCR recovery (% of filtered cells)", color = "TCR mapping", fill = "TCR mapping") + 
  theme_classic() + Axis_themes + ylim(10,80)


ggsave("Coverage_summary_filtered.png",plot = last_plot(), device = "png", width = 3.2, height = 2.8, units = "in")
ggsave("Coverage_summary_filtered.eps",plot = last_plot(), device = "eps", width = 3.2, height = 2.8, units = "in", family = "ArialMT")


# plot specific clones ----------------------------------------------------

A_CDR <- c("CASSTDRANTEVFF","CASSQSSYEQYF","CASRDDGSYNSPLYF","CASSQDDNYAEQFF")

CDR3_A_graph <- arrange(Meta_data_1k,TRB_count) %>% filter(TRB_CDR3 %in% A_CDR) %>% 
  mutate(TRB_CDR3 = factor(TRB_CDR3, levels = c("CASSTDRANTEVFF","CASSQSSYEQYF","CASRDDGSYNSPLYF","CASSQDDNYAEQFF")))

A_color <- c("#fdbe85","#fd8d3c","#e6550d","#a63603")

CDR3_B_graph <- arrange(Meta_data,TRB_expan_count) %>% filter(TRB_CDR3 %in% B_CDR) %>% 
  mutate(TRB_CDR3 = factor(TRB_CDR3, levels = c("CASSLSSPSGNTLYF", "CASSQDGGNYAEQFF","CASSLSEYEQYF","CASSQDWGDTYAEQFF","CASSYTEQDTQYF","CASSPGQGEQAPLF")))

B_color <- c("#c6dbef","#9ecae1","#6baed6","#4292c6","#2171b5","#084594")

TSNE_TRB_CDR3_A <- ggplot(arrange(Meta_data_1k,TRB_count), aes(x = Tsne.x, y = Tsne.y)) + 
  geom_point(color = "grey87", size = 0.5) + 
  geom_point(data = CDR3_A_graph, 
             aes(color = TRB_CDR3), size = 1) + Axis_themes + theme(axis.text.y = element_blank(), 
                                                                    axis.text.x = element_blank(), 
                                                                    axis.ticks.x= element_blank(),
                                                                    axis.ticks.y= element_blank()) + labs(x = "tSNE 1", y = "tSNE 2", color = "TCR clonotypes")
  #scale_color_manual(values = A_color) +  #labs(color = "Group A TRB CDR3") 

ggsave("TSNE_TRB_specific.png",plot = TSNE_TRB_CDR3_A, device = "png",width = 4, height = 3.2, units = "in")
TSNE_TRB_CDR3_B <- ggplot(arrange(Meta_data,TRB_expan_count), aes(x = Tsne.x, y = Tsne.y)) + 
  geom_point(color = "grey87", size = 0.5) + 
  geom_point(data = CDR3_B_graph, 
             aes(color = TRB_CDR3), size = 1) +
  scale_color_manual(values = B_color) + TSNE_theme #labs(color = "Group B TRB CDR3")

# RNAseqQC ----------------------------------------------------------------
RNASeq_QC <- VlnPlot(E7_lim_1k, c("nGene", "nUMI", "reads", "percent.mito"), group.by = "orig.ident", 
                     x.lab.rot = TRUE, nCol = 2, size.title.use = 10, point.size.use = 0, do.return = TRUE,return.plotlist = TRUE)

QC_titles <- list("Number of genes","Number of UMI","Number of Reads","Percent of Mitochondrial Genes")

RNASeq_QC <- map2(RNASeq_QC, QC_titles, function(test,QC_titles) test + labs(title = QC_titles, x = "") + 
                    Axis_themes + theme(plot.margin = unit(c(0, 0, 0, 0), "in")))

RNASeq_QC <- map(RNASeq_QC, function(test) test + geom_boxplot(fill = 'white',width = 0.2, outlier.shape = NA))

save_plot("RNASeq_QC.png", plot_grid(plotlist = RNASeq_QC), base_height = 4, unit = "in")
save_plot("RNASeq_QC.pdf", plot_grid(plotlist = RNASeq_QC), base_height = 4, unit = "in", family = "ArialMT")


# TCR_QC ------------------------------------------------------------------

# Clean the data for combination ------------------------------------------

TRA_sliced <- TRA_df %>% group_by(BC) %>% slice(1) %>% rename(TRA_CDR3 = CDR3, 
                                                                   UMI_count_TRA = UMI_count,
                                                                   TRA_reads = nReads) %>% ungroup()

TRB_sliced <- TRB_df %>% group_by(BC) %>% slice(1) %>% rename(TRB_CDR3 = CDR3,
                                                                   UMI_count_TRB = UMI_count,
                                                                   TRB_reads = nReads) %>% ungroup()


# test out putting together the Meta_data_1k and TRAB data ------------------

Meta_data_1k <- Meta_data_1k %>% left_join(TRA_sliced %>% select(BC,TRA_CDR3,UMI_count_TRA,TRA_reads), by = c("BC","TRA_CDR3"))

Meta_data_1k <- Meta_data_1k %>% left_join(TRB_sliced %>% select(BC,TRB_CDR3,UMI_count_TRB,TRB_reads), by = c("BC","TRB_CDR3"))


# making graphs of the Reads and UMI count and reads ----------------------
Color_map_range<- brewer.pal(9,"Blues")
Color_map <- colorRampPalette(c(Color_map_range[c(9,1)]))

my_breaks <- log10(10^(c(0,seq(5))))
my_breaks_label <- scientific(10^(c(0,seq(5))))

QC_list <- vector("list",4)

names(QC_list) <- c("TRA_reads","TRB_reads","TRA_UMI","TRB_UMI")

QC_list[[1]] <- ggplot(Meta_data_1k, aes(x = orig.ident.2, y = log10(TRA_reads), fill = orig.ident)) + geom_violin(scale = "width") +
  scale_y_continuous(breaks = my_breaks, labels = my_breaks_label) + labs(x = '',y = 'Number of reads', title = 'TCRA reads per UMI') + theme_classic() + Axis_themes + 
  theme(legend.position = 'none',axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5))

QC_list[[2]] <- ggplot(Meta_data_1k, aes(x = orig.ident.2, y = log10(TRB_reads), fill = orig.ident)) + geom_violin(scale = "width") +
  scale_y_continuous(breaks = my_breaks, labels = my_breaks_label) + labs(x = '',y = 'Number of reads', title = 'TCRB reads per UMI') + theme_classic() + Axis_themes + 
  theme(legend.position = 'none',axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5))
#hmm I don't think this looks particularly good. Since it's hard to see the low end of the scale. I think it just has to be chnaged. and
#that means I might have to change the others as well. or I have to change the UMI 

QC_list[[3]] <- ggplot(Meta_data_1k %>% filter(!is.na(UMI_count_TRA)), aes(x = orig.ident.2, fill = factor(UMI_count_TRA, ordered = TRUE))) + 
  geom_bar(position = position_fill(reverse = TRUE)) + scale_fill_viridis_d("UMI count", direction = 1, option = "plasma") + 
  labs(x = '',y = 'Fraction of cells (0-1)', title = 'Number of TCRA UMI per cell') + theme_classic() + Axis_themes + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1), plot.title = element_text(hjust = 0.5))


QC_list[[4]] <- ggplot(Meta_data_1k %>% filter(!is.na(UMI_count_TRB)), aes(x = orig.ident.2, fill = factor(UMI_count_TRB, ordered = TRUE))) + 
  geom_bar(position = position_fill(reverse = TRUE)) + scale_fill_viridis_d("UMI count", direction = 1, option = "plasma") + 
  labs(x = '',y = 'Fraction of cells (0-1)', title = 'Number of TCRB UMI per cell') + theme_classic() + Axis_themes + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1), plot.title = element_text(hjust = 0.5))

QC_list <- map(QC_list, function(test) test + geom_boxplot(fill = 'white',width = 0.2, outlier.shape = NA))

save_plot("TCR_QC_2.pdf", plot_grid(plotlist = QC_list[1:2]), base_height = 2, base_width = 4, unit = "in", family = "ArialMT")

#These are I think the best possible way to show them. SOme of the data is a bit less than preferable. 

# output into saved format ------------------------------------------------

QC_list <- map(QC_list, `+`, theme(plot.margin = unit(c(0, 0, 0, 0), "in")))

library(cowplot)

save_plot("TCR_QC.png", plot_grid(plotlist = QC_list), base_height = 4, unit = "in")
save_plot("TCR_QC.eps", plot_grid(plotlist = QC_list), base_height = 4, unit = "in", family = "ArialMT")

#the resulting plots are not perfect, but need illustrator touches to actually get the scales to be correct.

#need to look at what the scale is for the Seurat output.

# plotting genes for CD8 T cell subtypes ----------------------------------

Naive_genes <- c("CCR7","TCF7","SELL","CD44")
Eff_Em_genes <- c("GZMB","PRF1","ID2","GZMA","KLRG1")
Activation_exhaustion <- c("PDCD1","LAG3","TIGIT","CTLA4")

gene_list <- list(Naive_genes,Eff_Em_genes,Activation_exhaustion)
names(gene_list) <- c("Naive","Eff_em","Act_exh")

Selected_T_cell_genes <- map(gene_list,function(gene_list) FeaturePlot(E7_lim_1k, gene_list, pt.size = 0.5, do.return = TRUE))
names(Selected_T_cell_genes) <- c("Naive","Eff_em","Act_exh")

Selected_T_cell_genes <- map(Selected_T_cell_genes,function(Selected_T_cell_genes) map(Selected_T_cell_genes,`+`,TSNE_theme))

plot_grid(plotlist = Selected_T_cell_genes[["Naive"]], nrow = 1)

save_plot("Naive_genes.png", plot_grid(plotlist = Selected_T_cell_genes[["Naive"]],nrow = 1),ncol = 4, base_height = 1.2, base_aspect_ratio = 0.9, unit = "in")
save_plot("Naive_genes.eps", plot_grid(plotlist = Selected_T_cell_genes[["Naive"]],nrow = 1),ncol = 4, base_height = 1.2, base_aspect_ratio = 0.9, unit = "in", family = "ArialMT")

save_plot("Eff_Em_genes.png", plot_grid(plotlist = Selected_T_cell_genes[["Eff_em"]],nrow = 1),ncol = 5, base_height = 1.2, base_aspect_ratio = 0.9, unit = "in")
save_plot("Eff_Em_genes.eps", plot_grid(plotlist = Selected_T_cell_genes[["Eff_em"]],nrow = 1),ncol = 5, base_height = 1.2, base_aspect_ratio = 0.9, unit = "in", family = "ArialMT")

save_plot("Act_exh_genes.png", plot_grid(plotlist = Selected_T_cell_genes[["Act_exh"]],nrow = 1),ncol = 4, base_height = 1.2, base_aspect_ratio = 0.9, unit = "in")
save_plot("Act_exh_genes.eps", plot_grid(plotlist = Selected_T_cell_genes[["Act_exh"]],nrow = 1),ncol = 4, base_height = 1.2, base_aspect_ratio = 0.9, unit = "in", family = "ArialMT")

#After some discussion, I think the idea is that the selected genes are good for just a general showing of known
#genes and where they are, but doens't acutlaly directly maps expansion to their expression. Another idea is to make
#a heatmap of the gene expression by clonotypes, and order the clonotypes by the size of their expansion.


# Heatmap of genes by expanded clonotypes ---------------------------------

#first, let's get a set of cells that we want, in order
exp_ordered_cells <- Meta_data_1k %>% filter(!is.na(TRB_CDR3)) %>% arrange(desc(TRB_count),TRB_CDR3) %>% select(cell_id)

#then we want the some sort of ordered genes as well, I think just getting the cyto_v_myc genes would be fine, but we'd probably
#want to filter out for the ones that are significnat.

Cyto_v_Myc %>% filter(avg_logFC>0) %>% count(significance)

#I think there's a pretty even distribution of the two, so I think we can go ahead and make those settings.

genes_1_2 <- Cyto_v_Myc %>% filter(significance == "significant") %>% arrange(desc(avg_logFC)) %>% select(genes) %>% filter(!(genes == "CCR7"))

#CCR 7 showed up here, mostly as an artifact of the fact that there are some CCR7 expression in the Myc-activated cells,
#which could show as an intermediate between naive and activated cells

# we have not ran a DE test for Naive cells. so i think I should run it to get the genes that would differentiate between the two.
levels(E7_lim_1k@ident)

Naive_v_all <- DE_test(E7_lim_1k, group_1 = "Naive",group_2 = c("Myc_activated","Cytotoxic"))

#examine the table and see if the number of genes seem okay. 

genes_3 <- Naive_v_all %>% filter(significance == "significant", avg_logFC > 0) %>% arrange(desc(avg_logFC)) %>% 
  filter(row_number(-avg_logFC)<=30) %>% select(genes)

heatmap_genes <- c(genes_1_2$genes,genes_3$genes)

class(heatmap_genes)

#let's try to fetch data with these parameters, and see if we can ge the dataset that we want.

#we need to scale all genes in the dataset, so let's create a new seuart object here.

Heatmap_seurat <- ScaleData(E7_lim_1k, vars.to.regress = c('nUMI', 'percent.mito'), 
                            model.use = 'poisson', genes.use = heatmap_genes)



test <- FetchData(Heatmap_seurat, vars.all = heatmap_genes, cells.use = exp_ordered_cells$cell_id, use.scaled = TRUE)

#set a min and max to make sure the ranges show

min <- -1
max <- 1
test[which(test<min)] = min
test[which(test>max)] = max

#the plots don't look particularly segregated. This might not be super surprising, since there is quite a lot of 
#heterogeneity within the effector and myc- activated functions. I think A better way to look at the data might
#be trying to do a clonal Dotplot instead.


# Clonal Dotplot ----------------------------------------------------------
#Collapsed_cluster has the correct set of parameters that we want to get included in the overall
#Dot plot, but I do think maybe we can re-run that.

View(Collapsed_cluster_list)


#To do with the the correct DotPlot, we should go and and subset the correct Seurat

Dot_plot_seurat <- SetAllIdent(E7_lim_1k, id = "TRB_CDR3")
Dot_plot_seurat <- SubsetData(Dot_plot_seurat, ident.use = Collapsed_cluster_list, subset.raw = TRUE)

Dot_plot_seurat@ident <- factor(Dot_plot_seurat@ident, levels = Collapsed_cluster_list)

Dot_plot_genes <- c(Eff_Em_genes,Activation_exhaustion,Naive_genes)

Dot_plot_seurat <- Seurat_process(Dot_plot_seurat)

Dot_plot_seurat

ScaleData(E7_lim_1k, vars.to.regress = c('nUMI', 'percent.mito'), 
          model.use = 'poisson', genes.use = heatmap_genes)

test <- DotPlot(Dot_plot_seurat, genes.plot = Dot_plot_genes, cols.use = c(plasma(10)[1],plasma(10)[10]),plot.legend = TRUE, do.return = TRUE) + Axis_themes

test2 <- DotPlot(Dot_plot_seurat, genes.plot = unique(heatmap_genes), cols.use = c(plasma(10)[1],plasma(10)[10]),plot.legend = TRUE, do.return = TRUE) + Axis_themes

#I wonder if we use the module scores, it would be helpful.

Naive_genes
Eff_Em_genes
Activation_exhaustion

Dot_plot_seurat <- AddModuleScore(Dot_plot_seurat,genes.list = list(c(Naive_genes[1:3])), ctrl.size = 5, enrich.name = "Naive_Module_test")

Dot_plot_seurat <- AddModuleScore(Dot_plot_seurat,genes.list = list(c(Eff_Em_genes)), ctrl.size = 5, enrich.name = "Eff_Em")

Dot_plot_seurat <- AddModuleScore(Dot_plot_seurat,genes.list = list(c(Activation_exhaustion)), ctrl.size = 5, enrich.name = "Act_exh")

View(Dot_plot_seurat@meta.data)

Dot_plot_meta <- Dot_plot_seurat@meta.data %>% group_by(TRB_CDR3) %>% summarise(Naive = mean(Naive_Module_test1), Eff_Em = mean(Eff_Em1), Act_exh = mean(Act_exh1))

Dot_plot_meta <- Dot_plot_meta %>% mutate(TRB_CDR3 = factor(TRB_CDR3,levels = Collapsed_cluster_list)) %>% arrange(TRB_CDR3)

# Send it out as matrix.

pheatmap(as.matrix(Dot_plot_meta[2:4]), cluster_cols = FALSE, cluster_rows = FALSE)

#COnferred with Brinda, and decided that the best thing to do here is to change the dot plot with the 
#larger number of genes into a heatmap. That way we could 

# making new Z score and heatmap ------------------------------------------

expanded_cells <- Meta_data_1k %>% filter(TRB_CDR3 %in% Collapsed_cluster_list) %>% mutate(TRB_CDR3 = factor(TRB_CDR3,levels = Collapsed_cluster_list)) %>% arrange(TRB_CDR3) %>% select(cell_id, TRB_CDR3)
Heatmap_data <- FetchData(Heatmap_seurat, vars.all = heatmap_genes, cells.use = expanded_cells$cell_id)
Heatmap_data <- Heatmap_data %>% mutate(cell_id = expanded_cells$cell_id, TRB_CDR3 = expanded_cells$TRB_CDR3)
Heatmap_data <- Heatmap_data %>% select(cell_id, TRB_CDR3, everything())

Heatmap_data_collapse <- Heatmap_data %>% select(-cell_id) %>% group_by(TRB_CDR3) %>% summarise_all(funs(mean))

Heatmap_data_collapse_z <- Heatmap_data_collapse %>% mutate_if(is.double, scale)

Heatmap_data_collapse_z_matrix <- as.matrix(Heatmap_data_collapse_z[-1])
rownames(Heatmap_data_collapse_z_matrix) <-  Heatmap_data_collapse_z$TRB_CDR3

#then take scaling

min <- -2.5
max <- 2.5
Heatmap_data_collapse_z_matrix[which(Heatmap_data_collapse_z_matrix<min)] <- min
Heatmap_data_collapse_z_matrix[which(Heatmap_data_collapse_z_matrix>max)] <-  max

pheatmap(Heatmap_data_collapse_z_matrix, cluster_cols = FALSE, cluster_rows = FALSE)

#the heatmap works reasonably well. I think the next thing to do is to include some actual Naive cells on the bottom,
#that way that side of teh graphs will also clean up nicely. 

#let's try to get the singletons in the naive clusters
naive_single <- Meta_data_1k %>% filter(res.2 == 9, TRB_count<=1) %>% select(TRB_CDR3, cell_id)

#To make it more easily plotted, we'll randomly sample 10-15 of these cells, and add that to the end.

sampled_single<- naive_single %>% sample_n(15)

Heatmap_data <- FetchData(Heatmap_seurat, vars.all = heatmap_genes, cells.use = c(expanded_cells$cell_id,sampled_single$cell_id))

Heatmap_data <- as_tibble(Heatmap_data) %>% mutate(cell_id = c(expanded_cells$cell_id,sampled_single$cell_id), TRB_CDR3 = c(as.character(expanded_cells$TRB_CDR3),sampled_single$TRB_CDR3)) %>% 
  select(cell_id, TRB_CDR3, everything())

Heatmap_data <- Heatmap_data %>% mutate(TRB_CDR3 = factor(TRB_CDR3, levels = c(Collapsed_cluster_list,sampled_single$TRB_CDR3)))

Heatmap_data_collapse <- Heatmap_data %>% select(-cell_id) %>% group_by(TRB_CDR3) %>% summarise_all(funs(mean))

Heatmap_data_collapse_z <- Heatmap_data_collapse %>% mutate_if(is.double, scale)

Heatmap_data_collapse_z_matrix <- as.matrix(Heatmap_data_collapse_z[-1])
rownames(Heatmap_data_collapse_z_matrix) <-  Heatmap_data_collapse_z$TRB_CDR3

min <- -2
max <- 2
Heatmap_data_collapse_z_matrix[which(Heatmap_data_collapse_z_matrix<min)] <- min
Heatmap_data_collapse_z_matrix[which(Heatmap_data_collapse_z_matrix>max)] <-  max

pheatmap(Heatmap_data_collapse_z_matrix, cluster_cols = TRUE, cluster_rows = FALSE)

#the heatmap looks pretty good, but I think we can make it slightly more 


# continue to clean the heatmap -------------------------------------------
length(heatmap_genes)

#somehow ended up with a set of genes that tend to correlate pretty well.
#but shalek still want to try to get some known signatures, and then make modules out of them
#in order to get good markers. I think the 10.1016/j.cell.2016.08.052 paper from regev lab might be helpful
#the issue is we probably have to look at more than a few of the modules, in order to get good scoring.

#I think after talking to Chris we decided that we just have to go with the modules from the paper
#and then put clean up the color annotation

#make sample annotations

sample_row <- data.frame(Clonal_exp = rep(c("Expanded","Singleton"),c(2,2)))
rownames(sample_col) <- colnames(total_join_3_animals_matrix)
sample_col$cluster <- rep(c("cytotoxic","Myc activated"),c(5))

TR_heatmap_plot(total_join_3_animals_matrix,row_clust = FALSE,color_levels = magma(20), na_col = "grey40",
                annotation_col= sample_col)
mat_colors <- list(sample = c(gg_color_hue(8)[8],gg_color_hue(4)), cluster = Dark2_4[1:2])
names(mat_colors$sample) <- c(unique(as.character(sample_col$sample)))
names(mat_colors$cluster) <- c(unique(sample_col$cluster))



TRB_3_animal<- pheatmap(total_join_3_animals_matrix,cluster_rows = FALSE,cluster_cols = FALSE,color = magma(20), na_col = "grey40",
                        annotation_col= sample_col,border_color =NA, annotation_colors = mat_colors[1:7],
                        show_colnames = FALSE, width=8, height=2, fontsize = 6)


# importing gene lists ----------------------------------------------------
module_lists <- read_csv("NIHMS813125_cluster.csv")

module_lists <- module_lists %>% rename(gene = x) %>% mutate(cluster = factor(as.character(cluster)))
module_lists <- module_lists %>% select(gene,cluster,everything())
module_lists <- module_lists %>% mutate(cluster = fct_relevel(cluster,"10",after = Inf))
#there are three thousand genes, and they are grouped into clusters as defined in the paper.
#we should pick a couple to use 

module_vector <- levels(module_lists$cluster)
cluster_names <- str_c("cluster",module_vector)

cluster_genes <- map(module_vector, function(module_vector) module_lists %>% filter(cluster == module_vector) %>% 
      select(gene))

names(cluster_genes) <- cluster_names

#One thing that we'll have to look at is whether we can just first get an idea of how much overlap we have
#between the cluster and our dataset

All_genes <- rownames(Heatmap_seurat@data)

cluster_genes <- map(cluster_genes, function(cluster_genes) intersect(cluster_genes$gene, All_genes))

#There is still a decent amount of overlap. can continue.

#before continuing, let's just scale everything first, in case that comes up useful later.
Heatmap_seurat <- ScaleData(E7_lim_1k, vars.to.regress = c('nUMI', 'percent.mito'), 
                            model.use = 'poisson', genes.use = NULL)

#then let's add the modules, then let them calclulate teh scores.

Heatmap_seurat <- AddModuleScore(Heatmap_seurat,genes.list = cluster_genes, ctrl.size = 10, enrich.name = "Cluster")

View(Heatmap_seurat@meta.data)
#seems like the scores have been added on there correctly. Now let's try to pull the TRB list that were in figure 3
# and then 

Heatmap_module <- Heatmap_seurat@meta.data %>% filter(TRB_CDR3 %in% c(as.character(expanded_cells$TRB_CDR3),sampled_single$TRB_CDR3)) %>% 
  mutate(TRB_CDR3 = factor(TRB_CDR3, levels = c(Collapsed_cluster_list,sampled_single$TRB_CDR3))) %>% arrange(TRB_CDR3) %>% 
  select(cell_id, TRB_CDR3,starts_with("Cluster"))


Heatmap_module_collapse <- Heatmap_module %>% select(-cell_id,-cluster) %>% group_by(TRB_CDR3) %>% 
  summarise_all(funs(mean)) %>% mutate_if(is.double, scale)

View(Heatmap_module_collapse)

Heatmap_module_collapse_Matrix <- as.matrix(Heatmap_module_collapse[-1])
rownames(Heatmap_module_collapse_Matrix) <-  Heatmap_module_collapse$TRB_CDR3

min <- -2.5
max <- 2.5
Heatmap_module_collapse_Matrix[which(Heatmap_module_collapse_Matrix<min)] <- min
Heatmap_module_collapse_Matrix[which(Heatmap_module_collapse_Matrix>max)] <-  max

pheatmap(Heatmap_module_collapse_Matrix, cluster_cols = TRUE, cluster_rows = FALSE)

#The heatmap here looks pretty good as well. I think in the end, we'll have to show both. so let's focus on cleaning
#up a version that has the correct annotations

expanded_length <- length(Collapsed_cluster_list)
singleton_length <- length(sampled_single$TRB_CDR3)

sample_row <- data.frame(Clonal_exp = rep(c("Expanded","Singleton"),c(expanded_length,singleton_length)))
rownames(sample_row) <- rownames(Heatmap_module_collapse_Matrix)
#sample_col$cluster <- rep(c("cytotoxic","Myc activated"),c(5))

col_colors <- list(Clonal_exp = c(viridis(10)[1],viridis(10)[10]))
names(col_colors$Clonal_exp) <- c(unique(as.character(sample_row$Clonal_exp)))



#now replot the heatmap with color row annotation
pheatmap(Heatmap_module_collapse_Matrix, cluster_cols = TRUE, cluster_rows = FALSE,
         annotation_row= sample_row, border_color =NA, annotation_colors = col_colors,
         show_colnames = TRUE, width=3.2, height=6, fontsize = 6, filename= "module_heatmap.pdf")
# t cell expansion elements -------------------------------------------)----

#one of the review comments is that we need to add more clarification around if the dark color dots on the tsne plot of TRB expansion
#are correlating to multipel clones or just one clone. I think one way to do that is a dotplot of expansion on x axis, and then count on
#y axis.

View(Meta_data_1k)

ggplot(Meta_data_1k %>% filter(TRB_count>1) %>% distinct(TRB_CDR3,TRB_count), aes(x =TRB_count)) + geom_dotplot(method = 'histodot', binwidth = 1, aes(fill = factor(log2(TRB_count)))) + 
  scale_fill_viridis(direction = -1,discrete = TRUE) + ylim(0,0.12)

#this funciton to generate the x-y values manually. the geom_dotplot function just doesn't do what I need it do, given the
#specific scale and colors that I would need.
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

temp <- Meta_data_1k %>% filter(TRB_count>=1) %>% distinct(TRB_CDR3,TRB_count)
temp <- temp %>% count(TRB_count)

breaks <- 2^(1:9)

my_breaks_color <- log2(2^(0:9))
my_breaks_label_color <- 2^(c(0,seq(9)))

ggplot_dot(x = temp$TRB_count, count = temp$n, pch = 21, color = "grey87") + 
  scale_x_log10(breaks = breaks, limits = c(2,512)) + scale_fill_viridis(direction = -1,breaks = my_breaks_color, 
                                                                         labels = my_breaks_label_color,limits = c(0,log2(512))) + 
  labs(x = "TRB Clonal size", y = "Count of unique clones", fill = "TRB Clonal size") + Axis_themes + ylim(0,35)

ggsave("dotplot.png",plot = last_plot(), device = "png",width = 8, height = 2, units = "in")
ggsave("TRB_expan_dotplot.eps",plot = last_plot(), device = "eps",width = 8, height = 2, units = "in", family = "ArialMT")


#so far looks pretty decent, as long as we , but need to make sure we clean up the labels.


# addressing PCA clusters versus the hand-defined ones --------------------

#Q2 also ask the question regarding how exactly are we defining the difference in the PCA based clusters,
#and the manual defined ones between myc and cytotoxic ones.
#I think in the text we'll have to talk about a couple of the ones we generate with regards
#to the overall clusters being called

#the basic plots are done in E7_1k_TSNE plots, and saved as such.

#Another thing taht I think would correlate well, might be to break the clusters by expansion level
#but also just where the clonotypes fall on all these clsuters.

#clonotype to res.2 clusters

res.2_TR_hmp <- TR_cluster_heatmap(seurat_meta = Meta_data_1k %>% filter(TRB_CDR3 %in% Collapsed_cluster_list), cluster_var = res.2)
res.2_TR_hmp <- TR_heatmap_spread(res.2_TR_hmp,cluster_var = res.2)
TR_heatmap_plot(res.2_TR_hmp,name = "TRB_heatmap_test.png")

pheatmap(res.2_TR_hmp)
#the results confirm that cluster 4 is the one with the clearest separation of speicific clones
#and I think when we sort that out a bit more, we'll be able to get a better looking heatmap

res.2_cluster_cols <- hclust(dist(t(res.2_TR_hmp)))
plot(res.2_cluster_cols, main = "Unsorted Dendrogram", xlab = "", sub = "")

library(dendsort)

sort_hclust <- function(...) as.hclust(dendsort(as.dendrogram(...)))

res.2_cluster_cols <- sort_hclust(res.2_cluster_cols)
plot(res.2_cluster_cols, main = "sorted Dendrogram", xlab = "", sub = "")



pheatmap(res.2_TR_hmp, cluster_cols = res.2_cluster_cols)
res.2_cluster_rows <- hclust(dist(res.2_TR_hmp))

plot(res.2_cluster_rows, main = "Unsorted Dendrogram", xlab = "", sub = "")
res.2_cluster_rows <- sort_hclust(res.2_cluster_rows)

pheatmap(res.2_TR_hmp, cluster_cols = res.2_cluster_cols,cluster_rows = res.2_cluster_rows)
#let's first clean up by cleaning the clonotypes so they overlap with the other graphs that we have shown before

#then let's clean up the color, and see if we can use the classic pHE
require(grid)
pdf(file = "cluster_v_clonotypes.pdf" ,width = 6, height = 6)
setHook("grid.newpage", function() pushViewport(viewport(x=1,y=1,width=0.9, height=0.9, name="vp", just=c("right","top"))), action="prepend")
pheatmap(res.2_TR_hmp, color = cividis(20), border_color =NA, drop_levels = TRUE, fontsize = 6, cluster_rows = res.2_cluster_rows, cluster_cols = res.2_cluster_cols)
setHook("grid.newpage", NULL, "replace")
grid.text("cluster", y=-0.07, x = 0.35, gp=gpar(fontsize=16))
grid.text("TRB clones", x=-0.07, rot=90, gp=gpar(fontsize=16))
dev.off()

#Then let's also make a plot of average expansion number for each cluster.

View(Meta_data_1k)
temp <- Meta_data_1k %>% select(TRB_CDR3, res.2,TRB_count) %>% filter(!is.na(TRB_CDR3)) %>% distinct()

ggplot(temp, aes(x = res.2, y =  log2(TRB_count), fill = res.2)) + geom_violin(scale = 'width') + theme_classic() +
  Axis_themes + labs(x = 'cluster',y = 'log2(TRB clonal size)', fill = 'cluster') + 
  scale_y_continuous(breaks = my_breaks_color, labels = my_breaks_label_color)
  
ggsave("violin_clonalsize.png",plot = last_plot(), device = "png",width = 4, height = 3, units = "in")
ggsave("violin_clonalsize.eps",plot = last_plot(), device = "eps",width = 4, height = 3, units = "in", family = "ArialMT")

#We should try to make a 3 plot 
ggsave("TSNE_TRB_expan.png",plot = TSNE_TRB_expan + theme(legend.position = "none"), device = "png",width = 2.8, height = 3.2, units = "in")



E7s_TSNE_list <- list(TSNE_TRB_expan, E7_1k_TSNE[[1]], collapsed_identiy_TSNE)

E7s_TSNE_list <- map(E7s_TSNE_list,`+`,No_axis_labels + No_legend)

plot_grid(plotlist = E7s_TSNE_list, nrow = 1) 

save_plot("E7s_TSNE_list.png", plot_grid(plotlist = E7s_TSNE_list, nrow = 1), base_height = 2.2, base_width = 5.9, unit = "in")
save_plot("E7s_TSNE_list.eps", plot_grid(plotlist = E7s_TSNE_list, nrow = 1), base_height = 2.2, base_width = 5.9, unit = "in",family = "ArialMT")


# Redoing figure 3 panel D ------------------------------------------------
#I'm getting very tired of writing and re-doing everything 200 times. I just don't think Chris
#understands the dfificulty in explaining the figure, and actually making it cleanly flow.
#adding all the reviewer comment is making it much more difficult to actually write the paper.

#have to take the list of clones, and re-do where the plots.

clone_list <- c(Collapsed_cluster_list,sampled_single$TRB_CDR3)

Figure3d <- TR_cluster_heatmap(Meta_data_1k %>% filter(TRB_CDR3 %in% clone_list), cluster_var = res.3, clone_cutoff = 1)
Figure3d <- TR_heatmap_spread(Figure3d,cluster_var = res.3)

Figure3d <- Figure3d[order(match(rownames(Figure3d),clone_list)), ]

sample_col_2 <- data.frame(cluster = c("cytotoxic","Myc activated","Naive","HSP High"))
rownames(sample_col_2) <- colnames(Figure3d)

col_colors_2 <- list(cluster = Dark2_4)
names(col_colors_2$cluster) <- c(as.character(sample_col_2$cluster))

pheatmap(Figure3d,color = magma(20), annotation_row= sample_row, annotation_colors = c(col_colors,col_colors_2),
         annotation_col = sample_col_2,border_color =NA, drop_levels = TRUE, 
         cluster_rows = FALSE, cluster_cols = FALSE,width=3.2, height=6, fontsize = 6, filename= "cluster_heatmap_updated.pdf")

# I think also redo the gene plot per cluster -----------------------------

heatmap_genes <- c(genes_1_2$genes[!(genes_1_2$genes == "TCF7")],genes_3$genes)

pheatmap(Heatmap_data_collapse_z_matrix, cluster_cols = TRUE, cluster_rows = FALSE,
         border_color =NA,width=7.5, height=6, fontsize = 6, filename= "cluster_heatmap_updated_2.pdf")


# Update notes ------------------------------------------------------------

#It's possible that people might ask to see a separation of genes. I wonder if I can just provide that as
#a text file, instead of heatmap, which would be much easier.

#then I think Shalek still wants to see clonotype matched DE. which is still going to take some time.

