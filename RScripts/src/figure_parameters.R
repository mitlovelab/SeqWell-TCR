#The purpose of this script is to include functions/codes for common parameters
#that outputs useful figures, that are well configured, so we minimized the amount of time
#playing around with resizing and illustrator problems

#this is for small TSNE plots. For large TSNE/UMAP plots, we would format using axis themes

require(ggplot2)

TSNE_theme <- theme_bw()+theme(axis.text.y = element_blank(), 
                               axis.text.x = element_blank(), 
                               axis.ticks.x= element_blank(),
                               axis.ticks.y= element_blank(),
                               axis.title.x= element_blank(),
                               axis.title.y= element_blank(),
                               strip.text.x = element_text(size = 8),
                               panel.grid.major = element_blank(),
                               panel.grid.minor = element_blank(),
                               panel.border = element_blank(),
                               panel.background = element_rect(colour = "black", fill = NA,size=1),
                               legend.position = "none",
                               plot.title = element_text(hjust = 0.5, size = 8))

Axis_themes <- theme(plot.title = element_text(size = 8),
                     axis.title = element_text(size = 8), 
                     axis.text = element_text(size = 6),
                     legend.text = element_text(size =6),
                     legend.title = element_text(size = 8),
                     strip.text.x = element_text(size = 8))

#for sizing, letter is 8.5 x 11 inch. so really we can work within 8x10 parameter
#288 pts is about 4 inches. 72 pt/inch

#primarily use ggsave to save our figures into different file format. need plot_list/cowplot to do last_plot() before
#we can actually push the figures into the devices.

#There are a lot of font inconsistency. From here on out, just use arial font. However, in saving for EPS file,
#we need to specify Family = "ArialMT" in order to make sure that the font turns out correct. The default font
#output for eps is Helvetica. But for pdf, arial will work. ggplot already plots automatically in Arial



