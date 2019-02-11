# Import required packages

library(ggplot2)
library(monocle)

#--------------------------------------------------------------------------------------
# Setup the working directory and loading r object.

setwd('/home/parashar/Data/scRNASeq_workshop/data/neuron_1k_v2')
getwd()

load('neuro.v2.Robj')

#--------------------------------------------------------------------------------------
# Setting up Monocle CellDatSet

mono.neuro.v2 <-  importCDS(otherCDS = neuro.v2)
mono.neuro.v2 <- estimateSizeFactors(mono.neuro.v2)
mono.neuro.v2 <- estimateDispersions(mono.neuro.v2)

#--------------------------------------------------------------------------------------
# Setting up genes for reducing dimensions using most variable genes identified by Seurat

mono.neuro.v2 <- setOrderingFilter(mono.neuro.v2, neuro.v2@var.genes)
mono.neuro.v2 <- reduceDimension(mono.neuro.v2, max_components = 2,
                                 reduction_method = 'DDRTree', verbose=TRUE)
mono.neuro.v2 <- orderCells(mono.neuro.v2)

#--------------------------------------------------------------------------------------
# Trajectory visualization

plot_cell_trajectory(mono.neuro.v2)
plot_cell_trajectory(mono.neuro.v2, color_by = 'res.0.8')
plot_cell_trajectory(mono.neuro.v2, color_by = "res.0.8") + facet_wrap(~State, nrow = 2)

plot_cell_trajectory(mono.neuro.v2, color_by = "Pseudotime")

top_marker <- markers %>% group_by(cluster) %>% top_n(1, avg_logFC)
plot_genes_in_pseudotime(mono.neuro.v2[top_marker$gene], color_by = "res.0.8")
