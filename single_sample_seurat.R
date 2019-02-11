# Import required packages

library(ggplot2)
library(Seurat)
library(dplyr)

#--------------------------------------------------------------------------------------
# Setup the working directory

setwd('/home/parashar/Data/scRNASeq_workshop/data/neuron_1k_v2')
getwd()

#--------------------------------------------------------------------------------------
# Load 10x data

data.tenx.v2 <- Read10X(data.dir = 'seurat_matrix')
neuro.v2 <- CreateSeuratObject(raw.data = data.tenx.v2, min.cells = 10, min.genes = 200, 
                               project = "10x_v2")

#--------------------------------------------------------------------------------------
# Summarize data

neuro.v2
VlnPlot(object = neuro.v2, features.plot = c("nGene", "nUMI"), nCol = 2)

#--------------------------------------------------------------------------------------
# Fetching all gene names and mitochondrial gene names

gene.names <- rownames(neuro.v2@data)

mito.genes <- grep(pattern = "^MT-", x = rownames(neuro.v2@data), value = TRUE)
mito.genes

mito.genes <- grep(pattern = "^mt-", x=gene.names, value = TRUE)
mito.genes
length(x=mito.genes)

#--------------------------------------------------------------------------------------
# Calculate percentage of mitochondrial UMIs per cell

mito.umi.per.cell <- Matrix::colSums(neuro.v2@raw.data[mito.genes, ])
head(mito.umi.per.cell)
total.umi.per.cell <-Matrix::colSums(neuro.v2@raw.data)
head(total.umi.per.cell)
percent.mito <- 100*mito.umi.per.cell/total.umi.per.cell
head(percent.mito)

#--------------------------------------------------------------------------------------
# Calculate percentage of ribosomal UMIs per cell

ribo.small.genes <- grep(pattern = "^Rps", x = rownames(neuro.v2@data), value = TRUE)
ribo.small.genes
ribo.large.genes <- grep(pattern = "^Rpl", x = rownames(neuro.v2@data), value = TRUE)
ribo.large.genes

ribo.genes <- c(ribo.small.genes, ribo.large.genes)
ribo.genes
length(ribo.genes)

ribo.umi.per.cell <- Matrix::colSums(neuro.v2@raw.data[ribo.genes, ])
head(ribo.umi.per.cell)
percent.ribo <- 100*ribo.umi.per.cell/total.umi.per.cell
head(percent.ribo)

#--------------------------------------------------------------------------------------
# Add mitochondrial and ribosomal percentages as metadata to the Seurat object

head(neuro.v2@meta.data)

neuro.v2 <- AddMetaData(object = neuro.v2, metadata = percent.mito, col.name = "percent.mito")
neuro.v2 <- AddMetaData(object = neuro.v2, metadata = percent.ribo, col.name = "percent.ribo")

head(neuro.v2@meta.data)

VlnPlot(object = neuro.v2, features.plot = c("nGene", "nUMI", "percent.mito", "percent.ribo"), nCol = 4)

#--------------------------------------------------------------------------------------
# Filtering poor quality cells

GenePlot(object = neuro.v2, gene1 = "nUMI", gene2 = "nGene")

ggplot(neuro.v2@meta.data, aes(x=nUMI, y=nGene)) + geom_point()
ggplot(neuro.v2@meta.data, aes(x=nGene, y=percent.mito)) + geom_point()
ggplot(neuro.v2@meta.data, aes(x=nGene, y=percent.ribo)) + geom_point()

neuro.v2 <- FilterCells(object = neuro.v2, subset.names = c("nGene", "nUMI", "percent.mito"), 
                        low.thresholds = c(200, 500, -Inf), high.thresholds = c(4500, 20000, 20))
neuro.v2
ggplot(neuro.v2@meta.data, aes(x=nUMI, y=nGene)) + geom_point()

#--------------------------------------------------------------------------------------
# Additional filtering

filter.test1 <- neuro.v2@meta.data$nGene < 1000
head(filter.test1)
table(filter.test1)

filter.test2 <- neuro.v2@meta.data$nUMI > 15000
head(filter.test2)
table(filter.test2)

ggplot(neuro.v2@meta.data, aes(x=nUMI)) + geom_histogram(binwidth = 500, color="black", fill="grey")

wierd.cells <- filter.test1 & filter.test2
table(wierd.cells)
neuro.v2@meta.data[wierd.cells, ]

neuro.v2@meta.data['wierd'] <- wierd.cells
ggplot(neuro.v2@meta.data, aes(x=nUMI, y=nGene)) + geom_point(aes(color=wierd))

ggplot(neuro.v2@meta.data, aes(x=percent.mito, y=nGene)) + geom_point(aes(color=wierd))
ggplot(neuro.v2@meta.data, aes(x=percent.ribo, y=nGene)) + geom_point(aes(color=wierd))

remw.neuro.v2 <- SubsetData(neuro.v2, cells.use = rownames(neuro.v2@meta.data[!wierd.cells, ]))
remw.neuro.v2
neuro.v2

#--------------------------------------------------------------------------------------
# Data Normalization

neuro.v2 <- NormalizeData(object = neuro.v2, normalization.method = "LogNormalize", 
                          scale.factor = 1000)

#--------------------------------------------------------------------------------------
# Identification of highly variable genes

neuro.v2 <- FindVariableGenes(object = neuro.v2, mean.function = ExpMean, dispersion.function = LogVMR, 
                              x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
length(x = neuro.v2@var.genes)
neuro.v2 <- FindVariableGenes(object = neuro.v2, mean.function = ExpMean, dispersion.function = LogVMR, 
                              x.low.cutoff = 0.1, x.high.cutoff = 2, y.cutoff = 1)
length(x = neuro.v2@var.genes)

#--------------------------------------------------------------------------------------
# (Optional step): Some wierd genes are highly variable, remove them

gene.names <- rownames(neuro.v2@data)
length(gene.names)

genes.exclude.idx <- grep(pattern = "^Rps|^Rpl|^mt-", x = gene.names)
genes.exclude.idx
length(genes.exclude.idx)
gene.names[genes.exclude.idx]

genes.keep <- gene.names[-genes.exclude]
length(genes.keep)

cells.keep <- colnames(neuro.v2@data)
length(cells.keep)

head(neuro.v2@meta.data)

neuro.v2 <- CreateSeuratObject(neuro.v2@raw.data[genes.keep, cells.keep], project = "10x_v2")
neuro.v2@meta.data['wierd'] <- neuro.v2@meta.data$nGene < 1000 & neuro.v2@meta.data$nUMI > 15000
neuro.v2
head(neuro.v2@meta.data)

neuro.v2 <- NormalizeData(object = neuro.v2, normalization.method = "LogNormalize", 
                          scale.factor = 1000)
neuro.v2 <- FindVariableGenes(object = neuro.v2, mean.function = ExpMean, dispersion.function = LogVMR, 
                              x.low.cutoff = 0.1, x.high.cutoff = 2, y.cutoff = 1)
length(x = neuro.v2@var.genes)

#--------------------------------------------------------------------------------------
# Cell cycle scoring

cc.genes <- readLines(con = "../regev_lab_cell_cycle_genes.txt")
length(cc.genes)
cc.genes

letters.first <- substring(cc.genes, 1, 1)
letters.first
letters.rest <- substring(cc.genes, 2, 100)
letters.rest
cc.genes <- paste0(letters.first, tolower(letters.rest))
cc.genes

s.genes <- cc.genes[1:43]
s.genes
g2m.genes <- cc.genes[44:97]
g2m.genes

neuro.v2 <- CellCycleScoring(object = neuro.v2, s.genes = s.genes, g2m.genes = g2m.genes)
head(neuro.v2@meta.data)

table(neuro.v2@meta.data$Phase)
ggplot(neuro.v2@meta.data, aes(x=S.Score)) + geom_histogram(binwidth = 0.01, color="black", fill="grey")
ggplot(neuro.v2@meta.data, aes(x=G2M.Score)) + geom_histogram(binwidth = 0.01, color="black", fill="grey")

ggplot(neuro.v2@meta.data, aes(x=S.Score, y=G2M.Score)) + geom_point(aes(color=Phase))

neuro.v2@meta.data$CC.Difference <- neuro.v2@meta.data$S.Score - neuro.v2@meta.data$G2M.Score
ggplot(neuro.v2@meta.data, aes(x=S.Score, y=G2M.Score)) + geom_point(aes(color=CC.Difference)) + scale_colour_gradient2()

#--------------------------------------------------------------------------------------
# Data Scaling

pbmc <- ScaleData(object = neuro.v2, vars.to.regress = c("nGene"))
#pbmc <- ScaleData(object = neuro.v2, vars.to.regress = c("nGene", "S.Score", "G2M.Score"))
#pbmc <- ScaleData(object = neuro.v2, vars.to.regress = c("nGene", "S.Score", "G2M.Score"))
pbmc <- ScaleData(object = neuro.v2, vars.to.regress = c("nGene"))





