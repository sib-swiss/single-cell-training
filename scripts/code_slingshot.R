library(SingleCellExperiment)
library(BiocSingular)
library(slingshot)
library(ggplot2)
library(ggbeeswarm)
set.seed(3)
# Load the SCE object
deng_SCE <- readRDS("data/deng_dataset/deng-reads.rds")

# Perform the first steps of the analysis
# The deng_SCE object contains cells that were isolated at different stages of
# mouse embryogenesis, from the zygote stage to the late blastula.

# What class is the deng_SCE object, and how is it organized?
class(deng_SCE)

structure(deng_SCE)
# How many mouse cells are at each stage?
table(deng_SCE$cell_type1)
table(deng_SCE$cell_type2)
# Re-order the levels of the factor storing the cell developmental stage.
deng_SCE$cell_type2 <- factor(deng_SCE$cell_type2,
                              levels = c("zy", "early2cell", "mid2cell", "late2cell",
                                         "4cell", "8cell", "16cell", "earlyblast", "midblast",
                                         "lateblast"))

# Run PCA on Deng data. Use the runPCA function from the SingleCellExperiment package.
deng_SCE <- runPCA(deng_SCE, ncomponents = 50)

# Use the reducedDim function to access the PCA and store the results.
pca <- SingleCellExperiment::reducedDim(deng_SCE, "PCA")

# Describe how the PCA is stored in a matrix. Why does it have this structure?
head(pca)

# Add PCA data to the deng_SCE object.
deng_SCE$PC1 <- pca[, 1]
deng_SCE$PC2 <- pca[, 2]

# Plot PC biplot with cells colored by cell_type2.
# colData(deng_SCE) accesses the cell metadata DataFrame object for deng_SCE.
# Look at Figure 1A of the paper as a comparison to your PC biplot.
ggplot(as.data.frame(colData(deng_SCE)), aes(x = PC1, y = PC2, color = cell_type2)) +
  geom_point(size=2, shape=20) +
  theme_classic() +
  xlab("PC1") + ylab("PC2") + ggtitle("PC biplot")

# PCA is a simple approach and can be good to compare to more complex algorithms
# designed to capture differentiation processes. As a simple measure of pseudotime
# we can use the coordinates of PC1.
# Plot PC1 vs cell_type2.
deng_SCE$pseudotime_PC1 <- rank(deng_SCE$PC1)  # rank cells by their PC1 score

# Create a jitter plot
# library(ggbeeswarm)
ggplot(as.data.frame(colData(deng_SCE)), aes(x = pseudotime_PC1, y = cell_type2,
                                             colour = cell_type2)) +
  geom_quasirandom(groupOnX = FALSE) +
  theme_classic() +
  xlab("PC1") + ylab("Timepoint") +
  ggtitle("Cells ordered by first principal component")

# Read the Slingshot documentation (?slingshot) and then run Slingshot below.
# Given your understanding of the algorithm and the documentation, what is one
# major set of parameters we omitted here when running Slingshot?
sce <- slingshot(deng_SCE, reducedDim = 'PCA')  # no clusters
# Check how the slingshot object looks like
SlingshotDataSet(sce)

# Plot PC1 vs PC2 colored by Slingshot pseudotime.
colors <- rainbow(50, alpha = 1)
plot(reducedDims(sce)$PCA, col = colors[cut(sce$slingPseudotime_1,breaks=50)], pch=16, asp = 1)
lines(SlingshotDataSet(sce), lwd=2)

# Plot Slingshot pseudotime vs cell stage.
ggplot(as.data.frame(colData(deng_SCE)), aes(x = sce$slingPseudotime_1, y = cell_type2,
                                             colour = cell_type2)) +
  geom_quasirandom(groupOnX = FALSE) +
  theme_classic() +
  xlab("Slingshot pseudotime") + ylab("Timepoint") +
  ggtitle("Cells ordered by Slingshot pseudotime")

# Cluster cells using the Seurat workflow below.
gcdata <- CreateSeuratObject(counts = counts(deng_SCE), project = "slingshot")
gcdata <- NormalizeData(object = gcdata, normalization.method = "LogNormalize",
                        scale.factor = 10000)
gcdata <- FindVariableFeatures(object = gcdata, mean.function = ExpMean, dispersion.function = LogVMR)
topgenes <- VariableFeatures(gcdata)

gcdata <- ScaleData(object = gcdata, do.center = T, do.scale = F)
gcdata <- RunPCA(object = gcdata, pc.genes = gcdata@var.genes, do.print = TRUE, pcs.print = 1:5,
                 genes.print = 5)
gcdata <- FindNeighbors(gcdata, reduction="pca", dims = 1:5)

gcdata <- FindClusters(object = gcdata,
                       resolution = 0.6)

# Add clustering information from Seurat to the deng_SCE object
# Then run Slingshot using these cluster assignments.
colData(deng_SCE)$Seurat_clusters <- as.character(Idents(gcdata))  # go from factor to character
deng_SCE <- slingshot(deng_SCE, clusterLabels = 'Seurat_clusters', reducedDim = 'PCA')
?slingshot
head(colData(deng_SCE))
# Check how the slingshot object has evolved
SlingshotDataSet(deng_SCE)


# Plot PC1 vs PC2 colored by Slingshot pseudotime (or pseudolines).
colors <- rainbow(50, alpha = 1)
plot(reducedDims(deng_SCE)$PCA, col = colors[cut(deng_SCE$slingPseudotime_2,breaks=50)], pch=16, asp = 1)
lines(SlingshotDataSet(deng_SCE), lwd=2)

# Plot Slingshot pseudotime vs cell stage.
ggplot(data.frame(cell_type2 = deng_SCE$cell_type2, slingPseudotime_1 = deng_SCE$slingPseudotime_1), aes(x = slingPseudotime_1, y = cell_type2,
                                             colour = cell_type2)) +
  geom_quasirandom(groupOnX = FALSE) +
  theme_classic() +
  xlab("Slingshot pseudotime") + ylab("Timepoint") +
  ggtitle("Cells ordered by Slingshot pseudotime")


ggplot(data.frame(cell_type2 = deng_SCE$cell_type2, slingPseudotime_2 = deng_SCE$slingPseudotime_2), aes(x = slingPseudotime_2, y = cell_type2,
                                             colour = cell_type2)) +
  geom_quasirandom(groupOnX = FALSE) +
  theme_classic() +
  xlab("Slingshot pseudotime") + ylab("Timepoint") +
  ggtitle("Cells ordered by Slingshot pseudotime")


dim(colData(deng_SCE))


# Plot PC1 vs PC2 colored by Slingshot lineages, with color corresponding to seurat clusters

plot(reducedDim(deng_SCE), col = colData(deng_SCE)$Seurat_clusters, pch = 16, cex = 1)
lines(SlingshotDataSet(deng_SCE), lwd = 2, type = 'lineages', col = 'black')

# Plot PC1 vs PC2 colored by Slingshot lineages, with color corresponding to celltypes
plot(reducedDim(deng_SCE), col = colData(deng_SCE)$cell_type2, pch = 16, cex = 1)
lines(SlingshotDataSet(deng_SCE), lwd = 2, type = 'lineages', col = 'black')
legend("bottomright",col=1:length(colData(deng_SCE)$cell_type2),legend = levels(factor(colData(deng_SCE)$cell_type2)),pch=20)

# Save current progress.
# save(deng_SCE, file = Rda.slingshot.path)
# To load the data, run the following command.
# load(Rda.slingshot.path)
