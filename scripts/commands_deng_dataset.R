deng_SCE <- readRDS("data/deng_dataset/deng-reads.rds")

deng_SCE$cell_type2 <- factor(deng_SCE$cell_type2,
                              levels = c("zy",
                                         "early2cell",
                                         "mid2cell",
                                         "late2cell",
                                         "4cell",
                                         "8cell",
                                         "16cell",
                                         "earlyblast",
                                         "midblast",
                                         "lateblast"))

deng_SCE <- BiocSingular::runPCA(deng_SCE, ncomponents = 50)

sce <- slingshot::slingshot(deng_SCE, reducedDim = 'PCA')

PCAplot_slingshot <- function(sce, draw_lines = TRUE, variable = NULL, legend = FALSE, ...){
  # set palette for factorial variables
  palf <- colorRampPalette(RColorBrewer::brewer.pal(8, "Set2"))
  # set palette for numeric variables
  paln <- colorRampPalette(RColorBrewer::brewer.pal(9, "Blues"))
  # extract pca from SingleCellExperiment object
  pca <- SingleCellExperiment::reducedDims(sce)$PCA
  
  if(is.null(variable)){
    col <- "black"
  }
  if(is.character(variable)){
    variable <- as.factor(variable)
  }
  if(is.factor(variable)){
    colpal <- palf(length(levels(variable)))
    colors <- colpal[variable]
  }
  if(is.numeric(variable)){
    colpal <- paln(50)
    colors <- colpal[cut(variable,breaks=50)]
  }
  
  # draw the plot
  plot(pca, bg = colors, pch = 21)
  # draw lines
  if(draw_lines){
    slingshot::lines(slingshot::SlingshotDataSet(sce), lwd = 2, ... )
  }
  # add legend
  if(legend & is.factor(variable)){
    legend("bottomright", pt.bg = colpal,legend = levels(variable),pch=21)
    
  }
}

PCAplot_slingshot(sce, variable = sce$slingPseudotime_1, draw_lines = TRUE)

gcdata <- Seurat::CreateSeuratObject(counts = SingleCellExperiment::counts(deng_SCE),
                                     project = "slingshot")

gcdata <- Seurat::NormalizeData(object = gcdata,
                                normalization.method = "LogNormalize",
                                scale.factor = 10000)

gcdata <- Seurat::FindVariableFeatures(object = gcdata,
                                       mean.function = ExpMean,
                                       dispersion.function = LogVMR)

gcdata <- Seurat::ScaleData(object = gcdata,
                            do.center = T,
                            do.scale = F)

gcdata <- Seurat::RunPCA(object = gcdata,
                         pc.genes = gcdata@var.genes)

gcdata <- Seurat::FindNeighbors(gcdata,
                                reduction="pca",
                                dims = 1:5)

# clustering with resolution of 0.6
gcdata <- Seurat::FindClusters(object = gcdata,
                               resolution = 0.6)

gcdata <- Seurat::RunUMAP(gcdata, dims = 1:5)
png("umap_gcdata.png")
Seurat::DimPlot(gcdata)
dev.off()

deng_SCE$Seurat_clusters <- as.character(Idents(gcdata))  # go from factor to character

deng_SCE <- slingshot::slingshot(deng_SCE,
                                 clusterLabels = 'Seurat_clusters',
                                 reducedDim = 'PCA')

PCAplot_slingshot(deng_SCE, variable = deng_SCE$slingPseudotime_2)

PCAplot_slingshot(deng_SCE,
                  variable = deng_SCE$Seurat_clusters,
                  type = 'lineages',
                  col = 'black',
                  legend = TRUE)

PCAplot_slingshot(deng_SCE,
                  variable = deng_SCE$cell_type2,
                  type = 'lineages',
                  col = 'black',
                  legend = TRUE)
