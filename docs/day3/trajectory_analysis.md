## Learning outcomes

**After having completed this chapter you will be able to:**

## Material

[:fontawesome-solid-file-pdf: Download the presentation](../assets/pdf/sequencing_technologies.pdf){: .md-button }

## Exercises

### Trajectory analysis using Slingshot

> :material-transit-connection-variant: This part uses the `Deng` dataset

!!! bug
    Paper citation needed

Read in data. It is an object of class `SingleCellExperiment`.

```R
deng_SCE <- readRDS("data/deng-reads.rds")
```

Perform the first steps of the analysis. The deng_SCE object contains cells that were isolated at different stages of mouse embryogenesis, from the zygote stage to the late blastula.

The levels of the cell type are in alphabetical order. We now change the level order for plotting in cell cycle order:

```R
deng_SCE$cell_type2 <- factor(deng_SCE$cell_type2,
                              levels = c("zy", "early2cell", "mid2cell", "late2cell",
                                         "4cell", "8cell", "16cell", "earlyblast", "midblast",
                                         "lateblast"))
```

We can run a PCA directly on the object of class `SingleCellExperiment` with the function `runPCA`:

```R
deng_SCE <- BiocSingular::runPCA(deng_SCE, ncomponents = 50)
```

Read the Slingshot documentation (`?slingshot::slingshot`) and then run Slingshot below.

```R
sce <- slingshot::slingshot(deng_SCE, reducedDim = 'PCA')
```

**Exercise:** Given your understanding of the algorithm and the documentation, what is one
major set of parameters we omitted here when running Slingshot?

??? done "Answer"
    We didn't set the parameter `clusterLabels`

Here is a custom function to plot the PCA based on a `slingshot` object. Run it in the console to add it to your global environment:

```R
PCAplot_slingshot <- function(sce, draw_lines = TRUE, variable = NULL, legend = FALSE, ...){
  palf <- colorRampPalette(RColorBrewer::brewer.pal(8, "Set2"))
  paln <- colorRampPalette(RColorBrewer::brewer.pal(9, "Blues"))
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
  plot(pca, bg = colors, pch = 21)
  if(draw_lines){
    slingshot::lines(slingshot::SlingshotDataSet(sce), lwd = 2, ... )
  }
  if(legend & is.factor(variable)){
    legend("bottomright", pt.bg = colpal,legend = levels(variable),pch=21)

  }
}
```

Have a look at the PCA with the slingshot pseudotime line:

```R
PCAplot_slingshot(sce, variable = sce$slingPseudotime_1, draw_lines = TRUE)
```

Also have a look at pseudotime versus cell type:

```R
ggplot(as.data.frame(colData(deng_SCE)), aes(x = sce$slingPseudotime_1, y = cell_type2,
                                             colour = cell_type2)) +
  ggbeeswarm::geom_quasirandom(groupOnX = FALSE) +
  theme_classic() +
  xlab("Slingshot pseudotime") + ylab("Timepoint") +
  ggtitle("Cells ordered by Slingshot pseudotime")
```

This already looks pretty good. Let's see whether we can improve it. First we generate clusters by using `Seurat`:

```R
gcdata <- Seurat::CreateSeuratObject(counts = counts(deng_SCE), project = "slingshot")
gcdata <- Seurat::NormalizeData(object = gcdata, normalization.method = "LogNormalize",
                        scale.factor = 10000)
gcdata <- Seurat::FindVariableFeatures(object = gcdata, mean.function = ExpMean, dispersion.function = LogVMR)
gcdata <- Seurat::ScaleData(object = gcdata, do.center = T, do.scale = F)
gcdata <- Seurat::RunPCA(object = gcdata, pc.genes = gcdata@var.genes, do.print = TRUE, pcs.print = 1:5,
                 genes.print = 5)
gcdata <- Seurat::FindNeighbors(gcdata, reduction="pca", dims = 1:5)

gcdata <- Seurat::FindClusters(object = gcdata,
                       resolution = 0.6)
```

Now we can add these clusters to the `slingshot` function:

```R
deng_SCE$Seurat_clusters <- as.character(Idents(gcdata))  # go from factor to character
deng_SCE <- slingshot::slingshot(deng_SCE, clusterLabels = 'Seurat_clusters', reducedDim = 'PCA')
```

There have been added two `slingPseudotime` columns:

```R
head(colData(deng_SCE))
```

Let's see whether things have improved:

```R
PCAplot_slingshot(deng_SCE, variable = deng_SCE$slingPseudotime_2)
```

```R
ggplot(as.data.frame(colData(deng_SCE)), aes(x = slingPseudotime_2, y = cell_type2,
                                             colour = cell_type2)) +
  ggbeeswarm::geom_quasirandom(groupOnX = FALSE) +
  theme_classic() +
  xlab("Slingshot pseudotime") + ylab("Timepoint") +
  ggtitle("Cells ordered by Slingshot pseudotime")
```

Particularly the later stages, separation seems to improve. Since we have included the Seurat clustering, we can plot the PCA, with colors according to these clusters:

```R
PCAplot_slingshot(deng_SCE, variable = deng_SCE$Seurat_clusters, type = 'lineages', col = 'black', legend = TRUE)
PCAplot_slingshot(deng_SCE, variable = deng_SCE$cell_type2, type = 'lineages', col = 'black', legend = TRUE)
```
