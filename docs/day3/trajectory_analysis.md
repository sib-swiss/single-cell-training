## Learning outcomes

**After having completed this chapter you will be able to:**

## Material

[:fontawesome-solid-file-pdf: Download the presentation](../assets/pdf/sequencing_technologies.pdf){: .md-button }

- `slingshot` [vignette](https://www.bioconductor.org/packages/release/bioc/vignettes/slingshot/inst/doc/vignette.html)

## Exercises

Load the following packages:

```R
library(SingleCellExperiment)
library(BiocSingular)
library(slingshot)
library(ggplot)
library(ggbeeswarm)
```

### Trajectory analysis using Slingshot

> :material-transit-connection-variant: This part uses the `Deng` dataset

Read in data. It is an object of class `SingleCellExperiment`.

```R
deng_SCE <- readRDS("data/deng_dataset/deng-reads.rds")
```

Perform the first steps of the analysis. The deng_SCE object contains cells that were isolated at different stages of mouse embryogenesis, from the zygote stage to the late blastula.

The levels of the cell type are in alphabetical order. We now change the level order for plotting in developmental order:

```R
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
```

Have a look at the PCA with the slingshot pseudotime line:

```R
PCAplot_slingshot(sce, variable = sce$slingPseudotime_1, draw_lines = TRUE)
```

Also have a look at pseudotime versus cell type:

```R
ggplot(as.data.frame(colData(deng_SCE)), aes(x = sce$slingPseudotime_1,
                                             y = cell_type2,
                                             colour = cell_type2)) +
  ggbeeswarm::geom_quasirandom(groupOnX = FALSE) +
  theme_classic() +
  xlab("Slingshot pseudotime") + ylab("Timepoint") +
  ggtitle("Cells ordered by Slingshot pseudotime")
```

This already looks pretty good. Let's see whether we can improve it. First we generate clusters by using `Seurat`:

```R
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
                                dims = 1:5)

# clustering with resolution of 0.6
gcdata <- Seurat::FindClusters(object = gcdata,
                               resolution = 0.6)
```

**Exercise:** Have a look at the UMAP and color it according to the clustering (we used a resolution of 0.6). Does it look acceptable to you?

!!! note
    The UMAP has not yet been generated for this object.

??? done "Answer"

    First run the UMAP and after that generate the plot:

    ```R
    gcdata <- Seurat::RunUMAP(gcdata, dims = 1:5)
    Seurat::DimPlot(gcdata)
    ```

    Returns:

    <figure>
     <img src="../../assets/images/umap_gcdata.png" width="400"/>
    </figure>


Now we can add these clusters to the `slingshot` function:

```R
deng_SCE$Seurat_clusters <- as.character(Idents(gcdata))  # go from factor to character

deng_SCE <- slingshot::slingshot(deng_SCE,
                                 clusterLabels = 'Seurat_clusters',
                                 reducedDim = 'PCA')
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
ggplot(as.data.frame(colData(deng_SCE)), aes(x = slingPseudotime_2,
                                             y = cell_type2,
                                             colour = cell_type2)) +
  ggbeeswarm::geom_quasirandom(groupOnX = FALSE) +
  theme_classic() +
  xlab("Slingshot pseudotime") + ylab("Timepoint") +
  ggtitle("Cells ordered by Slingshot pseudotime")
```

Particularly the later stages, separation seems to improve. Since we have included the Seurat clustering, we can plot the PCA, with colors according to these clusters:

```R
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
```
