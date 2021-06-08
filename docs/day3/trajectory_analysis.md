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
library(ggplot2)
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

Use the `reducedDim` function to access the PCA and store the results.

```R
pca <- SingleCellExperiment::reducedDim(deng_SCE, "PCA")
```

Describe how the PCA is stored in a matrix. Why does it have this structure?

```R
head(pca)
```

Add PCA data to the deng_SCE object.

```R
deng_SCE$PC1 <- pca[, 1]
deng_SCE$PC2 <- pca[, 2]
```

Plot PC biplot with cells colored by cell_type2.
`colData(deng_SCE)` accesses the cell metadata `DataFrame` object for `deng_SCE`.
Look at Figure 1A of the [paper](https://science.sciencemag.org/content/343/6167/193) as a comparison to your PC biplot.

```R
ggplot(as.data.frame(colData(deng_SCE)), aes(x = PC1, y = PC2, color = cell_type2)) +
  geom_point(size=2, shape=20) +
  theme_classic() +
  xlab("PC1") + ylab("PC2") + ggtitle("PC biplot")
```

PCA is a simple approach and can be good to compare to more complex algorithms
designed to capture differentiation processes. As a simple measure of pseudotime
we can use the coordinates of PC1.
Plot PC1 vs cell_type2.

```R
deng_SCE$pseudotime_PC1 <- rank(deng_SCE$PC1)  # rank cells by their PC1 score
```

Create a jitter plot

```
ggplot(as.data.frame(colData(deng_SCE)), aes(x = pseudotime_PC1, y = cell_type2,
                                             colour = cell_type2)) +
  ggbeeswarm::geom_quasirandom(groupOnX = FALSE) +
  theme_classic() +
  xlab("PC1") + ylab("Timepoint") +
  ggtitle("Cells ordered by first principal component")
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
    lines(slingshot::SlingshotDataSet(sce), lwd = 2, ... )
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
                                reduction = "pca",
                                dims = 1:5)

# clustering with resolution of 0.6
gcdata <- Seurat::FindClusters(object = gcdata,
                               resolution = 0.6)
```

<!-- **Exercise:** Have a look at the UMAP and color it according to the clustering (we used a resolution of 0.6). Does it look acceptable to you?

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
    </figure> -->


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

Check how the slingshot object has evolved

```R
SlingshotDataSet(deng_SCE)
```

Plot PC1 versus PC2 colored by slingshot pseudotime:

```R
PCAplot_slingshot(deng_SCE, variable = deng_SCE$slingPseudotime_2)
```

Plot Slingshot pseudotime vs cell stage.

```R
ggplot(data.frame(cell_type2 = deng_SCE$cell_type2,
                  slingPseudotime_1 = deng_SCE$slingPseudotime_1),
        aes(x = slingPseudotime_1, y = cell_type2,
        colour = cell_type2)) +
  ggbeeswarm::geom_quasirandom(groupOnX = FALSE) +
  theme_classic() +
  xlab("Slingshot pseudotime") + ylab("Timepoint") +
  ggtitle("Cells ordered by Slingshot pseudotime")

ggplot(data.frame(cell_type2 = deng_SCE$cell_type2,
                  slingPseudotime_2 = deng_SCE$slingPseudotime_2),
        aes(x = slingPseudotime_2, y = cell_type2,
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

Clear your environment:

```R
rm(list = ls())
gc()
.rs.restartR()
```

### Trajectory analysis with `monocle3`

> :fontawesome-solid-ribbon: This part uses the `gbm` dataset

This part showcases how you can use `monocle3` to perform a trajectory analysis. First load the `gbm` dataset:

```R
gbm <- readRDS("gbm_day3.rds")
```

Load the required package into your environment:

```R
library(monocle3)
```

Generate a `monocle3` object (with class `cell_data_set`) from our `Seurat` object:

```R
feature_names <- as.data.frame(rownames(gbm))
rownames(feature_names) <- rownames(gbm)
colnames(feature_names) <- "gene_short_name"
gbm_monocl <- monocle3::new_cell_data_set(gbm@assays$RNA@counts,
                              cell_metadata = gbm@meta.data,
                              gene_metadata = feature_names)
```

We pre-process the newly created object. What does it involve? Check:

```R
?preprocess_cds
```

Preprocess the dataset:

```R
gbm_monocl <- monocle3::preprocess_cds(gbm_monocl)
```

And check out the elbow plot:

```R
monocle3::plot_pc_variance_explained(gbm_monocl)
```

Perform UMAP using the implementation in the `monocle3` package and its default parameters:

```R
gbm_monocl <- monocle3::reduce_dimension(gbm_monocl, reduction_method = "UMAP")
```

Plot the `monocle3` UMAP coloring cells according to the cluster ID ran with `Seurat`:

```R
monocle3::plot_cells(gbm_monocl, color_cells_by = "RNA_snn_res.0.2")
monocle3::plot_cells(gbm_monocl, genes = "PMP2") # to plot expression level of a gene
```

Cluster cells using `monocle3`'s clustering function:

```R
gbm_monocl <- monocle3::cluster_cells(gbm_monocl, resolution=0.00025)
p1 <- monocle3::plot_cells(gbm_monocl, label_cell_groups = F)
p2 <- monocle3::plot_cells(gbm_monocl, color_cells_by = "RNA_snn_res.0.2", label_cell_groups = F)
cowplot::plot_grid(p1, p2, ncol = 2) # Are there differences?
```

learn graph (i.e. identify trajectory) using `monocle3` UMAP and clustering:

```R
gbm_monocl <- monocle3::learn_graph(gbm_monocl)
monocle3::plot_cells(gbm_monocl)
monocle3::plot_cells(gbm_monocl, color_cells_by = "RNA_snn_res.0.2")
```

Replace `monocle3` cluster id by `Seurat` cluster id if we want to keep the same information:

```R
gbm_monocl@clusters$UMAP$clusters <- colData(gbm_monocl)$RNA_snn_res.0.2
names(gbm_monocl@clusters$UMAP$clusters) <- rownames(colData(gbm_monocl))
gbm_monocl <- monocle3::learn_graph(gbm_monocl)
monocle3::plot_cells(gbm_monocl, label_cell_groups = F)
```

Select the "initial" cells to calculate pseudotime. A pop up window will open and you need to click on the "initial" cells (one node per trajectory), then click "Done".

```R
gbm_monocl<-monocle3::order_cells(gbm_monocl)#
```

```R
monocle3::plot_cells(gbm_monocl,
           color_cells_by = "pseudotime",
           label_cell_groups=F,
           label_leaves=F,
           label_branch_points=FALSE,
           graph_label_size=1.5, cell_size = 1)
```

Plot a gene's expression vs pseudotime

```R
plot_genes_in_pseudotime(subset(gbm_monocl, rowData(gbm_monocl)$gene_short_name=="PMP2"))
```
