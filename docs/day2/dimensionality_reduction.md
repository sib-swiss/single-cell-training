## Learning outcomes

**After having completed this chapter you will be able to:**

## Material

[:fontawesome-solid-file-pdf: Download the presentation](../assets/pdf/sequencing_technologies.pdf){: .md-button }

## Exercises

> :fontawesome-solid-ribbon: This chapter uses the `gbm` dataset

### Dimensionality reduction using Seurat

Load the `gbm` dataset you have created yesterday:

```R
gbm <- readRDS("gbm_day1.rds")
```

And load the following packages:

```R
library(Seurat)
library(clustree)
```

Once the data is normalized, scaled and variable features have been identified, we can start to reduce the dimensionality of the data.
For the PCA, by default, only the previously determined variable features are used as input, but can be defined using features argument if you wish to specify a vector of genes. The PCA will only be run on the variable features, that you can check with `VariableFeatures(gbm)`.

```R
gbm <- Seurat::RunPCA(gbm)
```

To view the PCA plot:

```R
Seurat::DimPlot(gbm, reduction = "pca")
```

We can colour the PCA plot according to any factor that is present in `@meta.data`. For example we can take the column `Phase` (i.e. predicted cell cycle phase):

```R
Seurat::DimPlot(gbm, reduction = "pca", group.by = "Phase")
```

!!! note
    Coming back to the cell cycle analysis, we can check the distribution of the different cell cycle phases over the PCA, and eventually regress it out using the `ScaleData()` function. But here, the PCA doesn't seem to cluster according to the cell cycle phase.

We can generate heatmaps according to the correlations with the different dimensions of our PCA:

```R
Seurat::DimHeatmap(gbm, dims = 1:12, cells = 500, balanced = TRUE)
```

The elblowplot can help you in determining how many PCs to use for downstream analysis such as UMAP:

```R
Seurat::ElbowPlot(gbm, ndims = 40)
```

The elblow plot ranks principle components based on the percentage of variance explained by each one. Where we observe an "elblow" or flattening curve, the majority of true signal is captured by this number of PCs, eg around 25 PCs for the gbm dataset.

Including too many PCs usually does not affect much the result, while including too few PCs can affect the results very much.

UMAP: The goal of these algorithms is to learn the underlying manifold of the data in order to place similar cells together in low-dimensional space.

```R
gbm <- Seurat::RunUMAP(gbm, dims = 1:25)
```

To view the UMAP plot:

```R
Seurat::DimPlot(gbm, reduction = "umap")
```

Cells can be coloured according to cell cycle phase.
Is there a group of cells than contains a high proportion of cells in G2/M phase?

```R
Seurat::DimPlot(gbm, reduction = "umap", group.by = "Phase")
```

### Clustering

The method implemented in Seurat first constructs a KNN graph based on the euclidean distance in PCA space, and refine the edge weights between any two cells based on the shared overlap in their local neighborhoods (Jaccard similarity). This step is performed using the `FindNeighbors()` function, and takes as input the previously defined dimensionality of the dataset.

```R
gbm <- Seurat::FindNeighbors(gbm, dims = 1:25)
```

To cluster the cells, Seurat next implements modularity optimization techniques such as the Louvain algorithm (default) or SLM [SLM, Blondel et al., Journal of Statistical Mechanics], to iteratively group cells together, with the goal of optimizing the standard modularity function. The `FindClusters()` function implements this procedure, and contains a resolution parameter that sets the ‘granularity’ of the downstream clustering, with increased values leading to a greater number of clusters.

```R
gbm <- Seurat::FindClusters(gbm, resolution = seq(0.1, 0.8, by=0.1))
```

Cluster id of each cell is added to the metadata object, as a new column for each resolution tested:

```R
head(gbm@meta.data)
```

To view how clusters sub-divide at increasing resolution:

```R
library(clustree)
clustree::clustree(gbm@meta.data[,grep("RNA_snn_res", colnames(gbm@meta.data))],
                   prefix = "RNA_snn_res.")
```

You can view the UMAP coloring each cell according to a cluster id like this:

```R
Seurat::DimPlot(gbm, group.by = "RNA_snn_res.0.1")
```

**Exercise:** Visualise clustering based on a few more resolutions. Taking the clustering and the UMAP plots into account what do you consider as a good resolution to perform the clustering?

??? done "Answer"
    Of course, there is no 'optimal' resolution, but based on resolution of 0.2, it seems that clustering fits the UMAP well:

    ```R
    Seurat::DimPlot(gbm, group.by = "RNA_snn_res.0.2")
    ```

    <figure>
      <img src="../../assets/images/UMAP_res.0.2.png" width="400"/>
    </figure>

### Save the dataset and clear environment

Now, save the dataset so you can use it later today:

```R
saveRDS(gbm, "gbm_day2_part1.rds")
```

Clear your environment:

```R
rm(list = ls())
gc()
.rs.restartR()
```
