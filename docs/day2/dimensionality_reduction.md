
## Material

[:fontawesome-solid-file-pdf: Download the presentation](../assets/pdf/scRNAseq_RM_dimreduction_only.pdf){: .md-button }

- Making sense of [PCA](https://stats.stackexchange.com/questions/2691/making-sense-of-principal-component-analysis-eigenvectors-eigenvalues)
- Understanding [t-SNE](https://distill.pub/2016/misread-tsne/)
- [t-SNE explained](https://www.youtube.com/watch?v=NEaUSP4YerM) by Josh Starmer
- Understanding [UMAP](https://pair-code.github.io/understanding-umap/)
- [Video](https://www.youtube.com/watch?v=nq6iPZVUxZU) by one of the UMAP authors
- More info on [UMAP parameters](https://umap.scikit-tda.org/parameters.html)


## Exercises

> :fontawesome-solid-ribbon: This chapter uses the `gbm` dataset

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

**Exercise:** Try to change:

**A.** The number of neighbors used for the calculation of the UMAP. Which is the parameter to change and how did it affect the output. What is the default ? In which situation would you lower/increase this ?

**B.** The number of dims to extremes dims = 1:5 or dims = 1:50 how
did it affect the output ? In your opinion better few PCAs too much or too few ?
Why does dims = 1:100 not work ? When would more precision be needed ?

??? done "Answer"
    **Answer A**

    ```R 
    gbm <- Seurat::RunUMAP(gbm, dims = 1:25,n.neighbors = 50)
    ```

    It can be of interest to change the number of neighbors if one has subset the data (for instance in the situation where you would only consider the t-cells inyour data set), then maybe the number of neighbors in a cluster would anyway be most of the time lower than 30 then 30 is too much. In the other extreme where your dataset is extremely big an increase in the number of neighbors can be considered.

    **Answer B** 

    ```R
    gbm <- Seurat::RunUMAP(gbm, dims = 1:5)
    ```

    ```R
    gbm <- Seurat::RunUMAP(gbm, dims = 1:50) 
    ```
    Taking dims = 1:100 does not work as in the step RunPCA by default only 50pcs are calculated, so the maximum that we can consider in further steps are 50, if more precision makes sense, for instance, if the genes that is of interest for your study is not present when the RunPCA was calculated, then an increase in the number of components calculated at start might be interesting tobe considered. Taking too few PCs we have a « blob » everything looks connected. Too many PCs tends to separate everything. Personally it is more interesting for me too have maybe 2 clusters separated of epithelial cells that I then group for further downstream analysis rather than having very distinct cells being clustered together. So I would rather take the « elbow » of the elbow plota bit further to the right.
