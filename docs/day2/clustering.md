## Material

[:fontawesome-solid-file-pdf: Download the presentation](../assets/pdf/scRNAseq_clustering_RM_november.pdf){: .md-button }

- Evaluation of [clustering methods](https://f1000research.com/articles/7-1141)

## Exercises

The method implemented in Seurat first constructs a SNN graph based on the euclidean distance in PCA space, and refine the edge weights between any two cells based on the shared overlap in their local neighborhoods (Jaccard similarity). This step is performed using the `FindNeighbors()` function, and takes as input the previously defined dimensionality of the dataset.

!!! note
    We use the integrated object (`seu_int`) and the assay `integrated`. Unsure? Check `DefaultAssay(seu_int)`, and set it by `DefaultAssay(seu_int) <- "integrated"`. 

```R
seu_int <- Seurat::FindNeighbors(seu_int, dims = 1:25)
```

To cluster the cells, Seurat next implements modularity optimization techniques such as the Louvain algorithm (default) or SLM [SLM, Blondel et al., Journal of Statistical Mechanics], to iteratively group cells together, with the goal of optimizing the standard modularity function. The `FindClusters()` function implements this procedure, and contains a resolution parameter that sets the ‘granularity’ of the downstream clustering, with increased values leading to a greater number of clusters.

```R
seu_int <- Seurat::FindClusters(seu_int, resolution = seq(0.1, 0.8, by=0.1))
```

Cluster id of each cell is added to the metadata object, as a new column for each resolution tested:

```R
head(seu_int@meta.data)
```

To view how clusters sub-divide at increasing resolution:

```R
library(clustree)
clustree::clustree(seu_int@meta.data[,grep("integrated_snn_res", colnames(seu_int@meta.data))],
                   prefix = "integrated_snn_res.")
```

You can view the UMAP coloring each cell according to a cluster id like this:

```R
Seurat::DimPlot(seu_int, group.by = "integrated_snn_res.0.1")
```

**Exercise:** Visualise clustering based on a few more resolutions. Taking the clustering and the UMAP plots into account what do you consider as a good resolution to perform the clustering?

??? done "Answer"
    Of course, there is no 'optimal' resolution, but based on resolution of 0.3, the tree stays relatively stable for a few resolution steps, and it seems that clustering fits the UMAP well:

    ```R
    Seurat::DimPlot(seu_int, group.by = "integrated_snn_res.0.3")
    ```

    <figure>
      <img src="../../assets/images/UMAP_res_0.3.png" width="400"/>
    </figure>

**Exercise:** When do the number of neighbors need to be changed? How does changing the method of clustering in `FindClusters` affect the output? Which parameter should be changed?

??? done "Answer"
    As FindClusters is an unsupervised clustering method supposedly telling yousomething about your UMAP plot, the two should go along. If one has reasons to change the number of neighbors in the UMAP function, here the same parameter should be adapted.
    
    The method can be changed with algorithm = 2,3 or 4

### Save the dataset and clear environment

Now, save the dataset so you can use it later today:

```R
saveRDS(seu_int, "seu_int_day2_part1.rds")
```

Clear your environment:

```R
rm(list = ls())
gc()
.rs.restartR()
```
