## Learning outcomes

**After having completed this chapter you will be able to:**

## Material

[:fontawesome-solid-file-pdf: Download the presentation](../assets/pdf/sequencing_technologies.pdf){: .md-button }

## Exercises

The `gbm` dataset does not contain any samples, treatments or methods to integrate. Therefore for these exercises we will use a different dataset that is described in [Comprehensive Integration of Single CellData](https://www.biorxiv.org/content/10.1101/460147v1). It is a dataset comprising of four different single cell experiment performed by using four different methods.

```R
pancreas.data <- readRDS(file = "data/pancreas_v3_files/pancreas_expression_matrix.rds")
metadata <- readRDS(file = "data/pancreas_v3_files/pancreas_metadata.rds")
```

Create a Seurat object with all datasets.

```R
pancreas <- Seurat::CreateSeuratObject(pancreas.data, meta.data = metadata)
```

The object `pancreas` is now of class `Seurat` and comparable with the object `gbm` that we have used in the previous exercises.

**Exercise:** Have a look at the object. How many cells are in there? And how many features? What kind of information is in the `meta.data` slot?

??? done "Answer"
    Just by running `pancreas`, we get the following information:

    ```
    An object of class Seurat
    34363 features across 6321 samples within 1 assay
    Active assay: RNA (34363 features, 0 variable features)
    ```

    Running `head(pancreas@meta.data)` gives you:

    ```
                  orig.ident nCount_RNA nFeature_RNA   tech celltype
    D101_5     SeuratProject   4615.810         1986 celseq    gamma
    D101_43    SeuratProject  11711.506         3942 celseq    gamma
    D101_93    SeuratProject   5567.659         2418 celseq    gamma
    D102_4     SeuratProject   6804.533         2846 celseq    gamma
    D172444_23 SeuratProject   5541.101         2436 celseq    gamma
    D172444_68 SeuratProject   4301.892         2015 celseq    gamma
    ```

    So, apparently there is also information in there about the technology (`tech`) and cell type annotation (`celltype`)

Now, we can repeat what we have learned in the previous chapters. Let's assume we don't have to filter for e.g. mitochondrial UMIs and number of features, and we can directly proceed to the normalization and scaling.

**Exercise:** To perform normalization, scaling, PCA and UMAP, run the following functions (with sensible parameters) on the `pancreas` object:

- `Seurat::NormalizeData`
- `Seurat::FindVariableFeatures`
- `Seurat::ScaleData`
- `Seurat::RunPCA`
- `Seurat::RunUMAP`

??? done "Answer"

    ```R
    pancreas <- Seurat::NormalizeData(pancreas)
    pancreas <- Seurat::FindVariableFeatures(pancreas, selection.method = "vst", nfeatures = 2000)
    pancreas <- Seurat::ScaleData(pancreas)
    pancreas <- Seurat::RunPCA(pancreas, npcs = 30)
    pancreas <- Seurat::RunUMAP(pancreas, reduction = "pca", dims = 1:30)
    ```

Now we plot the UMAP based on technology. There is clearly clustering according to technology:

```R
Seurat::DimPlot(pancreas, reduction = "umap", group.by = "tech")
```

**Exercise:** Generate the same UMAP plot but now grouped by `celltype`. Does cell type group correctly together?

??? done "Answer"
    Generate the plot:
    ```R
    Seurat::DimPlot(pancreas, reduction = "umap", group.by = "celltype")
    ```
    This returns:

    <figure>
      <img src="../../assets/images/umap_celltype.png" width="400"/>
    </figure>

    This shows that within a techology cell types cluster together, but not between technology (this is e.g. very clear if you look at the clusters annotated as "alpha")

To perform the integration, we split the combined object into a list, with each dataset as an element. We perform standard preprocessing (log-normalization), and identify variable features individually for each dataset based on a variance stabilizing transformation (`"vst"`).

```R
pancreas.list <- Seurat::SplitObject(pancreas, split.by = "tech")

for (i in 1:length(pancreas.list)) {
    pancreas.list[[i]] <- Seurat::NormalizeData(pancreas.list[[i]], verbose = FALSE)
    pancreas.list[[i]] <- Seurat::FindVariableFeatures(pancreas.list[[i]], selection.method = "vst", nfeatures = 2000,
        verbose = FALSE)
}
```

After this, we prepare the integration by selecting integration anchors:

```R
pancreas.anchors <- Seurat::FindIntegrationAnchors(object.list = pancreas.list, dims = 1:30)
```

And finally perform the integration:

```R
pancreas.integrated <- Seurat::IntegrateData(anchorset = pancreas.anchors, dims = 1:30)
```

After running `IntegrateData`, the `Seurat` object will contain an additional element of class `Assay` with the integrated (or ‘batch-corrected’) expression matrix. This new `Assay` is called `integrated`, and exists next to the already existing `RNA` element with class `Assay`.

!!! warning
    Use the `Assay` `integrated` **only** for clustering and visualisation. It will give unexpected results during e.g. differential gene expression analysis. Therefore, use the `RNA` element for other analyses.

We can then use this new integrated matrix for clustering and visualization. Now, we can scale the integrated data, run PCA, and visualize the results with UMAP.

!!! note
    No need to re-run `FindVariableFeatures`, these were automatically set by calling `IntegrateData`.

First, switch the default `Assay` to `integrated` (in stead of `RNA`).

```R
Seurat::DefaultAssay(pancreas.integrated) <- "integrated"
```

**Exercise:** In order to redo the clustering, scale the integrated data, run the PCA and the UMAP again (using the function `ScaleData`, `RunPCA` and `RunUMAP`). After that, generate the same two UMAP plots (grouped by `tech` and by `celltype`). Did the integration perform well?

??? done "Answer"

    Performing the scaling, PCA and UMAP:

    ```R
    pancreas.integrated <- Seurat::ScaleData(pancreas.integrated, verbose = FALSE)
    pancreas.integrated <- Seurat::RunPCA(pancreas.integrated, npcs = 30, verbose = FALSE)
    pancreas.integrated <- Seurat::RunUMAP(pancreas.integrated, reduction = "pca", dims = 1:30)
    ```

    Plotting the UMAP:

    ```R
    Seurat::DimPlot(pancreas.integrated, reduction = "umap", group.by = "tech")
    Seurat::DimPlot(pancreas.integrated, reduction = "umap", group.by = "celltype", label = TRUE, repel = TRUE)
    ```
    Returning:

    <figure>
      <img src="../../assets/images/umap_integrated_tech.png" width="400"/>
    </figure>

    <figure>
      <img src="../../assets/images/umap_integrated_celltype.png" width="400"/>
    </figure>

    So, yes, integration performed well. Clustering is now not according to technology, but according to cell type.  

Finally, store the integrated dataset as an `.rds` file. We will use it tomorrow:

```R
saveRDS(pancreas.integrated, "pancreas.integrated.rds")
```
