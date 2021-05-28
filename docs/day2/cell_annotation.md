
## Learning outcomes

**After having completed this chapter you will be able to:**

## Material

[:fontawesome-solid-file-pdf: Download the presentation](../assets/pdf/sequencing_technologies.pdf){: .md-button }

## Exercises

> :material-zodiac-cancer: This chapter uses the `gbm` dataset

Load the `gbm` dataset you have created earlier today:

```R
gbm <- readRDS("gbm_day2_part1.rds")
```

In the last exercise we saw that probably clustering at a resolution of 0.2 gave the most sensible results. Let's therefore set the default identity of each cell based on this clustering:

```R
gbm <- Seurat::SetIdent(gbm, value = gbm$RNA_snn_res.0.2)
```

!!! note
    From now on, grouping (e.g. for plotting) is done by the active identity (set at `@active.ident`) by default.

Based on the UMAP we have generated, we can visualize expression for a gene in each cluster:

```R
Seurat::FeaturePlot(gbm, "PMP2")
```

Based on expression of sets of genes you can do a manual cell type annotation. If you know the marker genes for some cell types, you can check whether they are up-regulated in one or the other cluster. Here we have some marker genes for two different cell types:

```R
immune_genes<-c("GZMA", "CD3E", "CD3D")
microglia_genes<-c("CCL4", "CCL3", "P2RY12", "C1QB", "CSF1ER", "CY3CR1")
```

Let's have a look at the expression of the three immune genes:

```R
Seurat::FeaturePlot(gbm, immune_genes, ncol=2)
```

These cells are almost all in cluster 6. Which becomes clearer when looking at the violin plot:

```R
Seurat::VlnPlot(gbm,
                features = immune_genes
                ncol = 2)
```

**Exercise:** Have a look at the microglia genes as well. Which cluster contains probably microglial cells?

??? done "Answer"
    Running

    ```R
    Seurat::FeaturePlot(gbm, microglia_genes, ncol=2)
    ```

    Returns:

    <figure>
      <img src="../../assets/images/featureplots_microglia.png" width="400"/>
    </figure>

    Corresponding mainly to cluster 1 (and cluster 4):

    ```R
    Seurat::VlnPlot(gbm,
                    features = microglia_genes,
                    ncol = 2)
    ```

    <figure>
      <img src="../../assets/images/violinplots_microglia.png" width="300"/>
    </figure>

We can also automate this with the function `AddModuleScore`. For each cell, an expression score for a group of genes is calcuated:

```R
gbm <- Seurat::AddModuleScore(gbm,
                              features = list(immune_genes),
                              name = "immune_genes")
```

**Exercise:** After running `AddModuleScore`, a column was added to `gbm@meta.data`.

**A.** What is the name of that column? What kind of data is in there?

**B.** Generate a UMAP with color accoding to this column and a violinplot grouped by cluster. Is this according to what we saw in the previous exercise?

??? done "Answer"

    **A.** The new column is called `immune_genes1`. It contains the module score for each cell (which is basically the average expression of the set of genes).

    **B.** You can plot the UMAP with

    ```R
    Seurat::FeaturePlot(gbm, "immune_genes1")
    ```

    Returning:

    <figure>
      <img src="../../assets/images/UMAP_immune_genes_modules.png" width="300"/>
    </figure>

    ```R
    Seurat::VlnPlot(gbm,
                    "immune_genes1")
    ```

    Which indeed shows these genes are mainly expressed in cluster 6:

    <figure>
      <img src="../../assets/images/violinplot_immune_genes_modules.png" width="300"/>
    </figure>

### Cell type annotation using `SingleR`

To do a fully automated annoation, we need a reference dataset of primary cells. Any reference could be used. The package `scRNAseq` in Bioconductor includes several scRNAseq datasets that can be used as reference to `SingleR`. One could also use a reference made of bulk RNA seq data. Here we are using the Human Primary Cell Atlas dataset from `celldex`. Check out what's in there:

```R
hpca.se <- celldex::HumanPrimaryCellAtlasData()
class(hpca.se)
table(hpca.se$label.main)
```

Now `SingleR` compares our normalized count data to a reference set, and finds the most probable annation:

```R
gbm_SingleR <- SingleR::SingleR(test = Seurat::GetAssayData(gbm, slot = "data"),
                                ref = hpca.se,
                                labels = hpca.se$label.main)
```

See what's in there by using `head`:

```R
head(gbm_SingleR)
```

In order to visualize it in our UMAP, we have to add the annotation to `gbm@meta.data`:

```R
gbm$SingleR_annot <- gbm_SingleR$labels

Seurat::DimPlot(gbm, group.by = "SingleR_annot", label = T, repel = T)
```

**Exercise:** Compare our manual annotation (based on the set of immune genes) to the annotation with `SingleR`. Do they correspond?

??? done "Answer"
    We can have a look at the mean module score for each `SingleR` annotation like this:

    ```R
    mean_scores <- tapply(gbm$immune_genes1, gbm$SingleR_annot, mean)
    mean_scores[order(mean_scores, decreasing = TRUE)[1:6]]
    ```

    Returning:

    ```
        T_cells             NK_cell                  DC            Monocyte          Macrophage        Chondrocytes
    1.015605822         0.679592135         0.017018871        -0.005414908        -0.011926460        -0.011975024
    ```

    Showing that T-cells and NK-cells have a high module score based on our set of immune genes, which makes sense.

    Of course, it was also already clear from the UMAP plots that cluster 6 (the cluster with the high module score for the immune genes) contained the T-cells and NK-cells.

### Save the dataset and clear environment

Now, save the dataset so you can use it tomorrow:

```R
saveRDS(gbm, "gbm_day2_part2.rds")
```

Clear your environment:

```R
rm(list = ls())
gc()
.rs.restartR()
```
