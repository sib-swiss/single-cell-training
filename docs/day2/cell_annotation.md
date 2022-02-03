
## Material

[:fontawesome-solid-file-pdf: Download the presentation](../assets/pdf/cell_annotation_Day2_scRNAseq_112021.pdf){: .md-button }

- [Review on automated cell annotation](https://www.sciencedirect.com/science/article/pii/S2001037021000192)

## Exercises

Load the `seu_int` dataset you have created earlier today:

```R
seu_int <- readRDS("seu_int_day2_part1.rds")
```

And load the following packages:

```R
library(celldex)
library(SingleR)
```

In the last exercise we saw that probably clustering at a resolution of 0.3 gave the most sensible results. Let's therefore set the default identity of each cell based on this clustering:

```R
seu_int <- Seurat::SetIdent(seu_int, value = seu_int$integrated_snn_res.0.3)
```

!!! note
    From now on, grouping (e.g. for plotting) is done by the active identity (set at `@active.ident`) by default.

During cell annotation we will use the original count data (not the integrated data):

```R
DefaultAssay(seu_int) <- "RNA"
```


Based on the UMAP we have generated, we can visualize expression for a gene in each cluster:

```R
Seurat::FeaturePlot(seu_int, "HBA1")
```

Based on expression of sets of genes you can do a manual cell type annotation. If you know the marker genes for some cell types, you can check whether they are up-regulated in one or the other cluster. Here we have some marker genes for two different cell types:

```R
tcell_genes <- c("IL7R", "LTB", "TRAC", "CD3D")
monocyte_genes <- c("CD14", "CST3", "CD68", "CTSS")
```

Let's have a look at the expression of the four T cell genes:

```R
Seurat::FeaturePlot(seu_int, tcell_genes, ncol=2)
```

These cells are almost all in cluster 0 and 8. Which becomes clearer when looking at the violin plot:

```R
Seurat::VlnPlot(seu_int,
                features = tcell_genes,
                ncol = 2)
```

**Exercise:** Have a look at the monocyte genes as well. Which clusters contain probably monocytes?

??? done "Answer"
    Running

    ```R
    Seurat::FeaturePlot(seu_int, monocyte_genes, ncol=2)
    ```

    Returns:

    <figure>
      <img src="../../assets/images/featureplots_monocytes.png" width="500"/>
    </figure>

    Corresponding mainly to cluster 2 and 9:

    ```R
    Seurat::VlnPlot(seu_int,
                    features = monocyte_genes,
                    ncol = 2)
    ```

    <figure>
      <img src="../../assets/images/violinplots_monocytes.png" width="500"/>
    </figure>

We can also automate this with the function `AddModuleScore`. For each cell, an expression score for a group of genes is calcuated:

```R
seu_int <- Seurat::AddModuleScore(seu_int,
                              features = list(tcell_genes),
                              name = "tcell_genes")
```

**Exercise:** After running `AddModuleScore`, a column was added to `seu_int@meta.data`.

**A.** What is the name of that column? What kind of data is in there?

**B.** Generate a UMAP with color accoding to this column and a violinplot grouped by cluster. Is this according to what we saw in the previous exercise?

??? done "Answer"

    **A.** The new column is called `tcell_genes1`. It contains the module score for each cell (which is basically the average expression of the set of genes).

    **B.** You can plot the UMAP with

    ```R
    Seurat::FeaturePlot(seu_int, "tcell_genes1")
    ```

    Returning:

    <figure>
      <img src="../../assets/images/UMAP_immune_genes_modules.png" width="500"/>
    </figure>

    ```R
    Seurat::VlnPlot(seu_int,
                    "tcell_genes1")
    ```

    Which indeed shows these genes are mainly expressed in cluster 6:

    <figure>
      <img src="../../assets/images/violinplot_immune_genes_modules.png" width="500"/>
    </figure>

### Annotating cells according to cycling phase

```R
s.genes <- Seurat::cc.genes.updated.2019$s.genes
g2m.genes <- Seurat::cc.genes.updated.2019$g2m.genes
```

```R
seu_int <- Seurat::CellCycleScoring(seu_int,
                                     s.features = s.genes,
                                     g2m.features = g2m.genes)
```

```R
Seurat::DimPlot(seu_int, group.by = "Phase")
```


### Cell type annotation using `SingleR`

To do a fully automated annoation, we need a reference dataset of primary cells. Any reference could be used. The package `scRNAseq` in Bioconductor includes several scRNAseq datasets that can be used as reference to `SingleR`. One could also use a reference made of bulk RNA seq data. Here we are using the a hematopoietic reference dataset from `celldex`. Check out what's in there:

```R
ref <- celldex::NovershternHematopoieticData()
class(ref)
table(ref$label.main)
```

!!! note
    You can find more information on different reference datasets at the [`celldex` documentation](https://bioconductor.org/packages/3.14/data/experiment/vignettes/celldex/inst/doc/userguide.html)

Now `SingleR` compares our normalized count data to a reference set, and finds the most probable annation:

```R
seu_int_SingleR <- SingleR::SingleR(test = Seurat::GetAssayData(seu_int, slot = "data"),
                                ref = ref,
                                labels = ref$label.main)
```

See what's in there by using `head`:

```R
head(seu_int_SingleR)
```

Visualize singleR score quality scores:

```R
SingleR::plotScoreHeatmap(seu_int_SingleR)
SingleR::plotDeltaDistribution(seu_int_SingleR)
```

There are some annotations that contain only a few cells. They are usually not of interest, and they clogg our plots. Therefore we remove them from the annotation:

```R
singleR_labels <- seu_int_SingleR$labels
t <- table(singleR_labels)
other <- names(t)[t < 10]
singleR_labels[singleR_labels %in% other] <- NA
```

In order to visualize it in our UMAP, we have to add the annotation to `seu_int@meta.data`:

```R
seu_int$SingleR_annot <- singleR_labels
```

We can plot the annotations in the UMAP. Here, we use a different package for plotting (`dittoSeq`) as it has a bit better default coloring, and some other plotting functionality we will use later on.

```R
dittoSeq::dittoDimPlot(seu_int, "SingleR_annot", size = 0.7)
```

We can check out how many cells per sample we have for each annotated cell type:

```R
dittoSeq::dittoBarPlot(seu_int, var = "SingleR_annot", group.by = "orig.ident")
```

**Exercise:** Compare our manual annotation (based on the set of T cell genes) to the annotation with `SingleR`. Do they correspond?

!!! hint
    You can for example use the plotting function `dittoBarPlot` to visualize the cell types according to cluster (use `integrated_snn_res.0.3` in stead of `orig.ident`))

??? done "Answer"
    We can have a look at the mean module score for each `SingleR` annotation like this:

    ```R
    dittoSeq::dittoBarPlot(seu_int, 
                       var = "SingleR_annot", 
                       group.by = "integrated_snn_res.0.3")
    ```

    This returns:

    <figure>
      <img src="../../assets/images/barplot_cluster_singler.png" width="700"/>
    </figure>

    Here, you can see that cluster 0 and 8 contain cells annotated as T cells (CD4+ and CD8+).
    



### Save the dataset and clear environment

Now, save the dataset so you can use it tomorrow:

```R
saveRDS(seu_int, "seu_int_day2_part2.rds")
```

Clear your environment:

```R
rm(list = ls())
gc()
.rs.restartR()
```
