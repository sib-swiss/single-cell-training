
## Material

[:fontawesome-solid-file-pdf: Download the presentation](../assets/pdf/scRNAseq_Day3_112021.pdf){: .md-button }

- More information on [pseudobulk analysis](https://bioconductor.org/books/release/OSCA/multi-sample-comparisons.html#differential-expression-between-condition)
- [Muscat](https://bioconductor.org/packages/devel/bioc/vignettes/muscat/inst/doc/analysis.html) for pseudobulk DGE.
- [Paper](https://www.nature.com/articles/nmeth.4612) on the robustness of different differential expression analysis methods

## Exercises

### Find all markers for each cluster

> :fontawesome-solid-ribbon: This part uses the `gbm` dataset

Load the `gbm` dataset you have created yesterday:

```R
gbm <- readRDS("gbm_day2_part2.rds")
```

And load the following packages:

```R
library(Seurat)
library(edgeR)
library(limma)
```

The function `FindAllMarkers` performs a Wilcoxon plot to determine the genes differentially expressed between each cluster and the rest of the cells. Other types of tests than the Wilcoxon test are available. Check it out by running `?Seurat::FindAllMarkers`.

Now run analysis:

```R
de_genes <- Seurat::FindAllMarkers(gbm,  min.pct = 0.25)
```

!!! note "Time for coffee"
    This takes a while. Have a break.

We are usually only interested in significant marker genes, so let's filter based on an adjusted p-value smaller than 0.05:

```R
de_genes <- subset(de_genes, de_genes$p_val_adj < 0.05)
View(de_genes)
```

**Exercise:** What are significant marker genes in cluster 6? Are the immune genes in there?

!!! hint
    You can re-load the vector with immune genes with:

    ```R
    immune_genes <- c("GZMA", "CD3E", "CD3D")
    ```

??? done "Answer"
    Running

    ```R
    de_genes[de_genes$gene %in% immune_genes,]
    ```

    Returns:
    ```
                  p_val avg_log2FC pct.1 pct.2     p_val_adj cluster gene
    GZMA 4.015094e-221   2.441962 0.331 0.005 9.781973e-217       6 GZMA
    CD3E  0.000000e+00   2.113342 0.474 0.006  0.000000e+00       6 CD3E
    CD3D  0.000000e+00   2.507228 0.519 0.006  0.000000e+00       6 CD3D
    ```

    So, yes, the immune genes are highly significant markers for cluster 6.

### Differential expression between clusters

> :fontawesome-solid-ribbon: This part uses the `gbm` dataset

The `FindMarkers` function allows to test for differential gene expression analysis specifically between 2 clusters, i.e. perform pairwise comparisons, eg between cells of cluster 0 vs cluster 2, or between cells annotated as astrocytes and macrophages.

First we can set the default cell identity to the cell types defined by `SingleR`:

```R
gbm <- Seurat::SetIdent(gbm, value = "SingleR_annot")
```

Run the differential gene expression analysis (runs for a couple of minutes):

```R
DEG_astro_vs_macro <- Seurat::FindMarkers(gbm,
                                           ident.1 = "Astrocyte",
                                           ident.2 = "Macrophage",
                                           group.by = gbm$SingleR_annot,
                                           test.use = "wilcox")
```

**Exercise:** What is the top 10 differentially expressed genes? What does the sign (i.e. positive or negative) mean in the log fold change values? Check your answer by generating a violin plot of a top differentially expressed gene.

??? done "Answer"
    You can look at the top 10 differentially expressed genes with:

    ```R
    top_order <- order(DEG_astro_vs_macro$p_val_adj)
    DEG_astro_vs_macro[top_order[1:10],]
    ```

    Returning:

    ```
              p_val avg_log2FC pct.1 pct.2 p_val_adj
    SLC2A5       0  -1.481971 0.080 0.656         0
    TNFRSF1B     0  -1.529484 0.041 0.694         0
    CAMK2N1      0   2.622397 0.907 0.106         0
    C1QA         0  -4.843247 0.088 0.964         0
    C1QC         0  -4.374756 0.075 0.943         0
    C1QB         0  -4.940663 0.099 0.946         0
    LAPTM5       0  -3.685947 0.075 0.987         0
    MARCKSL1     0   2.690161 0.911 0.221         0
    CNN3         0   2.940659 0.818 0.110         0
    CD53         0  -1.861605 0.035 0.804         0
    ```

    For an explanation of the log fold change have a look at `?Seurat::FindMarkers`. At **Value** it says:

    > `avg_logFC`: log fold-chage of the average expression between the two groups. Positive values indicate that the gene is more highly expressed in the first group

    Plotting the top gene `SLC2A5`:

    ```R
    Seurat::VlnPlot(gbm, features = "SLC2A5")
    ```

    Returning:

    <figure>
      <img src="../../assets/images/violinplot_SLC2A5.png" width="400"/>
    </figure>

### Differential expression analysis including batch as covariates

> :fontawesome-solid-disease: This part uses the `pancreas` dataset

The Wilcoxon test implemented in `FindMarkers` does not allow to test for complex design (eg factorial experiments) or to include batch as a covariate.

We can use `edgeR` or `limma` which are designed for microarray or bulk RNA seq data and provide a design matrix that includes covariates for example.

We will go back to the pancreas cells sequenced with different technologies, analyze differentially expressed genes between 2 clusters of cells using the technologies as covariates. Similar approaches can be used to analyze differentially expressed genes between conditions, eg sick vs healthy, wild type versus knockout, etc, and including batches in the model if they are present.

We will load the `pancreas.integrated` object we have created yesterday:

```R
pancreas.integrated <- readRDS("pancreas.integrated.rds")
```

Since we will start wit differential gene expression, we set the default assay back to "RNA". Also, we set the default identity to the cell type:

```R
Seurat::DefaultAssay(pancreas.integrated) <- "RNA"
Seurat::Idents(pancreas.integrated) <- pancreas.integrated$celltype
```

Let's have a look at the UMAP (again), coloured by celltype:

```R
Seurat::DimPlot(pancreas.integrated)
```

Let's say we are specifically interested to test for differential gene expression between two cell types.

!!! note
    Here we could also test for e.g. healthy versus diseased within a celltype/cluster.

Now we will run differential expression analysis between cell type *delta* and *gamma* using the technology as a covariate by using `limma`.

First, we will subset the `pancreas.integrated` object, only leaving the *delta* and *gamma* cells:

```R
pancreas.dg <- subset(pancreas.integrated, idents = c("delta", "gamma"))
```

Get the count matrix and keep only genes that are expressed in at least one cell:

```R
counts <- Seurat::GetAssayData(pancreas.dg, slot = "counts")
counts <- counts[rowSums(counts) != 0,]
```

Generate a `DGEList` object to use as input for `limma`:

```R
dge <- edgeR::DGEList(counts = counts)
dge <- edgeR::calcNormFactors(dge)  
```

Generate a design matrix:

```R
design <- model.matrix(~ 0 + celltype + tech, data = pancreas.dg@meta.data)
colnames(design)<-c("delta", "gamma", "celseq2", "fluidigmc1", "smartseq2")
```

Specify which contrasts to check:

```R
contrast.mat <- limma::makeContrasts(delta - gamma,
                                     levels = design)
```

Now `limma` can perform the transformation with `voom`, fit the model, compute the contrasts and compute test statistics with `eBayes`:

```R
vm <- limma::voom(dge, design = design, plot = TRUE)
fit <- limma::lmFit(vm, design = design)
fit.contrasts <- limma::contrasts.fit(fit, contrast.mat)
fit.contrasts <- limma::eBayes(fit.contrasts)
```

We can use `topTable` to get the most significantly differentially expressed genes:
```R
limma::topTable(fit.contrasts, number = 10, sort.by = "P")
```

And we can check whether this corresponds to the counts by generating a violin plot:

```R
Seurat::VlnPlot(pancreas.dg, "PPY", split.by = "tech")
Seurat::VlnPlot(pancreas.dg, "RBP4", split.by = "tech")
```
