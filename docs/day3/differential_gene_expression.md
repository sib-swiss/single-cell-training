
## Material

[:fontawesome-solid-file-pdf: Download the presentation](../assets/pdf/scRNAseq_Day3_DE_022022.pdf){: .md-button }

- More information on [pseudobulk analysis](https://bioconductor.org/books/release/OSCA/multi-sample-comparisons.html#differential-expression-between-condition)
- [Muscat](https://bioconductor.org/packages/devel/bioc/vignettes/muscat/inst/doc/analysis.html) for pseudobulk DGE.
- [Paper](https://www.nature.com/articles/nmeth.4612) on the robustness of different differential expression analysis methods

## Exercises

### Find all markers for each cluster

Load the `seu_int` dataset you have created yesterday:

```R
seu_int <- readRDS("seu_int_day2_part2.rds")
```

And load the following packages (install them if they are missing):

```R
library(Seurat)
library(edgeR) # BiocManager::install("edgeR")
library(limma)
```

The function `FindAllMarkers` performs a Wilcoxon plot to determine the genes differentially expressed between each cluster and the rest of the cells. Other types of tests than the Wilcoxon test are available. Check it out by running `?Seurat::FindAllMarkers`.

Now run analysis:

```R
de_genes <- Seurat::FindAllMarkers(seu_int,  min.pct = 0.25,
                                   only.pos = TRUE)
```

!!! note "Time for coffee"
    This takes a while. Have a break.

Subset the table to only keep the significant genes, and you can save it as a csv file if you wish to explore it further. Then extract the top 3 markers per cluster:

```R
de_genes <- subset(de_genes, de_genes$p_val_adj<0.05)
write.csv(de_genes, "de_genes_FindAllMarkers.csv", row.names = F, quote = F)


library(dplyr)
top_specific_markers <- de_genes %>%
  group_by(cluster) %>%
  top_n(3, avg_log2FC)
```

And generate e.g. a dotplot:

```R
dittoSeq::dittoDotPlot(seu_int, vars = unique(top_specific_markers$gene), 
                       group.by = "integrated_snn_res.0.3")
```

<figure>
    <img src="../../assets/images/dotplot_degenes.png" width="700"/>
</figure>

**Exercise:** What are significant marker genes in cluster 0 and 8? Are the T cell genes in there?

!!! hint
    You can re-load the vector with immune genes with:

    ```R
    tcell_genes <- c("IL7R", "LTB", "TRAC", "CD3D")
    ```

??? done "Answer"
    Running

    ```R
    de_genes[de_genes$gene %in% tcell_genes,]
    ```

    Returns:
    ```
                   p_val avg_log2FC pct.1 pct.2     p_val_adj cluster gene
    CD3D    0.000000e+00  2.0432771 0.768 0.228  0.000000e+00       0 CD3D
    TRAC   1.157793e-287  1.6805543 0.618 0.205 2.161947e-283       0 TRAC
    LTB    1.072643e-266  1.5215326 0.757 0.395 2.002946e-262       0  LTB
    IL7R   5.112493e-211  1.5236657 0.438 0.114 9.546557e-207       0 IL7R
    LTB.7   6.554548e-36  1.1405474 0.672 0.465  1.223931e-31       7  LTB
    TRAC.8 7.568418e-117  1.8472689 0.759 0.273 1.413251e-112       8 TRAC
    CD3D.8 3.079377e-110  1.7144870 0.800 0.326 5.750121e-106       8 CD3D
    LTB.8   1.580808e-61  1.6529117 0.774 0.461  2.951843e-57       8  LTB
    IL7R.2  4.497526e-45  1.1489439 0.458 0.173  8.398231e-41       8 IL7R
    LTB.11  2.727014e-25  0.8193337 0.750 0.467  5.092153e-21      11  LTB
    ```

    So, yes, the T-cell genes are highly significant markers for cluster 0 and 8.

### Differential expression between groups of cells

The `FindMarkers` function allows to test for differential gene expression analysis specifically between 2 groups of cells, i.e. perform pairwise comparisons, eg between cells of cluster 0 vs cluster 2, or between cells annotated as T-cells and B-cells.

First we can set the default cell identity to the cell types defined by `SingleR`:

```R
seu_int <- Seurat::SetIdent(seu_int, value = "SingleR_annot")
```

Run the differential gene expression analysis and subset the table to keep the significant genes:

```R
deg_cd8_cd4 <- Seurat::FindMarkers(seu_int,
                                   ident.1 = "CD8+ T cells",
                                   ident.2 = "CD4+ T cells",
                                   group.by = seu_int$SingleR_annot,
                                   test.use = "wilcox")
deg_cd8_cd4 <- subset(deg_cd8_cd4, deg_cd8_cd4$p_val_adj<0.05)
```

**Exercise:** Are CD8A, CD8B and CD4 in there? What does the sign (i.e. positive or negative) mean in the log fold change values? Are they according to the CD8+ and CD4+ annotations? Check your answer by generating a violin plot of a top differentially expressed gene.

??? done "Answer"
    You can check out the results with:

    ```R
    View(deg_cd8_cd4)
    ```

    For an explanation of the log fold change have a look at `?Seurat::FindMarkers`. At **Value** it says:

    > `avg_logFC`: log fold-chage of the average expression between the two groups. Positive values indicate that the gene is more highly expressed in the first group

    To view CD8A, CD8B and CD4:

    ```R
    deg_cd8_cd4[c("CD4", "CD8A", "CD8B"),]
    ```

    Returning:

    ```
                p_val avg_log2FC pct.1 pct.2    p_val_adj
    CD4  1.070126e-13 -0.4000835 0.012 0.103 1.998246e-09
    CD8A 1.409306e-77  1.2956354 0.344 0.008 2.631597e-73
    CD8B 7.113148e-36  0.8536693 0.479 0.177 1.328238e-31
    ```

    Indeed, because we compared ident.1 = "CD8+ T cells" to ident.2 = "CD4+ T cells", a negative log2FC for the CD4 gene indicates a lower expression in CD8+ T-cells than in CD4+ T-cells, while a positive log2FC for the CD8A and CD8B genes indicates a higher expression in CD8+ T-cells.

    Plotting the genes in these two T-cell groups only:

    ```R
    Seurat::VlnPlot(seu_int, 
                features = c("CD4", "CD8A", "CD8B"),
                idents = c("CD8+ T cells", "CD4+ T cells"))
    ```

    Returning:

    <figure>
      <img src="../../assets/images/violinplot_CD8_CD4.png" width="600"/>
    </figure>

### Differential expression using `limma`

The Wilcoxon test implemented in `FindMarkers` does not allow you to test for complex design (eg factorial experiments) or to include batch as a covariate. It doesn't allow you to run paired-sample T tests for example. 

For more complex designs, we can use `edgeR` or `limma` which are designed for microarray or bulk RNA seq data and provide a design matrix that includes covariates for example, or sample IDs for paired analyses.

We will load an object containing only pro B cells, both from healthy tissues (PBMMC), and malignant tissues (ETV6-RUNX1). 

!!! warning 
    Please NOTE that in the original design of this data set, the healthy and malignant tissues were not  patient-matched, i.e. the real design was not the one of paired healthy and malignant tissues. However, for demonstration purposes, we will show you how to run a paired analysis, and do as if the PBMMC-1 and ETV6-RUNX1-1 samples both came from the same patient 1, the PBMMC-2 and ETV6-RUNX1-2 samples both came from the same patient 2, etc...


We can load the object and explore its UMAP and meta.data like this:

```R
proB <- readRDS("course_data/proB.rds")

DimPlot(proB, group.by = "orig.ident")

table(proB@meta.data$type)
# ETV6-RUNX1      PBMMC 
#      2000       1021

head(proB@meta.data)

```

!!! note 
    If you want to know how this pro-B cell subset is generated, have a look at the script [here](../../assets/scripts/generate_object_proB.R).

Since we will start with differential gene expression, we set the default assay back to "RNA". Also, we set the default identity to the cell type:

```R
Seurat::DefaultAssay(proB) <- "RNA"
Seurat::Idents(proB) <- proB$orig.ident
```

Let's have a look at the UMAP (again), coloured by celltype:

```R
Seurat::DimPlot(proB)
```

Let's say we are specifically interested to test for differential gene expression between the tumor and normal samples.

!!! note
    Here we could also test for e.g. healthy versus diseased within a celltype/cluster.

Now we will run differential expression analysis between cell type *delta* and *gamma* using the technology as a covariate by using `limma`.

Get the count matrix and keep only genes that are expressed in at least one cell (how many are left?):

```R
counts <- Seurat::GetAssayData(proB, slot = "counts")
counts <- counts[rowSums(counts) != 0,]
dim(counts)
```

Generate a `DGEList` object to use as input for `limma`:

```R
dge <- edgeR::DGEList(counts = counts)
dge <- edgeR::calcNormFactors(dge)  
```

Generate a design matrix, including patient ID to model for a paired analysis. If you need help to generate a design matrix, check out the very nice [edgeR User Guide](https://bioconductor.org/packages/release/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf), sections 3.3 and 3.4.
Extract the sample ID from the meta.data, then create the design matrix:

```R
proB$patient.id<-gsub("ETV6-RUNX1", "ETV6_RUNX1", proB$orig.ident)
proB$patient.id<-sapply(strsplit(proB$patient.id, "-"), '[', 2)

design <- model.matrix(~ 0 + type + patient.id , 
                        data = proB@meta.data)

head(design)

# change column names to more simple group names: 
colnames(design)[c(1:2)] <- make.names(c("ETV6-RUNX1", "PBMMC"))

```

Specify which contrast to analyse:

```R
contrast.mat <- limma::makeContrasts(ETV6.RUNX1 - PBMMC,
                                     levels = design)
```

Now `limma` can perform the transformation with `voom`, fit the model, compute the contrasts and compute test statistics with `eBayes`:

```R
vm <- limma::voom(dge, design = design, plot = TRUE)
fit <- limma::lmFit(vm, design = design)
fit.contrasts <- limma::contrasts.fit(fit, contrast.mat)
fit.contrasts <- limma::eBayes(fit.contrasts)
```

We can use `topTable` to get the most significantly differentially expressed genes, and save the full DE results to an object. How many genes are significant? Are you suprised by this number?
```R
limma::topTable(fit.contrasts, number = 10, sort.by = "P")

limma_de <- limma::topTable(fit.contrasts, number = Inf, sort.by = "P")
length(which(limma_de$adj.P.Val<0.05))

```

And we can check whether this corresponds to the counts by generating a violin plot, or a gene downregulated in tumor, or a gene upregulated in tumor:

```R
Seurat::VlnPlot(proB, "S100A9", split.by = "type")
Seurat::VlnPlot(proB, "SOCS2", split.by = "type")
```

We can run a similar analysis with `Seurat`, but this will not take into account the paired design. Run the code below. 

```R
tum_vs_norm <- Seurat::FindMarkers(proB, 
                                   ident.1 = "ETV6-RUNX1", 
                                   ident.2 = "PBMMC", 
                                   group.by = "type")
tum_vs_norm <- subset(tum_vs_norm, tum_vs_norm$p_val_adj<0.05)
```

How many genes are significant? How does the fold change of these genes compare to the fold change of the top genes found by limma?

Keep the `tum_vs_norm` object because we will use this output object for the enrichment analysis in the next section.


??? done "Answer"
    ```R
    dim(tum_vs_norm) 
    ```
    We find 1893 significant genes. If we merge the `FindMarkers` and the `limma` results, keep `limma`'s most significant genes and plot:
    ```R
    merge_limma_FindMarkers <- merge(tum_vs_norm, 
                                     limma_de, 
                                     by="row.names",
                                     all.x=T)
    merge_limma_FindMarkers <- subset(merge_limma_FindMarkers,
                                      merge_limma_FindMarkers$adj.P.Val<0.00001)

    par(mar=c(4,4,4,4))
    plot(merge_limma_FindMarkers$avg_log2FC,
    merge_limma_FindMarkers$logFC,
    xlab="log2FC Wilcoxon", ylab="log2FC limma",
    pch=15, cex=0.5)
    abline(a=0, b=1, col="red")
    ```
    Returning:

    <figure>
      <img src="../../assets/images/limma_vs_wilcoxon.png" width="600"/>
    </figure>