## Learning outcomes

**After having completed this chapter you will be able to:**

- Assign cell cycle phases to a single cell dataset
- Use the package `scater` evaluate cell quality based on reads originating from:
    - mitochondrial genes
    - ribosomal genes
    - dissociation-related genes
- Evaluate confounding effects on expression data by analyzing the explained variance

## Material

[:fontawesome-solid-file-pdf: Download the presentation](../assets/pdf/quality_control.pdf){: .md-button }

- [Paper](https://www.embopress.org/doi/full/10.15252/msb.20209946) on annotating more refined cell stages

## Exercises

> :fontawesome-solid-ribbon: This chapter uses the `gbm` dataset

### Cell cycle analysis

Cells can be captured in different cycling phases, which can be identified.

A list of cell cycle markers is described by:

> Tirosh I, Izar B, Prakadan SM, Wadsworth MH, Treacy D, Trombetta JJ, et al. Dissecting the multicellular ecosystem of metastatic melanoma by single-cell RNA-seq. Science. 2016;352:189–96.

The dataset is directly available from Seurat:

```R
Seurat::cc.genes.updated.2019
```

Extract the genes specific to the S phase and to the G2/M phase:

```R
s.genes <- Seurat::cc.genes.updated.2019$s.genes
g2m.genes <- Seurat::cc.genes.updated.2019$g2m.genes
```

The `CellCycleScoring()` function allows to assign cell cycle phase information
to each cell, stored in the metadata as the "S.Score", "G2M.Score" and "Phase" columns

```R
gbm <- Seurat::CellCycleScoring(gbm,
                                s.features = s.genes,
                                g2m.features = g2m.genes)
```

```R
head(gbm)
table(gbm$Phase)
#   G1  G2M    S
# 2887  711 1493
```

Visualize the distribution of cell cycle markers:

```R
Seurat::RidgePlot(gbm, features = c("PCNA", "MKI67"),
          group.by = "orig.ident",
          ncol = 2)
```

We can proceed with downstream analysis without removing cell cycle for example,
perform clustering, and come back to remove the effect of cell cycle if
we identify a cluster of cells which is
mostly composed of G2/M cells for example.

### Quality control with `scater`

For the below exercises we will use the following packages:

```R
library(scater)
library(SingleCellExperiment)
```

Scater includes different types of quality controls:

1. QC and filtering of cells
2. QC and filtering of features (genes)
3. QC of experimental variables

the scater package as well as other bioconductor packages, rely on
an object of the class `SingleCellExperiment`:

```R
cts <- Seurat::GetAssayData(gbm, slot = "counts")

gbm_sce <- SingleCellExperiment::SingleCellExperiment(
  assays = list(counts = cts),
  colData = gbm@meta.data,
  rowData = rownames(gbm)
)
class(gbm)
class(gbm_sce)
gbm_sce
```

#### Per cell QC

We can again check percent mitochondrial gene expression as well as dissociation protocol-related gene expression. Genes associated with the dissociation protocol, i.e. stress genes, can sometimes cause clustering of stressed cells apart from the other cells.
They have been described by:

> Van Den Brink SC, Sage F, Vértesy Á, Spanjaard B, Peterson-Maduro J, Baron CS, et al. Single-cell sequencing reveals dissociation-induced gene expression in tissue subpopulations. Nat Methods. 2017;14:935–6.

And are availabe in the file `data/dissocation_genes.txt`. Load these into a vector, and generate also a vector of ribosomal genes and mitochondrial genes:

```R
dissoc_genes <- readLines("data/gbm_dataset/dissocation_genes.txt")
ribo_genes <- rownames(gbm)[grep(pattern = "^RP[S|L]", rownames(gbm), perl = T)]
mito_genes <- rownames(gbm)[grep(pattern = "^MT-", rownames(gbm))]
```

`scater` calls the `addPerCellQC` function of the `scuttle` package to compute a number of quality control metrics for each cell and feature (i.e gene)

```R
gbm_sce <- scuttle::addPerCellQC(gbm_sce,
                        subsets=list(mito_genes=which(rownames(gbm_sce) %in% mito_genes),
                                     dissoc_genes=which(rownames(gbm_sce) %in% dissoc_genes),
                                     ribo_genes=which(rownames(gbm_sce) %in% ribo_genes)))
```



```R
SingleCellExperiment::colData(gbm_sce)
```

```R
scater::plotColData(gbm_sce, x = "sum", y="detected")
scater::plotColData(gbm_sce, x = "detected", y="subsets_mito_genes_percent")
scater::plotColData(gbm_sce, x = "detected", y="subsets_dissoc_genes_percent")
scater::plotColData(gbm_sce, x = "subsets_mito_genes_percent", y="subsets_ribo_genes_percent")
```

#### Highly expressed genes

On the gene level, we can look at a plot that shows the top (by default 50) most-expressed genes. Each row in the plot corresponds to a gene; each bar corresponds to the expression of a gene in a single cell; and the circle indicates the median expression of each gene, with which genes are sorted. We expect to see the “usual suspects”, i.e., mitochondrial genes, actin, ribosomal protein, MALAT1. If used, few spike-in transcripts may also be present here, though if all of the spike-ins are in the top 50, it suggests that too much spike-in RNA was added. A large number of pseudo-genes or predicted genes may indicate problems with alignment.

```R
scater::plotHighestExprs(gbm_sce, exprs_values = "counts", n = 30)
```

#### Find explanatory variables

Variable-level metrics are computed by the `getVarianceExplained()` function (after normalization, see below). This calculates the percentage of variance of each gene’s expression that is explained by each variable in the `colData` of the `SingleCellExperiment` object. We can then use this to determine which experimental factors are contributing most to the variance in expression. This is useful for diagnosing batch effects or to quickly verify that a treatment has an effect.

First, computing variance explained on the log-counts,
so that the statistics reflect changes in relative expression.

```R
gbm_sce <- scater::logNormCounts(gbm_sce)  # alternative to Seurat's normalization here using scater
```

In the gbm dataset, we only have 1 patient, so we cannot calculate the effect of experimental variables like sex or donor id, but in case we would have several variables, here is the method with the cell cycle phase as example:

```R
vars <- scater::getVarianceExplained(gbm_sce,
                             variables = "Phase")
head(vars)
```

A distribution of percentage variance explained by each gene is shown, and can indicate whether one or the other experimental variable has high contribution to the variance in the data:

```R
scater::plotExplanatoryVariables(vars)
```

If we think that the cell cycling has an effect on the analysis, and if we want to regress out this effect so that cycling cells are integrated into the rest of the cells and not clustering apart anymore, we can regress out the cell cycling phase at the moment of scaling the data using `ScaleData`. This might be slow to compute.

!!! warning
    The code below is an example, you don't need to run it.

```
# do not run
gbm_cc <- Seurat::ScaleData(gbm, vars.to.regress = c("S.Score", "G2M.Score"))
```

### Save the dataset and clear environment

Now, save the dataset so you can use it tomorrow:

```R
saveRDS(gbm, "gbm_day1.rds")
```

Clear your environment:

```R
rm(list = ls())
gc()
.rs.restartR()
```
