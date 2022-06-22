#!/usr/bin/env Rscript

# R script to install requirements for exercises -------------------------------

## a vector of packages to install (edit in this section) ----------------------
### packages could be either on CRAN or bioconductor

pkgs <- c("ggplot2", "BiocManager", "sctransform",
                   "devtools", "cowplot", "matrixStats",
                   "ggbeeswarm", "ggnewscale", "msigdbr",
                   "Seurat", "bit64", "Matrix.utils", "scater",
                   "AnnotationDbi",
                    "SingleR", "clusterProfiler", "celldex",
                    "dittoSeq", "DelayedArray",
                    "DelayedMatrixStats",
                    "limma", "SingleCellExperiment",
                    "SummarizedExperiment",
                    "slingshot", "batchelor",
                    "clustree", "edgeR")

### if packages need to be installed from github:
### devtools::install_github("namespace/repo")

## install Bioconductor --------------------------------------------------------
if (!require("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}

## install and check package loading -------------------------------------------
for (pkg in basename(pkgs)) {
    BiocManager::install(pkg, ask = FALSE, update = FALSE)

    if (! library(pkg, character.only = TRUE, logical.return = TRUE)) {
        write(paste0("Installation of package ",
                     pkg,
                     " exited with non-zero exit status"),
                     stdout())
        quit(status = 1, save = "no")
    }
}

# Monocle3:
devtools::install_github("cole-trapnell-lab/leidenbase")
devtools::install_github("cole-trapnell-lab/monocle3")
