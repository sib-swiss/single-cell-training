#!/usr/bin/env Rscript

# CRAN packages:
cran_packages <- c("ggplot2", "BiocManager", "sctransform",
                   "devtools", "cowplot", "matrixStats",
                   "ggbeeswarm", "ggnewscale", "msigdbr",
                   "Seurat", "bit64", "Matrix.utils", "scater")

for (p in cran_packages) {
    install.packages(p, repos = "http://cran.us.r-project.org");

    if (! library(p, character.only = TRUE, logical.return = TRUE)) {
        write(paste0("Installation of package ",
                     p,
                     " exited with non-zero exit status"),
                     stdout())
        quit(status = 1, save = "no")
    }
}

# Bioconductor packages:
bioc_packages <- c("AnnotationDbi",
                    "SingleR", "clusterProfiler", "celldex",
                    "dittoSeq", "DelayedArray",
                    "DelayedMatrixStats",
                    "limma", "SingleCellExperiment",
                    "SummarizedExperiment",
                    "slingshot", "batchelor",
                    "clustree", "edgeR")

for (p in bioc_packages) {
    BiocManager::install(p, ask = FALSE);

    if (! library(p, character.only = TRUE, logical.return = TRUE)) {
        write(paste0("Installation of package ",
                     p,
                     " exited with non-zero exit status"),
                     stdout())
        quit(status = 1, save = "no")
    }
}

# Monocle3:
devtools::install_github("cole-trapnell-lab/leidenbase")
devtools::install_github("cole-trapnell-lab/monocle3")
