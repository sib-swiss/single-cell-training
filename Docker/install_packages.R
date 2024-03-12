#!/usr/bin/env Rscript

# R script to install requirements for exercises -------------------------------

## a vector of packages to install (edit in this section) ----------------------
### packages could be either on CRAN or bioconductor

# force compilation from source for tidytree
# binary installation is only available for earlier versions


pkgs <- c("devtools", "remotes",
  "BiocManager", "dplyr",
  "scuttle",
  "cowplot",
  "ggbeeswarm", "ggnewscale", "msigdbr", "ggrastr",
  "bit64", "scater",
  "AnnotationDbi",
  "SingleR", "celldex",
  "dittoSeq", "DelayedArray",
  "DelayedMatrixStats",
  "limma", "SingleCellExperiment",
  "SummarizedExperiment",
  "slingshot", "batchelor",
  "clustree", "edgeR", "org.Hs.eg.db",
  "glmGamPoi"
)

## install Bioconductor --------------------------------------------------------
if (!require("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

## install and check package loading -------------------------------------------
for (pkg in basename(pkgs)) {
  BiocManager::install(pkg, ask = FALSE, update = FALSE)

  if (!library(pkg, character.only = TRUE, logical.return = TRUE)) {
    write(
      paste0(
        "Installation of package ",
        pkg,
        " exited with non-zero exit status"
      ),
      stdout()
    )
    quit(status = 1, save = "no")
  }
}
# Seurat from source
install.packages("Seurat", repos = "https://cran.rstudio.com", type = "source")

# clusterProfiler from repo
remotes::install_github("YuLab-SMU/yulab.utils")
remotes::install_github("YuLab-SMU/clusterProfiler")

# # installation old Matrix.utils because it's not on CRAN
# install.packages("https://cran.r-project.org/src/contrib/Archive/Matrix.utils/Matrix.utils_0.9.8.tar.gz", type = "source", repos = NULL)
# # adding this because Matrix might be downgraded by above command
# install.packages("Matrix")
# Monocle3:
devtools::install_github("cole-trapnell-lab/leidenbase")
devtools::install_github("cole-trapnell-lab/monocle3")
# Presto for fast wilcoxon rank sum test:
remotes::install_github("immunogenomics/presto")
# STACAS for integration:
remotes::install_github("carmonalab/STACAS")
