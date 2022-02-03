# Packages to install:

# CRAN packages:
install.packages(c("ggplot2", "BiocManager", "sctransform",
                   "devtools", "cowplot", "matrixStats",
                   "ggbeeswarm", "ggnewscale", "msigdbr"),
                   repos = "http://cran.us.r-project.org")

# Bioconductor packages:
BiocManager::install(c("scran", "scater", "SingleR",
                       "celldex", "BiocGenerics", "DelayedArray",
                       "DelayedMatrixStats",
                       "limma", "S4Vectors", "SingleCellExperiment",
                       "SummarizedExperiment", "batchelor", "Matrix.utils",
                       "slingshot",
                       "clustree"))

# Monocle3:
devtools::install_github("cole-trapnell-lab/leidenbase")
devtools::install_github("cole-trapnell-lab/monocle3")

install.packages(c("Seurat", "bit64"))
BiocManager::install("clusterProfiler")