# Packages to install:

# CRAN packages:
install.packages(c("ggplot2", "BiocManager", "Seurat", "sctransform", "devtools", "cowplot",
                   "ggbeeswarm"))

# Bioconductor packages:
BiocManager::install(c("scran", "clusterProfiler", "scater", "SingleR",
                       "celldex","BiocGenerics", "DelayedArray", "DelayedMatrixStats",
                       "limma", "S4Vectors", "SingleCellExperiment",
                       "SummarizedExperiment", "batchelor", "Matrix.utils",
                       "slingshot",
                       "clustree"))

# Monocle3:
devtools::install_github('cole-trapnell-lab/leidenbase')
devtools::install_github('cole-trapnell-lab/monocle3')
