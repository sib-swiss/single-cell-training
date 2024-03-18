# SIB's scRNAseq course
# 8724
# R version 4.3.2

# Correction of exercises

#####  ---- Day 1

# First step: cellranger count
# to generate filtered feature matrix counts
# (Use Terminal !)
# cellranger count --help

# --id = output folder name (user defined)
# --sample = Prefix of the filenames of FASTQs to select


# cellranger count \
# --id=ETV6-RUNX1_1 \
# --sample=ETV6-RUNX1_1 \
# --transcriptome=/data/cellranger_index \
# --fastqs=/home/rstudio/single_cell_course/course_data/reads \
# --localcores=4 


#### Analysis within R (R console)
setwd("/export/scratch/twyss/SIB_scRNAseq_course/March2024/data/")

# ----- Import cellranger output and generate Seurat object
library(Seurat) # v.5.0.2
library(ggplot2)
library(Matrix) # 1.6.5

# import Sample info and path to cellranger output for each sample:
?Read10X

sample_info <- read.csv("course_data/sample_info_course.csv")
View(sample_info)
datadirs <- file.path("course_data", "count_matrices", sample_info$SampleName,
                      "outs", "filtered_feature_bc_matrix")
names(datadirs) <- gsub("_", "-", sample_info$SampleName)
datadirs <- datadirs[1:3] # only start with 3 replicates of PBMMCs
datadirs 

# generate count matrix for these 3 samples
sparse_matrix <- Seurat::Read10X(data.dir = datadirs)
dim(sparse_matrix)
# [1] 33694  6946 # genes x cells

# each cell gets appended a prefix with SampleName (i.e names of vector datadirs)
head(colnames(sparse_matrix))
# [1] "PBMMC-1_AAACCTGCAGACGCAA-1" "PBMMC-1_AAACCTGTCATCACCC-1" "PBMMC-1_AAAGATGCATAAAGGT-1"
# [4] "PBMMC-1_AAAGCAAAGCAGCGTA-1" "PBMMC-1_AAAGCAACAATAACGA-1" "PBMMC-1_AAAGCAACATCAGTCA-1" 

tail(colnames(sparse_matrix))
# "PBMMC-3_TTTGTCAAGTACGTTC-1" "PBMMC-3_TTTGTCACAATGAAAC-1" "PBMMC-3_TTTGTCACATCTGGTA-1"
# [4] "PBMMC-3_TTTGTCAGTACAGCAG-1" "PBMMC-3_TTTGTCAGTGTGAATA-1" "PBMMC-3_TTTGTCAGTTCTCATT-1"

# check counts of 3 genes (within rownames of sparse_matrix):
sparse_matrix[c("PECAM1", "CD8A", "TSPAN1"), 1:30]

# create Seurat object with sparse_matrix:
?Seurat::CreateSeuratObject
# Do some low filtering: genes expressed in at least 3 cells,
# cells with at least 100 genes:
seu <- Seurat::CreateSeuratObject(counts = sparse_matrix,
                                  project = "pbmmc",
                                  min.cells = 3,
                                  min.features = 100)
seu # 
# 18673 features across 6946 samples within 1 assay 

# Checking the structure of the object:
View(seu)
# str(seu)
table(seu@active.ident)
head(seu@meta.data)
class(seu@meta.data)

# Variables stored in the columns of the meta.data can be used as normal variables
# in a data frame, to plot, etc
# Access meta.data variables either using $ or @meta.data$
hist(seu$nCount_RNA, 
     xlab="", 
     main="Number of reads per cell",
     col="grey")
hist(seu@meta.data$nCount_RNA,
     xlab="", 
     main="Number of reads per cell",
     col="lightblue")

# Seurat's built-in functions to create scatter plots of 1 variable versus the other
# ! Use the column names of seu@meta.data:
colnames(seu@meta.data)
Seurat::FeatureScatter(seu, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

# Simple QC plot of number of counts and detected genes per sample (by default identity
# available in meta.data)
Seurat::VlnPlot(seu, features = c("nCount_RNA",
                                  "nFeature_RNA"))
# check for ribosomal, mitochondrial and hemoglobin genes:
# mitochondrial genes
seu <- Seurat::PercentageFeatureSet(seu, 
                                    pattern = "^MT-", 
                                    col.name = "percent.mito")

# ribosomal genes
seu <- Seurat::PercentageFeatureSet(seu, 
                                    pattern = "^RP[SL]",
                                    col.name = "percent.ribo")

# hemoglobin genes (but not HBP)
seu <- Seurat::PercentageFeatureSet(seu,
                                    pattern = "^HB[^(P)]",
                                    col.name = "percent.globin")
# columns with percentage of each gene type added to meta.data as new columns:
head(seu@meta.data)

# generate violin plots:
Seurat::VlnPlot(seu, features = c("percent.mito",
                                  "percent.ribo",
                                  "percent.globin"))
# check negative correlation between hemoglobin genes and ribosomal genes:
Seurat::FeatureScatter(seu, 
                       feature1 = "percent.globin", 
                       feature2 = "percent.ribo")

# generate custom boxplot of most expressed genes across all cells:
most_expressed_boxplot <- function(object, ngenes = 20){ # 2 arguments, 1 seurat object and 
  # default of number of genes to plot
  
  # matrix of raw counts
  # cts <- Seurat::GetAssayData(object, assay = "RNA", slot = "counts")
   cts <- Seurat::GetAssayData(object, assay = "RNA", layer = "counts")
  
   # get percentage/cell
  cts <- t(cts)/colSums(cts)*100
  medians <- apply(cts, 2, median)
  
  # get top n genes
  most_expressed <- order(medians, decreasing = T)[ngenes:1]
  most_exp_matrix <- as.matrix((cts[,most_expressed]))
  
  # prepare for plotting
  most_exp_df <- stack(as.data.frame(most_exp_matrix))
  colnames(most_exp_df) <- c("perc_total", "gene")
  
  # boxplot with ggplot2
  boxplot <- ggplot(most_exp_df, aes(x=gene, y=perc_total)) +
    geom_boxplot() +
    coord_flip()
  return(boxplot)
}

most_expressed_boxplot(seu, 25)

# Decide on thresholds to filter out cells (as in publication of the data):
seu <- subset(seu, subset = nFeature_RNA > 200 & 
                nFeature_RNA < 5000 &
                percent.mito < 8)
seu

# Check distributions again:
Seurat::VlnPlot(seu, features = c("nFeature_RNA",
                                  "percent.mito"))

#### Normalization and scaling:
# check data before normalization:
Seurat::GetAssayData(seu)[1:10,1:10]  

# normalize (added to seurat object):
?Seurat::NormalizeData # log(1+ (count/colSums * scale.factor)) # natural logarithm
seu <- Seurat::NormalizeData(seu,
                             normalization.method = "LogNormalize",
                             scale.factor = 10000)
# check data after normalization:
Seurat::GetAssayData(seu, layer="data")[1:10,1:10]  

# Find variable genes (for PCA):
seu <- Seurat::FindVariableFeatures(seu,
                                    selection.method = "vst",
                                    nfeatures = 2000)
# what are the top 10 most variable genes?
top10 <- head(Seurat::VariableFeatures(seu), 10)
top10

# plot as scatterplot of average expression level (x-axis) versus variance (y-axis)
vf_plot <- Seurat::VariableFeaturePlot(seu) # store unlabeled plot in object
vf_plot
Seurat::LabelPoints(plot = vf_plot,
                    points = top10, repel = TRUE) # add gene labels of top variable genes

# scaling (for PCA), only the variable genes are scaled by default
seu <- Seurat::ScaleData(seu,
                         features = rownames(seu))

# Testing SCTransform (takes a few minutes)
if(FALSE) {
 seu_SCT <- Seurat::SCTransform(seu)
saveRDS(seu_SCT, "seu_SCT_day1.rds")
# import seurat object with SCTransform Day1:
} else {
seu_SCT<-readRDS("seu_SCT_day1.rds")
}
names(seu_SCT@assays)
# [1] "RNA" "SCT"

DefaultAssay(seu_SCT) <- "RNA"

# Save the seurat object for future use
saveRDS(seu, "seu_day1.rds")

# ! Clear your environment, as to not overload the online resources:

rm(list = ls())
gc()
.rs.restartR()


#### Bonus: DoubletFinder: only works with Seurat v4!!
library(DoubletFinder) # remotes::install_github('chris-mcginnis-ucsf/DoubletFinder')

# https://bioconductor.org/books/release/OSCA/doublet-detection.html

seurat_functions<-function(counts) {
  so <- CreateSeuratObject(counts = counts,
                           project = "pbmc", min.cells = 0, min.features = 100)
  so <- NormalizeData(object = so,
                      normalization.method = "LogNormalize", scale.factor = 10000)                 
  so <-FindVariableFeatures(so, selection.method = "vst",         
                            nfeatures = 2000, verbose = FALSE) 
  so <- ScaleData(so)
  so <- RunPCA(so, npcs = 30)
  print(ElbowPlot(so))
  so <- RunUMAP(so, dims=1:19)
  so<-FindNeighbors(so)
  so<-FindClusters(so, resolution=0.1)
  return(so)
}

# For Doublet Finder, generate a seurat object for each sample:
PBMMC1_so<-seurat_functions(Read10X("course_data/count_matrices/PBMMC_1/outs/filtered_feature_bc_matrix/"))
# PBMMC2_so<-seurat_functions(Read10X("course_data/count_matrices/PBMMC_2/outs/filtered_feature_bc_matrix/"))
# ...

# https://kb.10xgenomics.com/hc/en-us/articles/360001378811-What-is-the-maximum-number-of-cells-that-can-be-profiled-
percent.doublet<-0.016

# https://nbisweden.github.io/workshop-scRNAseq/labs/compiled/seurat/seurat_01_qc.html#Predict_doublets
# set a percentage of doublet rate at 1/1000 of the number of cells per sample:
sweep.res.list_pbmc1 <- paramSweep_v3(PBMMC1_so, PCs = 1:19, sct = FALSE)
sweep.stats_pbmc1 <- summarizeSweep(sweep.res.list_pbmc1, GT = FALSE)
bcmvn_pbmc1 <- find.pK(sweep.stats_pbmc1)
par(mar=c(6,6,6,6))
barplot(bcmvn_pbmc1$BCmetric, names.arg =
          bcmvn_pbmc1$pK, las=2)
## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
homotypic.prop <- modelHomotypic(PBMMC1_so$RNA_snn_res.0.1)           ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
nExp_poi <- round(percent.doublet*nrow(PBMMC1_so@meta.data))  ## Assuming 7% doublet formation rate - tailor for your dataset
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
PBMMC1_so <- doubletFinder_v3(PBMMC1_so, PCs = 1:10, pN = 0.24, pK = 0.005, 
                              nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
head(PBMMC1_so@meta.data)
PBMMC1_so <- doubletFinder_v3(PBMMC1_so, PCs = 1:10, pN = 0.24, pK = 0.005, 
                              nExp = nExp_poi.adj, reuse.pANN = "pANN_0.24_0.005_26", sct = FALSE)
head(PBMMC1_so@meta.data)

table(PBMMC1_so$DF.classifications_0.24_0.005_26,
      PBMMC1_so$DF.classifications_0.24_0.005_21)
###.       Doublet Singlet
#  Doublet      21      5
#  Singlet       0    1586

# Check number of expressed genes in singlets and doublets:
DimPlot(PBMMC1_so, group.by = "DF.classifications_0.24_0.005_26")


