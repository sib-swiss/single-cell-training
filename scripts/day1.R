# Single-cell RNA sequencing data analysis - SIB course, March 2021


#### ----- Day 1

# Load the necessary packages:
library(BiocManager)
library(Seurat)
library(clustree)
library(cowplot)
library(scran)
library(ggplot2)
library(clusterProfiler)
library(sctransform)
library(scater)
library(SingleR)
library(celldex)
library(devtools)
library(monocle3)
library(slingshot)
library(ggbeeswarm)
library(SingleCellExperiment)
library(DoubletFinder)
library(limma)
library(edgeR)

# set seed
set.seed(1234)

# Set working directory to your own home directory:
setwd("~")

# Import matrices from 10x data.
# origin of data:
# https://support.10xgenomics.com/single-cell-gene-expression/datasets/4.0.0/Parent_SC3v3_Human_Glioblastoma?
# To import this data to the server, on the Terminal:
# mkdir data
# cd data
# wget https://cf.10xgenomics.com/samples/cell-exp/4.0.0/Parent_SC3v3_Human_Glioblastoma/Parent_SC3v3_Human_Glioblastoma_filtered_feature_bc_matrix.tar.gz
# tar -xvf Parent_SC3v3_Human_Glioblastoma_filtered_feature_bc_matrix.tar.gz
# Typical structure of 10x filtered_feature_bc_matrix/:
# barcodes.tsv.gz = list of retained barcodes (i.e. cells)
# features.tsv.gz = list of genes in transcriptome
# matrix.mtx.gz = gene counts per cell matrix

# Seurat provides a generic function to read-in the 10x data by specifying the location of the files
gbm.data <- Read10X(data.dir = "~/data/filtered_feature_bc_matrix/")

# Initialize the Seurat object with the raw (non-normalized data).
gbm <- CreateSeuratObject(counts = gbm.data, project = "gbm",
                          min.cells = 3,
                          min.features = 100)
gbm
# An object of class Seurat
# 24363 features across 5506 samples within 1 assay
# Active assay: RNA (24363 features, 0 variable features)

# Structure of the raw count matrix, print rows that contain
# gene names PECAM1 etc and columns (=cells) 1 to 30:
gbm.data[c("PECAM1", "CD8A", "TSPAN1"), 1:30]

#### Standard pre-processing workflow:
## QC and selecting cells for further analysis
# We typically filter out cells with low (and too high) number of expressed genes, low number of
# detected molecules (i.e. UMI), and high proportion of mitochondrial gene expression.

# The [[ operator can add columns to object metadata, i.e. to add and store QC data
gbm[["percent.mt"]] <- PercentageFeatureSet(gbm, pattern = "^MT-")
# For mouse, the search pattern for mitochondrial genes is "^mt-"

# The same method can be used to check out the percentage of expression of ribosomal
# genes, by searching for the patterns "^RPS" or "^RPL"
ribo_genes<-rownames(gbm)[grep(pattern  = "^RP[S|L]", rownames(gbm), perl = T)]

# The seurat object contains several slots, one of them is the metadata
head(gbm@meta.data)
#                    orig.ident nCount_RNA nFeature_RNA percent.mt
# AAACCCAAGGCGATAC-1        gbm       2225          815  4.7191011
# AAACCCACAAGTCCCG-1        gbm      17882         4760  1.5098982
# AAACCCACAGATGCGA-1        gbm       8172         2319  6.4855605
# AAACCCACAGGTGAGT-1        gbm       9057         3395  3.1577785
# AAACCCAGTCTTGCGG-1        gbm       5612         2716  0.4276550
# AAACCCATCGATAACC-1        gbm       6254         1659  0.7195395

# nCount_RNA = number of UMIs
# nFeature_RNA = number of detected genes
# percent.mt = percent mitochondrial gene expression.

# Visualize QC metrics as a violin plot
VlnPlot(gbm, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(gbm, features="nFeature_RNA",pt.size=0)
VlnPlot(gbm, features="nCount_RNA",pt.size=0)
VlnPlot(gbm, features="percent.mt",pt.size=0)
# We can keep cells that express between 500 and 7500 cells, and that express
# less than 20% mitochondrial genes.
# These cut-offs are specific to each dataset and should be adapted.

# There is often a relationship between number of UMIs, number of detected genes and
# percent.mt

plot1 <- FeatureScatter(gbm, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(gbm, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2


# Subset the Seurat object to only keep retained cells:
gbm <- subset(gbm,
              subset = nFeature_RNA > 500 & nFeature_RNA < 7500 & percent.mt < 20)
# This removed about 500 cells:
gbm
# An object of class Seurat
# 24363 features across 5091 samples within 1 assay
# Active assay: RNA (24363 features, 0 variable features)



### Normalize the data
# After removing unwanted cells from the dataset, the next step is to normalize the data.
# By default, Seurat employs a global-scaling normalization method “LogNormalize” that
# normalizes the feature expression measurements for each cell by the total expression,
# multiplies this by a scale factor (10,000 by default), and log-transforms the result.
# Normalized values are stored in the "RNA" slot of the gbm object.
# Other normalization methods are available, check the help:
?NormalizeData

gbm <- NormalizeData(gbm,
                     normalization.method = "LogNormalize",
                     scale.factor = 10000)

### Using sctransform to normalize the data.
# The current recommendation is to use sctransform to regress out variables instead
# of ScaleData().
# https://doi.org/10.1186/s13059-019-1874-1
## For normalization, here is a bit of code to visually inspect the effect of normalization:
###
# Create a function to plot the mean-variance relationship
plotMeanVar = function( count_mat ){
  gene_attr <- data.frame(mean = rowMeans(count_mat),
                          detection_rate = rowMeans(count_mat > 0),
                          var = apply(count_mat, 1, var))
  gene_attr$log_mean <- log10(gene_attr$mean)
  gene_attr$log_var <- log10(gene_attr$var)
  rownames(gene_attr) <- rownames(count_mat)
  cell_attr <- data.frame(n_umi = colSums(count_mat),
                          n_gene = colSums(count_mat > 0))
  rownames(cell_attr) <- colnames(count_mat)
  g = ggplot(gene_attr, aes(log_mean, log_var)) + geom_point(alpha = 0.3, shape = 16) +
    geom_density_2d(size = 0.3) + geom_abline(intercept = 0, slope = 1, color = "red")
  return(g)
}

# Store the log-normalized and sctransform objects separately to compare:
gbm.norm <- NormalizeData(gbm,
                          normalization.method = "LogNormalize",
                          scale.factor = 10000)
gbm.sct <- SCTransform(object = gbm, verbose = T) # may last a couple of minutes
g1 = plotMeanVar(gbm@assays$RNA@counts) + ggtitle('raw data')
g2 = plotMeanVar(gbm.norm@assays$RNA@data) + ggtitle('log normalization')
g3 = plotMeanVar(gbm.sct@assays$SCT@data) + ggtitle('scTransform')
g1+g2+g3


### Identification of highly variable features (feature selection)
# We next calculate a subset of features that exhibit high cell-to-cell variation in the
# dataset (i.e, they are highly expressed in some cells, and lowly expressed in others).
# Focusing on these genes in downstream analysis helps to highlight biological signal
# in single-cell datasets.
# The procedure in Seurat models the mean-variance relationship inherent in single-cell
# data, and is implemented in the FindVariableFeatures() function.
# By default, 2,000 genes (features) per dataset are returned and these will be used in
# downstream analysis, like PCA.
gbm <- FindVariableFeatures(gbm,
                            selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(gbm), 10)

# plot variable features with labels for top 10 genes
LabelPoints(plot = VariableFeaturePlot(gbm),
            points = top10, repel = TRUE)

### Scaling the data
# Next, we apply scaling, a linear transformation that is a standard pre-processing
# step prior to dimensional reduction techniques like PCA. The ScaleData() function
# 1- shifts the expression of each gene, so that the mean expression across cells is 0
# 2- scales the expression of each gene, so that the variance across cells is 1
# This step gives equal weight in downstream analyses, so that highly-expressed genes do not dominate
# The results of this are stored in gbm[["RNA"]]@scale.data
gbm <- ScaleData(gbm,
                 features = rownames(gbm))
# This step can be used to remove unwanted sources of variance, eg
# mitochondrial gene expression level, number of detected genes, or cell cycle (see later).
# However, the current recommendation is to regress out variables using the
# sctransform method.

#### Cell cycle analysis

# Cells can be captured in different cycling phases, which can be identified.

# A list of cell cycle markers from Tirosh et al
# 2016 (https://science.sciencemag.org/content/352/6282/189) is loaded with Seurat.
cc.genes.updated.2019

# Extract the genes specific to the S phase and to the G2/M phase:
s.genes <- cc.genes.updated.2019$s.genes
g2m.genes <- cc.genes.updated.2019$g2m.genes

# The CellCycleScoring() function allows to assign cell cycle phase information
# to each cell, stored in the metadata as the "S.Score", "G2M.Score" and "Phase" columns
gbm <- CellCycleScoring(gbm, s.features = s.genes, g2m.features = g2m.genes)

head(gbm@meta.data)
table(gbm$Phase)
#   G1  G2M    S
# 2887  711 1493

# Visualize the distribution of cell cycle markers:
RidgePlot(gbm, features = c("PCNA", "MKI67"),
          group.by = "orig.ident",
          ncol = 2)

# We can proceed with downstream analysis without removing cell cycle for example,
# perform clustering, and come back to remove the effect of cell cycle if
# we identify a cluster of cells which is
# mostly composed of G2/M cells for example.

# If we think that the cell cycling has an effect on the analysis,
# and if we want to "remove" this effect so that cycling cells are
# integrated into the rest of the cells and not clustering apart anymore,
# we can regress out the cell cycling phase at the moment of normalizing
# the data using SCTransform.
# This might be slow to compute
# gbm_cc <- SCTransform(gbm, vars.to.regress = c("S.Score", "G2M.Score"))

# Seurat includes the function SCTransform() with the vars.to.regress parameter that
# allows to remove unwanted sources of variance as well.

### Alternative QC exploration using the scater package:
# https://bioconductor.org/packages/release/bioc/vignettes/scater/inst/doc/overview.html
# Scater includes different types of quality controls:
# 1. QC and filtering of cells
# 2. QC and filtering of features (genes)
# 3. QC of experimental variables

# the scater package as well as other bioconductor packages, rely on
# an object type called "SingleCellExperiment":
gbm_sce <- SingleCellExperiment(
  assays = list(counts = gbm@assays$RNA@counts), # use the raw counts
  colData = gbm@meta.data, # add the metadata info for each cell
  rowData = rownames(gbm)
)
gbm_sce

# We can again check percent mitochondrial gene expression as well as dissociation protocol-related
# gene expression.
# Genes associated with the dissociation protocol, i.e. stress genes, can
# sometimes cause clustering of stressed cells apart from the other cells.
# Described by Van Den Brink et al 2017 https://www.nature.com/articles/nmeth.4437#accession-codes
# https://www.nature.com/articles/s41467-019-11036-9#code-availability
## Mouse gene list is from Van Den Brink et al. R code in supplementary method under "In silico purification"
dissoc_genes <- c("Actg1__chr11","Ankrd1__chr19","Arid5a__chr1","Atf3__chr1","Atf4__chr15","Bag3__chr7","Bhlhe40__chr6",
                  "Brd2__chr17","Btg1__chr10","Btg2__chr1","Ccnl1__chr3","Ccrn4l__chr3","Cebpb__chr2","Cebpd__chr16",
                  "Cebpg__chr7","Csrnp1__chr9","Cxcl1__chr5","Cyr61__chr3","Dcn__chr10","Ddx3x__chrX","Ddx5__chr11",
                  "Des__chr1","Dnaja1__chr4","Dnajb1__chr8","Dnajb4__chr3","Dusp1__chr17","Dusp8__chr7",
                  "Egr1__chr18","Egr2__chr10","Eif1__chr11","Eif5__chr12","Erf__chr7","Errfi1__chr4","Fam132b__chr1",
                  "Fos__chr12","Fosb__chr7","Fosl2__chr5","Gadd45a__chr6","Gcc1__chr6","Gem__chr4","H3f3b__chr11",
                  "Hipk3__chr2","Hsp90aa1__chr12","Hsp90ab1__chr17","Hspa1a__chr17","Hspa1b__chr17","Hspa5__chr2",
                  "Hspa8__chr9","Hspb1__chr5","Hsph1__chr5","Id3__chr4","Idi1__chr13","Ier2__chr8","Ier3__chr17",
                  "Ifrd1__chr12","Il6__chr5","Irf1__chr11","Irf8__chr8","Itpkc__chr7","Jun__chr4","Junb__chr8",
                  "Jund__chr8","Klf2__chr8","Klf4__chr4","Klf6__chr13","Klf9__chr19","Litaf__chr16","Lmna__chr3",
                  "Maff__chr15","Mafk__chr5","Mcl1__chr3","Midn__chr10","Mir22hg__chr11","Mt1__chr8","Mt2__chr8",
                  "Myadm__chr7","Myc__chr15","Myd88__chr9","Nckap5l__chr15","Ncoa7__chr10","Nfkbia__chr12","Nfkbiz__chr16",
                  "Nop58__chr1","Nppc__chr1","Nr4a1__chr15","Odc1__chr12","Osgin1__chr8","Oxnad1__chr14","Pcf11__chr7",
                  "Pde4b__chr4","Per1__chr11","Phlda1__chr10","Pnp__chr14","Pnrc1__chr4","Ppp1cc__chr5","Ppp1r15a__chr7",
                  "Pxdc1__chr13","Rap1b__chr10","Rassf1__chr9","Rhob__chr12","Rhoh__chr5","Ripk1__chr13","Sat1__chrX",
                  "Sbno2__chr10","Sdc4__chr2","Serpine1__chr5","Skil__chr3","Slc10a6__chr5","Slc38a2__chr15",
                  "Slc41a1__chr1","Socs3__chr11","Sqstm1__chr11","Srf__chr17","Srsf5__chr12","Srsf7__chr17",
                  "Stat3__chr11","Tagln2__chr1","Tiparp__chr3","Tnfaip3__chr10","Tnfaip6__chr2","Tpm3__chr3",
                  "Tppp3__chr8","Tra2a__chr6","Tra2b__chr16","Trib1__chr15","Tubb4b__chr2","Tubb6__chr18",
                  "Ubc__chr5","Usp2__chr9","Wac__chr18","Zc3h12a__chr4","Zfand5__chr19","Zfp36__chr7","Zfp36l1__chr12",
                  "Zfp36l2__chr17","Zyx__chr6","Gadd45g__chr13","Hspe1__chr1","Ier5__chr1","Kcne4__chr1")
dissoc_genes <- toupper(sapply(dissoc_genes, function(x){
  strsplit(x, "__")[[1]][1]
}))
length(dissoc_genes) # 140
head(dissoc_genes)
# Actg1__chr11 Ankrd1__chr19  Arid5a__chr1    Atf3__chr1   Atf4__chr15    Bag3__chr7
# "ACTG1"      "ANKRD1"      "ARID5A"        "ATF3"        "ATF4"        "BAG3"

# scater calls the addPerCellQC function of the scuttle package to compute a number of quality control metrics
# for each cell and feature (i.e gene)
gbm_sce <- addPerCellQC(gbm_sce,
                        subsets=list(Mito=grep("^MT-", rownames(gbm_sce)),
                                     dissoc_genes=which(rownames(gbm_sce) %in% dissoc_genes),
                                     ribo_genes=which(rownames(gbm_sce) %in% ribo_genes)))
colnames(colData(gbm_sce))
# [1] "orig.ident"                    "nCount_RNA"                    "nFeature_RNA"
# [4] "percent.mt"                    "sum"                           "detected"
# [7] "subsets_Mito_sum"              "subsets_Mito_detected"         "subsets_Mito_percent"
# [10] "total"                         "sum"                           "detected"
# [13] "subsets_Mito_sum"              "subsets_Mito_detected"         "subsets_Mito_percent"
# [16] "subsets_dissoc_genes_sum"      "subsets_dissoc_genes_detected" "subsets_dissoc_genes_percent"
# [19]  "subsets_ribo_genes_sum"        "subsets_ribo_genes_detected"   "subsets_ribo_genes_percent"
# [22] "total"

# QC variables stored in the metadata can be plotted against each other using the plotColData() function
plotColData(gbm_sce, x = "sum", y="detected")

plotColData(gbm_sce, x = "detected", y="subsets_Mito_percent")
plotColData(gbm_sce, x = "detected", y="subsets_dissoc_genes_percent")
plotColData(gbm_sce, x = "subsets_Mito_percent", y="subsets_ribo_genes_percent")


# On the gene level, we can look at a plot that shows the top (by default 50) most-expressed genes.
# Each row in the plot corresponds to a gene; each bar corresponds to the expression of a gene in a single
# cell; and the circle indicates the median expression of each gene, with which genes are sorted. We expect to
# see the “usual suspects”, i.e., mitochondrial genes, actin, ribosomal protein, MALAT1. If used, few spike-in transcripts
# may also be present here, though if all of the spike-ins are in the top 50, it suggests that too much spike-in
# RNA was added. A large number of pseudo-genes or predicted genes may indicate problems with alignment.
plotHighestExprs(gbm_sce, exprs_values = "counts", n = 30)

# Variable-level metrics are computed by the getVarianceExplained() function (after normalization, see below).
# This calculates the percentage of variance of each gene’s expression that is explained by each variable in the
# colData of the SingleCellExperiment object. We can then use this to determine which experimental factors
# are contributing most to the variance in expression. This is useful for diagnosing batch effects or to
# quickly verify that a treatment has an effect.
#
# First, computing variance explained on the log-counts,
# so that the statistics reflect changes in relative expression.
gbm_sce <- logNormCounts(gbm_sce)  # alternative to Seurat's normalization here using scater

# In the gbm dataset, we only have 1 patient, so we cannot calculate the effect
# of experimental variables like sex or donor id, but in case we would have several variables,
# here is the method with the cell cycle phase as example:
vars <- getVarianceExplained(gbm_sce,
                             variables = "Phase"
                             #, variables=c("tissue", "detected", "sex", "age")
)
head(vars)
##                Phase
# AL627309.1 0.03755442
# AL627309.5 0.01721037
# AP006222.2 0.01993929
# LINC01409  0.07512823
# FAM87B     0.01490441
# LINC01128  0.05493752

# A distribution of percentage variance explained by each gene is shown,
# and can indicate whether one or the other experimental variable has high
# contribution to the variance in the data:
plotExplanatoryVariables(vars)

# save the gbm object so that you can import it again tomorrow:
saveRDS(gbm, "gbm.rds")
