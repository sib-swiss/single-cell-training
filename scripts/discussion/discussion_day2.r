# ------ Day 2

setwd("/export/scratch/twyss/SIB_scRNAseq_course/november2023/data/")

# Import seurat object, already filtered, normalized, scaled.
library(Seurat)

seu<-readRDS("seu_day1.rds")

# Run PCA and plot:
?RunPCA
seu <- Seurat::RunPCA(seu)

?DimPlot
Seurat::DimPlot(seu, reduction = "pca")
# Color PCA plot according to QC variables or gene expression:
Seurat::FeaturePlot(seu, reduction = "pca", features = "percent.globin")
Seurat::FeaturePlot(seu, reduction = "pca", features = "HBA1") # HBA1=hemoglobin subunit alpha 1
Seurat::FeaturePlot(seu, reduction = "pca", features = c("IL32", "TUBB"))

# Heatmap of gene expression of top genes contributing to each of the 12 first PCs:
Seurat::DimHeatmap(seu, dims = 1:12, cells = 500, balanced = TRUE)

# Elbowplot:
Seurat::ElbowPlot(seu, ndims = 40)

# If you want to run T-SNE:
?RunTSNE

# Run UMAP:
seu <- Seurat::RunUMAP(seu, dims = 1:25)

# Plot UMAP, by default coloring by sample id:
Idents(seu)
Seurat::DimPlot(seu, reduction = "umap")

# Exercises
# A. Color the dots in the UMAP according to a variable (e.g. percent.globin or HBA1). 
# Any idea where the erythrocytes probably are in the UMAP?
Seurat::FeaturePlot(seu, features = c("HBA1", 
                                      "percent.globin", 
                                      "IGKC", # top variable gene, immunoglobin kappa constant (antibody, plasma B cells)
                                      "percent.mito"))

FeaturePlot(seu,features=c("HBA1","percent.globin"), blend=TRUE)

# B. Change the number of neighbors used for the calculation of the UMAP. 
# Which is the parameter to change and how did it affect the output.
# What is the default ? In which situation would you lower/increase this ?
?RunUMAP
seu <- Seurat::RunUMAP(seu, dims = 1:25, n.neighbors = 5)
library(ggplot2)
Seurat::DimPlot(seu, reduction = "umap")  + labs(title="N-neighbors=5")

# C. The number of dims to extremes dims = 1:5 or dims = 1:50 how did it affect the output ?
# In your opinion better few PCAs too much or too few ? Why does dims = 1:100 not work? 
# When would more precision be needed?
seu <- Seurat::RunUMAP(seu, dims = 1:5)
Seurat::DimPlot(seu, reduction = "umap")

# Check if your small cluster of IGKC+ cells is still well defined:
Seurat::FeaturePlot(seu, features = "IGKC")
seu <- Seurat::RunUMAP(seu, dims = 1:50) 
Seurat::DimPlot(seu, reduction = "umap")

# Go back to a UMAP with 30 neighbors and 25 dims
seu <- Seurat::RunUMAP(seu, dims = 1:25)
Seurat::DimPlot(seu, reduction = "umap")

# ----- Integration: only if necessary (visual separation of samples on UMAP):

Seurat::DimPlot(seu, reduction = "umap")

# Method using functions from V4 still work
# Each sample has to be normalized independently before batch correction.
# Create list of seurat objects:
seu_list <- Seurat::SplitObject(seu, split.by = "orig.ident")

for (i in 1:length(seu_list)) {
  seu_list[[i]] <- Seurat::NormalizeData(seu_list[[i]])
  seu_list[[i]] <- Seurat::FindVariableFeatures(seu_list[[i]], selection.method = "vst", nfeatures = 2000,
                                                verbose = FALSE)
}

# Find anchors and integrate:
seu_anchors <- Seurat::FindIntegrationAnchors(object.list = seu_list, dims = 1:30)
seu_int <- Seurat::IntegrateData(anchorset = seu_anchors, dims = 1:30, 
                                 normalization.method = "LogNormalize")

names(seu_int@assays)
DefaultAssay(seu_int) # integrated
# Seurat v5: changed structure of RNA assay counts after integration...
seu_int@assays$RNA

# Need to re-scale, re-run PCA and UMAP
seu_int <- Seurat::ScaleData(seu_int)
seu_int <- Seurat::RunPCA(seu_int, npcs = 30)

seu_int <- Seurat::RunUMAP(seu_int, reduction = "pca", dims = 1:25)

Seurat::DimPlot(seu, reduction = "umap")
Seurat::DimPlot(seu_int, reduction = "umap")

saveRDS(seu_int, "seu_int_day2_part1.rds")

# Add a dimensional reduction to a Seurat object:
# the coordinates of each cell in the dim red. are stored in a data.frame in the seurat object:
head(seu_int@reductions$umap@cell.embeddings)

identical(rownames(seu@reductions$umap@cell.embeddings),
          rownames(seu_int@reductions$umap@cell.embeddings))

seu_int[["umap_noInt"]] <- Seurat::CreateDimReducObject(embeddings = seu@reductions$umap@cell.embeddings, 
                                                key = "UMAPnoInt_", assay = DefaultAssay(seu))
names(seu_int@reductions)

p1 <- Seurat::DimPlot(seu_int, reduction = "umap_noInt") + ggtitle("Before integration")
p2 <- Seurat::DimPlot(seu_int) + ggtitle("After integration")
cowplot::plot_grid(p1, p2, ncol = 2)

# New function for integration in Seurat V5:
seu_split <- seu
seu_split[["RNA"]] <- split(seu_split[["RNA"]], f = seu_split$orig.ident) # split RNA assay 
# according to batches:
# Because the data is split into layers, normalization and variable feature 
# identification is performed for each batch independently
seu_split
seu_split <- Seurat::NormalizeData(seu_split)
seu_split <- Seurat::FindVariableFeatures(seu_split)
seu_split <- Seurat::ScaleData(seu_split)
seu_split <- Seurat::RunPCA(seu_split)

seu_split <- Seurat::IntegrateLayers(
  object = seu_split, method = CCAIntegration,
  orig.reduction = "pca", new.reduction = "integrated.cca",
  verbose = FALSE
)
# Anchor-based CCA integration (method=CCAIntegration)
# Anchor-based RPCA integration (method=RPCAIntegration)
# Harmony (method=HarmonyIntegration)
# FastMNN (method= FastMNNIntegration)
# scVI (method=scVIIntegration)

seu_split <- Seurat::RunUMAP(seu_split, reduction = "integrated.cca", dims = 1:30, 
                             reduction.name = "umap.cca")
Seurat::DimPlot(seu_split, reduction = "umap.cca", group.by = "orig.ident", label=T)

# ---- Clustering at several resolutions:
head(seu_int@meta.data)

?FindNeighbors
?FindClusters
seu_int <- Seurat::FindNeighbors(seu_int, dims = 1:25)
seu_int <- Seurat::FindClusters(seu_int, resolution = seq(0.1, 0.8, by=0.1))

# Cluster ID of each cell is added to meta.data with the last
# resolution used as default ID:
head(seu_int@meta.data)
head(Idents(seu_int))
Idents(seu_int) <- "orig.ident"
DimPlot(seu_int)

DimPlot(seu_int,group.by="orig.ident")

# subdivision of clusters:
library(clustree)
clustree::clustree(seu_int@meta.data[,grep("integrated_snn_res", colnames(seu_int@meta.data))],
                   prefix = "integrated_snn_res.")

Seurat::Idents(seu_int)
Seurat::DimPlot(seu_int)
Seurat::DimPlot(seu_int, group.by = "integrated_snn_res.0.1")

# Which res. to choose? (use clustree also)
Seurat::DimPlot(seu_int, group.by = "integrated_snn_res.0.3")

# check distribution of cells per sample per cluster:
table(seu_int$integrated_snn_res.0.3,
      seu_int$orig.ident)

saveRDS(seu_int, "seu_int_day2_part2.rds")

####  Manual and automatic cell annotation:

seu_int<-readRDS("seu_int_day2_part2.rds")

# Manual with addModuleScore:

# Keep res 0.3 as the default cell ID:  
Idents(seu_int)<-"integrated_snn_res.0.3"
seu_int <- Seurat::SetIdent(seu_int, value = seu_int$integrated_snn_res.0.3)
# go back to RNA assay as default:
DefaultAssay(seu_int) <- "RNA"

# Manual annotation with marker genes:
tcell_genes <- c("IL7R", "LTB", "TRAC", "CD3D")
monocyte_genes <- c("CD14", "CST3", "CD68", "CTSS")

# calculate a score per cell for a group of genes:

Seurat::FeaturePlot(seu_int, tcell_genes, ncol=2)
# violin plots of distribution of gene expression per cluster:
Seurat::VlnPlot(seu_int,
                features = tcell_genes,
                ncol = 2) # cluster 0 and 8
# if you still doubt about which clustering resolution to retain:
Seurat::VlnPlot(seu_int,
                features = tcell_genes,
                group.by = "integrated_snn_res.0.1",
                ncol = 2)

Seurat::VlnPlot(seu_int,
                features = tcell_genes,
                group.by = "integrated_snn_res.0.1",
                split.by="orig.ident",
                ncol = 2)

# Monocyte genes:
Seurat::FeaturePlot(seu_int, monocyte_genes, ncol=2)
Seurat::VlnPlot(seu_int,
                features = monocyte_genes,
                ncol = 2) # cluster 2
DefaultAssay(seu_int)<-"RNA"
# calculate a score per cell for a group of genes:

# This works with Seurat v4:
# seu_int <- Seurat::AddModuleScore(seu_int,
#                       features = list(tcell_genes),
#                       name = "tcell_genes")

# need to join layers before running AddModuleScore
# this will probably be fixed in next releases
seu_int[["joined"]] <- JoinLayers(seu_int[["RNA"]])

seu_int <- Seurat::AddModuleScore(seu_int,
                                  features = list(tcell_genes),
                                  name = "tcell_genes",
                                  assay = "joined")

head(seu_int@meta.data)
Seurat::FeaturePlot(seu_int, "tcell_genes1", ncol=2)
Seurat::VlnPlot(seu_int,
                features = "tcell_genes1",
                ncol = 2) # cluster 0 and 8

# Use UCell for scoring gene signatures in single cells:

# remotes::install_github("carmonalab/UCell", ref="v1.3")
# For R version 4.2:
# BiocManager::install("UCell")
library(UCell)
set.seed(123)
signatures <- list(Immune = c("PTPRC"), 
                   Macrophage = c("CTSB", "C1QB", "LAPTM5",
                                   "TYROBP", "PSAP", "C1QA", "HLA-DRA", "CTSD", "NPC2", "FCER1G"), 
                   Tcell = c("CD3D","CD3E", "CD3G", "CD2"), 
                   Bcell = c("MS4A1", "CD79A", "CD79B", "CD19", "BANK1"),
                   Myeloid_cell = c("CD14", "LYZ", "CSF1R", "FCER1G", "SPI1", "LCK-"))

seu_int <- UCell::AddModuleScore_UCell(seu_int, features = signatures, name = NULL,
                                      ncores = 4, assay = "joined")
head(seu_int@meta.data)
Seurat::FeatureScatter(seu_int, feature1 = "tcell_genes1", feature2 = "Tcell",
                       group.by = "integrated_snn_res.0.3")


# identify cycling cells:
# extract the built-in genes for cell cycling:

s.genes <- Seurat::cc.genes.updated.2019$s.genes
g2m.genes <- Seurat::cc.genes.updated.2019$g2m.genes

# add cycling state to meta.data:
seu_int <- Seurat::CellCycleScoring(seu_int,
                                    s.features = s.genes,
                                    g2m.features = g2m.genes,
                                     assay="joined")

Seurat::DimPlot(seu_int, group.by = "Phase")

# Automatic annotation using reference gene expression profiles (with SingleR)
library(celldex)
library(SingleR)
library(dittoSeq) 

?NovershternHematopoieticData
# Download and cache the normalized expression values of 211 bulk human microarray 
# samples of sorted hematopoietic cell populations that can be found in GSE24759.
ref <- celldex::NovershternHematopoieticData()
# /home/twyss/.cache/R/ExperimentHub
# does not exist, create directory? (yes/no): yes
class(ref)
?SummarizedExperiment
table(ref$label.main)

?SingleR
# already ran:
seu_int_SingleR <- SingleR::SingleR(test = Seurat::GetAssayData(seu_int, assay = "joined"),
                                    ref = ref,
                                    labels = ref$label.main)
# the output contains annotation scores and assigned labels per cell:
head(seu_int_SingleR)
dim(colData(ref))

dim(seu_int_SingleR)

# assess cell type annotation
SingleR::plotScoreHeatmap(seu_int_SingleR)
SingleR::plotDeltaDistribution(seu_int_SingleR)

# remove low frequency cells from annotation
singleR_labels <- seu_int_SingleR$labels
t <- table(singleR_labels) # count the total number of cells per cell type
t
other <- names(t)[t < 10]
other
# assign to NA:
singleR_labels[singleR_labels %in% other] <- "none" # or NA as in course website

# add the singleR annotation to the integrated seurat's object meta.data:
seu_int$SingleR_annot <- singleR_labels

# UMAP with Seurat or dittoSeq:
Seurat::DimPlot(seu_int, group.by = "SingleR_annot")
dittoSeq::dittoDimPlot(seu_int, "SingleR_annot", size = 0.7) # can take Seurat, SingleCellExperiment, or SummarizedExperiment object.

# barplot of cell type numbers per sample:
dittoSeq::dittoBarPlot(seu_int, var = "SingleR_annot", group.by = "orig.ident")

# proportion of cell types per cluster at res. 0.3 (cluster 0 and 8 could be T cells)
dittoSeq::dittoBarPlot(seu_int, 
                       var = "SingleR_annot", 
                       group.by = "integrated_snn_res.0.3")

# You can export UMAP coordinates and cell clustering if needed:
write.csv(cbind(cell=rownames(seu_int@reductions$umap@cell.embeddings),
                seu_int@reductions$umap@cell.embeddings),
          "UMAP_coord.csv", row.names = F, quote = F)
write.csv(cbind(cell=rownames(seu_int@meta.data),
                res_0.3=seu_int$integrated_snn_res.0.3,
                singleR=seu_int$SingleR_annot),
          "cell_clustering_ID.csv", row.names = F, quote = F)

saveRDS(seu_int, "seu_int_day2_part2.rds")


# save data set and clear environment:
rm(list = ls())
gc()
.rs.restartR()


