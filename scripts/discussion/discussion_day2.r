# ------ Day 2

# setwd("/export/scratch/twyss/SIB_scRNAseq_course/March2024/data/")
setwd("/home/rstudio/single_cell_course/")

# Import seurat object, already filtered, normalized, scaled.
library(Seurat)
library(ggplot2)

seu<-readRDS("seu_day1.rds")

# Run PCA and plot:
?RunPCA
seu <- Seurat::RunPCA(seu)

?DimPlot # see also Value
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

# blend=T: visualize coexpression of 2 features:
FeaturePlot(seu,features=c("HBA1","percent.globin"), blend=TRUE)

# B. Change the number of neighbors used for the calculation of the UMAP. 
# Which is the parameter to change and how did it affect the output.
# What is the default ? In which situation would you lower/increase this ?
?RunUMAP

# Save initial UMAP in a ggplot object
p1<-Seurat::DimPlot(seu, reduction = "umap") + ggtitle("25 dims, N-neighbors=30")

# 5 neighbors:
seu <- Seurat::RunUMAP(seu, dims = 1:25, n.neighbors = 5)
Seurat::DimPlot(seu, reduction = "umap")  + labs(title="N-neighbors=5")
p2<-Seurat::DimPlot(seu, reduction = "umap")  + labs(title="25 dims, N-neighbors=5")

# C. The number of dims to extremes dims = 1:5 or dims = 1:50 how did it affect the output ?
# In your opinion better few PCAs too much or too few ? Why does dims = 1:100 not work? 
# When would more precision be needed?
seu <- Seurat::RunUMAP(seu, dims = 1:5)
Seurat::DimPlot(seu, reduction = "umap")
p3<-Seurat::DimPlot(seu, reduction = "umap")  + labs(title="5 dims, N-neighbors=30")

# Check if your small cluster of IGKC+ cells is still well defined:
Seurat::FeaturePlot(seu, features = "IGKC")
seu <- Seurat::RunUMAP(seu, dims = 1:50) 
Seurat::DimPlot(seu, reduction = "umap")
p4<-Seurat::DimPlot(seu, reduction = "umap")  + labs(title="50 dims, N-neighbors=30")

# Compare all:
cowplot::plot_grid(p1, p2, p3, p4, nrow = 2)

# Go back to a UMAP with 30 neighbors and 25 dims
seu <- Seurat::RunUMAP(seu, dims = 1:25)
Seurat::DimPlot(seu, reduction = "umap")

# save non-batch corrected umap coordinates, to later import
# it back into the seu object after batch correction
umap_noInt<-seu@reductions$umap@cell.embeddings

# ----- Integration: only if necessary (visual separation of samples on UMAP):

Seurat::DimPlot(seu, reduction = "umap")

# Method using functions from V4 still work

# Here, v5 method
# Splitting the RNA assay into layers corresponding to
# different samples:
seu[["RNA"]] <- split(seu[["RNA"]], f = seu$orig.ident)

# Run CCA integration (takes 30 seconds):
seu <- Seurat::IntegrateLayers(object = seu, method = CCAIntegration,
                               orig.reduction = "pca",
                               new.reduction = "integrated.cca",
                               verbose = FALSE)
# Other methods:
# Anchor-based CCA integration (method=CCAIntegration)
# Anchor-based RPCA integration (method=RPCAIntegration)
# Harmony (method=HarmonyIntegration)
# FastMNN (method= FastMNNIntegration)
# scVI (method=scVIIntegration)

# In the RNA assay, there are now 7 layers:
# 1 count and data layer for each sample, and a scale.data layer
View(seu@assays$RNA)
# The layers/samples are integrated and stored in a new reduction
# called "integrated.cca"
View(seu@reductions$integrated.cca@cell.embeddings)

# re-join layers after integration
seu[["RNA"]] <- JoinLayers(seu[["RNA"]])

# Use this new integrated matrix for clustering and visualization. 
# Re-run and visualize the results with UMAP.

# This over-writes the previous UMAP:
seu <- RunUMAP(seu, dims = 1:30, reduction = "integrated.cca")

Seurat::DimPlot(seu, reduction = "umap")

saveRDS(seu, "seu_day2-2.rds")

# Add a dimensional reduction to a Seurat object:
# the coordinates of each cell in the dim red. are stored 
# in a data.frame in the seurat object:
head(seu@reductions$umap@cell.embeddings)

# check that the cells have the same order:
identical(rownames(seu@reductions$umap@cell.embeddings),
          rownames(umap_noInt)) # TRUE

seu[["umap_noInt"]] <- Seurat::CreateDimReducObject(embeddings = umap_noInt, 
                                                key = "UMAPnoInt_", assay = DefaultAssay(seu))
names(seu@reductions)

p1 <- Seurat::DimPlot(seu, reduction = "umap_noInt") + ggtitle("Before integration")
p2 <- Seurat::DimPlot(seu) + ggtitle("After integration")
cowplot::plot_grid(p1, p2, ncol = 2)

  # ---- Clustering at several resolutions:
head(seu@meta.data)

?FindNeighbors
?FindClusters
seu <- Seurat::FindNeighbors(seu, dims = 1:25, reduction = "integrated.cca")
seu <- Seurat::FindClusters(seu, resolution = seq(0.1, 0.8, by=0.1))

# Cluster ID of each cell is added to meta.data with the last
# resolution used as default ID:
head(seu@meta.data)
head(Idents(seu))
Idents(seu) <- "orig.ident"
DimPlot(seu)

DimPlot(seu,group.by="orig.ident")

# subdivision of clusters:
library(clustree)
clustree::clustree(seu@meta.data[,grep("RNA_snn_res", colnames(seu@meta.data))],
                   prefix = "RNA_snn_res.")

Seurat::Idents(seu)
Seurat::DimPlot(seu)
Seurat::DimPlot(seu, group.by = "RNA_snn_res.0.1")

# Which res. to choose? (use clustree also)
Seurat::DimPlot(seu, group.by = "RNA_snn_res.0.1", label=T)
Seurat::DimPlot(seu, group.by = "RNA_snn_res.0.3", label=T)

# check distribution of cells per sample per cluster:
table(seu$RNA_snn_res.0.3,
      seu$orig.ident)

saveRDS(seu, "seu_day2_part2.rds")

####  Manual and automatic cell annotation:

seu<-readRDS("seu_day2_part2.rds")

# Manual with addModuleScore:

# Keep res 0.3 as the default cell ID:  
Idents(seu)<-"RNA_snn_res.0.3"
seu <- Seurat::SetIdent(seu, value = seu$RNA_snn_res.0.3)
# go back to RNA assay as default:
DefaultAssay(seu) <- "RNA"

# Manual annotation with marker genes:
tcell_genes <- c("IL7R", "LTB", "TRAC", "CD3D")
monocyte_genes <- c("CD14", "CST3", "CD68", "CTSS")

# calculate a score per cell for a group of genes:

Seurat::FeaturePlot(seu, tcell_genes, ncol=2)
# violin plots of distribution of gene expression per cluster:
Seurat::VlnPlot(seu,
                features = tcell_genes,
                ncol = 2) # cluster 0 and 8

Seurat::VlnPlot(seu,
                features = tcell_genes,
                split.by="orig.ident",
                ncol = 2)

# Monocyte genes:
Seurat::FeaturePlot(seu, monocyte_genes, ncol=2)
Seurat::VlnPlot(seu,
                features = monocyte_genes,
                ncol = 2) # cluster 2
DefaultAssay(seu)<-"RNA"
# calculate a score per cell for a group of genes:

# This works with Seurat v4:
# seu <- Seurat::AddModuleScore(seu,
#                       features = list(tcell_genes),
#                       name = "tcell_genes")
seu <- Seurat::AddModuleScore(seu,
                                  features = list(tcell_genes),
                                  name = "tcell_genes")

head(seu@meta.data)
Seurat::FeaturePlot(seu, "tcell_genes1", ncol=2)
Seurat::VlnPlot(seu,
                features = "tcell_genes1",
                ncol = 2) # cluster 0 and 8

# Use UCell for scoring gene signatures in single cells:

# For R version 4.2 or higher, install package:
# BiocManager::install("UCell")
library(UCell)
set.seed(123)
signatures <- list(Immune = c("PTPRC"), 
                   Macrophage = c("CTSB", "C1QB", "LAPTM5",
                                   "TYROBP", "PSAP", "C1QA", "HLA-DRA", "CTSD", "NPC2", "FCER1G"), 
                   Tcell = c("CD3D","CD3E", "CD3G", "CD2"), 
                   Bcell = c("MS4A1", "CD79A", "CD79B", "CD19", "BANK1"),
                   Myeloid_cell = c("CD14", "LYZ", "CSF1R", "FCER1G", "SPI1", "LCK-"))

seu <- UCell::AddModuleScore_UCell(seu, features = signatures, name = NULL,
                                      ncores = 4, assay = "RNA")
head(seu@meta.data)
Seurat::VlnPlot(seu, "Tcell") # UCell score ranges from 0 to 1
Seurat::FeatureScatter(seu, feature1 = "tcell_genes1", feature2 = "Tcell",
                       group.by = "RNA_snn_res.0.3")

# identify cycling cells:
# extract the built-in genes for cell cycling:

s.genes <- Seurat::cc.genes.updated.2019$s.genes
g2m.genes <- Seurat::cc.genes.updated.2019$g2m.genes

# add cycling state to meta.data:
seu <- Seurat::CellCycleScoring(seu,
                                    s.features = s.genes,
                                    g2m.features = g2m.genes,
                                     assay="RNA")
head(seu@meta.data)
Seurat::DimPlot(seu, group.by = "Phase")

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
# Dimensions of the metadata :
dim(colData(ref))
# Head of the metadata columns:
head(colData(ref))
table(ref$label.main)

?SingleR
# already ran:
seu_SingleR <- SingleR::SingleR(test = Seurat::GetAssayData(seu),
                                ref = ref,
                                labels = ref$label.main)
# the output contains annotation scores and assigned labels per cell:
head(seu_SingleR)

# assess cell type annotation
SingleR::plotScoreHeatmap(seu_SingleR)
SingleR::plotDeltaDistribution(seu_SingleR)

# remove low frequency cells from annotation
singleR_labels <- seu_SingleR$labels
t <- table(singleR_labels) # count the total number of cells per cell type
t
other <- names(t)[t < 10]
other
# assign to NA:
singleR_labels[singleR_labels %in% other] <- "none" # or NA as in course website

# add the singleR annotation to the integrated seurat's object meta.data:
seu$SingleR_annot <- singleR_labels

# UMAP with Seurat or dittoSeq:
Seurat::DimPlot(seu, group.by = "SingleR_annot")
dittoSeq::dittoDimPlot(seu, "SingleR_annot", size = 0.7) # can take Seurat, SingleCellExperiment, or SummarizedExperiment object.

# barplot of cell type numbers per sample:
dittoSeq::dittoBarPlot(seu, var = "SingleR_annot", group.by = "orig.ident")

# proportion of cell types per cluster at res. 0.3 (cluster 0 and 8 could be T cells)
dittoSeq::dittoBarPlot(seu, 
                       var = "SingleR_annot", 
                       group.by = "RNA_snn_res.0.3")

# You can export UMAP coordinates and cell clustering if needed:
write.csv(cbind(cell=rownames(seu@reductions$umap@cell.embeddings),
                seu@reductions$umap@cell.embeddings),
          "UMAP_coord.csv", row.names = F, quote = F)
write.csv(cbind(cell=rownames(seu@meta.data),
                res_0.3=seu$RNA_snn_res.0.3,
                singleR=seu$SingleR_annot),
          "cell_clustering_ID.csv", row.names = F, quote = F)

saveRDS(seu, "seu_day2_part2.rds")


# save data set and clear environment:
rm(list = ls())
gc()
.rs.restartR()



