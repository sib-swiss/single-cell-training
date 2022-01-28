sample_info <- read.csv("sample_info_course.csv")
datadirs <- file.path("course_data", "count_matrices", sample_info$SampleName,
                      "outs", "filtered_feature_bc_matrix")
names(datadirs) <- gsub("_", "-", sample_info$SampleName)
datadirs

datadirs <- datadirs[1:3]

library(Seurat)
sce_data <- Seurat::Read10X(data.dir = datadirs)
seu <- Seurat::CreateSeuratObject(counts = sce_data,
                                  project = "pbmmc",
                                  min.cells = 3,
                                  min.features = 100)

Seurat::VlnPlot(seu, features = c("nCount_RNA",
                                  "nFeature_RNA"))

seu <- Seurat::PercentageFeatureSet(seu, 
                                    pattern = "^MT-", 
                                    col.name = "percent.mito")

seu <- Seurat::PercentageFeatureSet(seu, 
                                    pattern = "^RP[SL]",
                                    col.name = "percent.ribo")


seu <- Seurat::PercentageFeatureSet(seu,
                                    pattern = "^HB[^(P)]",
                                    col.name = "percent.globin")

Seurat::VlnPlot(seu, features = c("percent.mito",
                                  "percent.ribo",
                                  "percent.globin"))

Seurat::FeatureScatter(seu, 
                       feature1 = "percent.globin", 
                       feature2 = "percent.ribo")


library(ggplot2)
library(Matrix)

most_expressed_boxplot <- function(object, ngenes = 20){
  
  # matrix of raw counts
  cts <- GetAssayData(seu, assay = "RNA", slot = "counts")
  
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

most_expressed_boxplot(seu, 20)

seu <- subset(seu, subset = nFeature_RNA > 200 & 
                nFeature_RNA < 5000 &
                percent.mito < 8)

VlnPlot(seu, features = c("nFeature_RNA",
                          "nCount_RNA",
                          "percent.mito"))

Seurat::GetAssay(seu)[1:10,1:10]  

seu <- Seurat::NormalizeData(seu,
                             normalization.method = "LogNormalize",
                             scale.factor = 10000)

seu <- Seurat::FindVariableFeatures(seu,
                                    selection.method = "vst",
                                    nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(Seurat::VariableFeatures(seu), 10)
top10

vf_plot <- Seurat::VariableFeaturePlot(seu)
Seurat::LabelPoints(plot = vf_plot,
                    points = top10, repel = TRUE)

seu <- Seurat::ScaleData(seu,
                         features = rownames(seu))

s.genes <- Seurat::cc.genes.updated.2019$s.genes
g2m.genes <- Seurat::cc.genes.updated.2019$g2m.genes

seu <- Seurat::CellCycleScoring(seu,
                                s.features = s.genes,
                                g2m.features = g2m.genes)

seu <- Seurat::RunPCA(seu)

seu <- Seurat::RunUMAP(seu, dims = 1:25)

Seurat::DimPlot(seu, reduction = "umap")

Seurat::DimPlot(seu, group.by = "Phase")

seu.list <- Seurat::SplitObject(seu, split.by = "orig.ident")

for (i in 1:length(seu.list)) {
  seu.list[[i]] <- Seurat::NormalizeData(seu.list[[i]])
  seu.list[[i]] <- Seurat::FindVariableFeatures(seu.list[[i]], selection.method = "vst", nfeatures = 2000,
                                                     verbose = FALSE)
}

seu.anchors <- Seurat::FindIntegrationAnchors(object.list = seu.list, dims = 1:30)

seu.integrated <- Seurat::IntegrateData(anchorset = seu.anchors, dims = 1:30)

Seurat::DefaultAssay(seu.integrated) <- "integrated"

seu.integrated <- Seurat::ScaleData(seu.integrated)
seu.integrated <- Seurat::RunPCA(seu.integrated, npcs = 30)
seu.integrated <- Seurat::RunUMAP(seu.integrated, reduction = "pca", dims = 1:30)

Seurat::DimPlot(seu.integrated, reduction = "umap", group.by = "orig.ident")
Seurat::DimPlot(seu.integrated, reduction = "umap", group.by = "Phase")

# install.packages("matrixStats", repos = "http://cran.us.r-project.org")
# BiocManager::install("celldex")
# BiocManager::install("SingleR")
library(celldex)
hpca.se <- celldex::HumanPrimaryCellAtlasData()
gbm_SingleR <- SingleR::SingleR(test = Seurat::GetAssayData(seu, slot = "data"),
                                  ref = hpca.se,
                                  labels = hpca.se$label.main)

seu.integrated$SingleR_annot <- gbm_SingleR$labels
seu$SingleR_annot <- gbm_SingleR$labels

Seurat::DimPlot(seu.integrated, group.by = "SingleR_annot", label = T, repel = T)
Seurat::DimPlot(seu, group.by = "SingleR_annot", label = T, repel = T)
t <- table(seu$orig.ident, seu$SingleR_annot)
lown <- colnames(t)[colSums(t) < 30]
seu.integrated$SingleR_annot[seu.integrated$SingleR_annot %in% lown] <- "other"
Seurat::DimPlot(seu.integrated, group.by = "SingleR_annot", label = T, repel = T)
t <- table(seu.integrated$orig.ident, seu.integrated$SingleR_annot)
tn <- t/rowSums(t)
barplot(t(tn))

Seurat::DefaultAssay(seu.integrated) <- "RNA"
Seurat::FeaturePlot(seu.integrated, features = "CD79A")

Seurat::FeaturePlot(seu.integrated, features = "CST3")
Seurat::FeaturePlot(seu.integrated, features = "CD3D")
Seurat::FeaturePlot(seu.integrated, features = "HBA1")
Seurat::FeaturePlot(seu.integrated, features = "CD34")
Seurat::FeaturePlot(seu.integrated, features = "MS4A1") # CD20
Seurat::FeaturePlot(seu.integrated, features = "SPN")
Seurat::FeaturePlot(seu.integrated, features = "KIT")
Seurat::DimPlot(seu.integrated, group.by = "Phase")

# BiocManager::install("dittoSeq")
library(dittoSeq)
dittoBarPlot(seu.integrated, var = "SingleR_annot", group.by = "orig.ident")
dittoBarPlot(seu.integrated, var = "Phase", group.by = "SingleR_annot")
dittoDimPlot(seu.integrated, var = "SingleR_annot", size = 0.25)

seu.integrated$pre_bcell <- 0
seu.integrated$pre_bcell[seu.integrated$SingleR_annot == "Pre-B_cell_CD34-"] <- 1
dittoDimPlot(seu.integrated, var = "pre_bcell", size = 0.25)
dittoDimPlot(seu.integrated, var = "CD34", 
             size = 0.5, min.color = "grey", max.color = "darkred")
dittoDimPlot(seu.integrated, var = "Phase", size = 0.25)

Seurat::DefaultAssay(seu.integrated) <- "integrated"
seu.integrated <- Seurat::FindNeighbors(seu.integrated, dims = 1:25)
seu.integrated <- Seurat::FindClusters(seu.integrated, resolution = seq(0.1, 0.8, by=0.1))

dittoDimPlot(seu.integrated, var = "integrated_snn_res.0.6", size = 0.25)
bcells <- subset(seu.integrated, subset = integrated_snn_res.0.6 %in% c(6, 11, 3, 12, 10))
bcells <- Seurat::FindVariableFeatures(bcells)
bcells <- Seurat::ScaleData(bcells)
bcells <- Seurat::RunPCA(bcells, npcs = 30)

bcells <- Seurat::RunUMAP(bcells, reduction = "pca", dims = 1:30)
dittoDimPlot(bcells, var = "integrated_snn_res.0.6")

dittoDotPlot(seu.integrated, vars = c("CD34", "MS4A1",
                                      "HBA1", "CD3D",
                                      "CST3", "CD79A"), 
             group.by = "SingleR_annot",
             min.color = "grey", max.color = "red")

# GMP: Granulocyte-monocyte progenitor 
# MEP: Megakaryocyte (platelets) erythrocyte (red blood cell) progenitor
# CMP: common myeloid progenitor
# BM: bone marrow
# Granulocytes: Basophil neutrophil
# Myeloid: = bone marrow: granulocytes + monocytes + erythrocytes + Mgk
# CD34+: hematopoietic stem cells

# in paper large part is annotated as erythrocytes, but they don't have a nucleus. 

# 

pca <- Embeddings(bcells,  reduction = "pca")

bcells$PC1 <- pca[, 2]
bcells$PC2 <- pca[, 3]

library(Seurat)
library(SingleCellExperiment)
ggplot(as.data.frame(bcells@meta.data), aes(x = PC1, y = PC2, color = integrated_snn_res.0.6)) +
  geom_point(size=2, shape=20) +
  theme_classic() +
  xlab("PC1") + ylab("PC2") + ggtitle("PC biplot")

