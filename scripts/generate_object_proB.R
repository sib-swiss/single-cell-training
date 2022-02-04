sample_info <- read.csv("course_data/sample_info_course.csv")
datadirs <- file.path("course_data", "count_matrices", sample_info$SampleName,
                      "outs", "filtered_feature_bc_matrix")
names(datadirs) <- gsub("_", "-", sample_info$SampleName)
datadirs

library(Seurat)
sce_data <- Seurat::Read10X(data.dir = datadirs)
all <- Seurat::CreateSeuratObject(counts = sce_data,
                                  project = "pbmmc",
                                  min.cells = 3,
                                  min.features = 100)
DefaultAssay(all) <- "RNA"
ref <- celldex::HumanPrimaryCellAtlasData()

all_SingleR <- SingleR::SingleR(test = Seurat::GetAssayData(all, slot = "data"),
                                    ref = ref,
                                    labels = ref$label.main)

all$SingleR_annot <- all_SingleR$labels

proB <- subset(all, subset = SingleR_annot == "Pro-B_cell_CD34+")
proB$orig.ident <- as.character(proB$orig.ident)
proB$type <- substring(proB$orig.ident, 1, nchar(proB$orig.ident) - 2)

# subset malignant cells to 2000
mal_proB <- colnames(proB)[proB$type == "ETV6-RUNX1"]
set.seed(1)
mal_proB_sub <- sample(mal_proB, 2000)

proB <- subset(proB, cells = c(colnames(proB)[proB$type == "PBMMC"],
                                mal_proB_sub))
DefaultAssay(proB) <- "RNA"

# proB <- Seurat::NormalizeData(proB)
# proB <- Seurat::FindVariableFeatures(proB)
# proB <- Seurat::ScaleData(proB, assay = "RNA")
# 
# s.genes <- Seurat::cc.genes.updated.2019$s.genes
# g2m.genes <- Seurat::cc.genes.updated.2019$g2m.genes
# 
# proB <- Seurat::CellCycleScoring(proB,
#                                      s.features = s.genes,
#                                      g2m.features = g2m.genes)
# 
# proB <- subset(proB, subset = Phase == "G1")
# DefaultAssay(proB) <- "RNA"

proB <- Seurat::NormalizeData(proB)
proB <- Seurat::FindVariableFeatures(proB)
proB <- Seurat::ScaleData(proB, assay = "RNA")
proB <- Seurat::RunPCA(proB, npcs = 30, assay = "RNA")
proB <- Seurat::RunUMAP(proB, reduction = "pca", dims = 1:30, assay = "RNA")
DimPlot(proB, group.by = "orig.ident")
DimPlot(proB, group.by = "type")
FeaturePlot(proB, features = "CD52")

saveRDS(proB, "proB.rds")
