sample_info <- read.csv("course_data/sample_info_course.csv")
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

saveRDS(seu, "seu_day1.rds")

# s.genes <- Seurat::cc.genes.updated.2019$s.genes
# g2m.genes <- Seurat::cc.genes.updated.2019$g2m.genes
# 
# seu <- Seurat::CellCycleScoring(seu,
#                                 s.features = s.genes,
#                                 g2m.features = g2m.genes)

seu <- readRDS("seu_day1.rds")

seu <- Seurat::RunPCA(seu)
Seurat::DimPlot(seu, reduction = "pca")
Seurat::FeaturePlot(seu, reduction = "pca", features = "percent.mito")
Seurat::FeaturePlot(seu, reduction = "pca", features = "percent.ribo")
Seurat::FeaturePlot(seu, reduction = "pca", features = "percent.globin")
Seurat::FeaturePlot(seu, reduction = "pca", features = "HBA1")

Seurat::DimHeatmap(seu, dims = 1:12, cells = 500, balanced = TRUE)


seu <- Seurat::RunUMAP(seu, dims = 1:25)

Seurat::DimPlot(seu, reduction = "umap")

Seurat::FeaturePlot(seu, features = c("HBA1", "percent.globin", "IGKC", "percent.mito"))

seu <- Seurat::RunUMAP(seu, dims = 1:25, n.neighbors = 5)
Seurat::DimPlot(seu, reduction = "umap")

seu <- Seurat::RunUMAP(seu, dims = 1:5)
Seurat::DimPlot(seu, reduction = "umap")

seu <- Seurat::RunUMAP(seu, dims = 1:25)

seu_list <- Seurat::SplitObject(seu, split.by = "orig.ident")

for (i in 1:length(seu_list)) {
  seu_list[[i]] <- Seurat::NormalizeData(seu_list[[i]])
  seu_list[[i]] <- Seurat::FindVariableFeatures(seu_list[[i]], selection.method = "vst", nfeatures = 2000,
                                                verbose = FALSE)
}

seu_anchors <- Seurat::FindIntegrationAnchors(object.list = seu_list, dims = 1:30)

seu_int <- Seurat::IntegrateData(anchorset = seu_anchors, dims = 1:30)

Seurat::DefaultAssay(seu_int) <- "integrated"

seu_int <- Seurat::ScaleData(seu_int)
seu_int <- Seurat::RunPCA(seu_int, npcs = 30)
seu_int <- Seurat::RunUMAP(seu_int, reduction = "pca", dims = 1:30)

Seurat::DimPlot(seu_int, reduction = "umap")


seu_int <- Seurat::FindNeighbors(seu_int, dims = 1:25)

seu_int <- Seurat::FindClusters(seu_int, resolution = seq(0.1, 0.8, by=0.1))

head(seu_int@meta.data)

library(clustree)
clustree::clustree(seu_int@meta.data[,grep("integrated_snn_res", colnames(seu_int@meta.data))],
                   prefix = "integrated_snn_res.")

Seurat::DimPlot(seu_int, group.by = "integrated_snn_res.0.3")

seu_int <- Seurat::SetIdent(seu_int, value = seu_int$integrated_snn_res.0.3)
DefaultAssay(seu_int) <- "RNA"

Seurat::FeaturePlot(seu_int, "HBA1")

immune_genes <- c("IL7R", "LTB", "TRAC", "CD3D")




Seurat::FeaturePlot(seu_int, immune_genes, ncol=2)

Seurat::VlnPlot(seu_int,
                features = immune_genes,
                ncol = 2)

monocyte_genes <- c("CD14", "CST3", "CD68", "CTSS")
Seurat::FeaturePlot(seu_int, monocyte_genes, ncol=2)

Seurat::VlnPlot(seu_int,
                features = monocyte_genes,
                ncol = 2)

seu_int <- Seurat::AddModuleScore(seu_int,
                                  features = list(immune_genes),
                                  name = "immune_genes")

Seurat::FeaturePlot(seu_int, "immune_genes1")

Seurat::VlnPlot(seu_int,
                features = "immune_genes1")


hpca.se <- celldex::HumanPrimaryCellAtlasData()
class(hpca.se)
table(hpca.se$label.main)

ref <- celldex::BlueprintEncodeData()
ref <- celldex::NovershternHematopoieticData()

seu_int_SingleR <- SingleR::SingleR(test = Seurat::GetAssayData(seu_int, slot = "data"),
                                    ref = ref,
                                    labels = ref$label.main)

library(SingleR)
plotScoreHeatmap(seu_int_SingleR)
plotDeltaDistribution(seu_int_SingleR)

singleR_labels <- seu_int_SingleR$labels
t <- table(singleR_labels)
other <- names(t)[t < 10]
singleR_labels[singleR_labels %in% other] <- NA

seu_int$SingleR_annot <- singleR_labels

library(dittoSeq)
dittoDimPlot(seu_int, "SingleR_annot", size = 0.7)
dittoDimPlot(seu_int, "integrated_snn_res.0.3", size = 0.7)
dittoBarPlot(seu_int, var = "SingleR_annot", group.by = "orig.ident")
dittoSeq::dittoBarPlot(seu_int, 
                       var = "SingleR_annot", 
                       group.by = "integrated_snn_res.0.3")


de_genes <- Seurat::FindAllMarkers(seu_int,  
                                   min.pct = 0.25,
                                   only.pos = TRUE)

de_genes <- subset(de_genes, de_genes$p_val_adj < 0.05 & de_genes$avg_log2FC > 0)
View(de_genes)

tcell_genes <- c("IL7R", "LTB", "TRAC", "CD3D")

de_genes[de_genes$gene %in% tcell_genes,]

seu_int <- Seurat::SetIdent(seu_int, value = "SingleR_annot")

deg_cd8_cd4 <- Seurat::FindMarkers(seu_int,
                                   ident.1 = "CD8+ T cells",
                                   ident.2 = "CD4+ T cells",
                                   group.by = seu_int$SingleR_annot,
                                   test.use = "wilcox")

deg_cd8_cd4[c("CD4", "CD8A", "CD8B"),]

Seurat::VlnPlot(seu_int, 
                features = c("CD4", "CD8A", "CD8B"),
                idents = c("CD8+ T cells", "CD4+ T cells"))

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

all_list <- Seurat::SplitObject(all, split.by = "orig.ident")

for (i in 1:length(all_list)) {
  all_list[[i]] <- Seurat::NormalizeData(all_list[[i]])
  all_list[[i]] <- Seurat::FindVariableFeatures(all_list[[i]], selection.method = "vst", nfeatures = 2000,
                                                verbose = FALSE)
}

all_anchors <- Seurat::FindIntegrationAnchors(object.list = all_list, dims = 1:30)

all_int <- Seurat::IntegrateData(anchorset = all_anchors, dims = 1:30)

Seurat::DefaultAssay(all_int) <- "integrated"

all_int <- Seurat::ScaleData(all_int)
all_int <- Seurat::RunPCA(all_int, npcs = 30)
all_int <- Seurat::RunUMAP(all_int, reduction = "pca", dims = 1:30)



DimPlot(all_int)

DefaultAssay(all_int) <- "RNA"
ref <- celldex::HumanPrimaryCellAtlasData()

all_int_SingleR <- SingleR::SingleR(test = Seurat::GetAssayData(all_int, slot = "data"),
                                    ref = ref,
                                    labels = ref$label.main)

singleR_labels <- all_int_SingleR$labels
t <- table(singleR_labels)
other <- names(t)[t < 10]
singleR_labels[singleR_labels %in% other] <- NA

all_int$SingleR_annot <- singleR_labels
DimPlot(all_int, group.by = "SingleR_annot")
dittoDimPlot(all_int, var = "SingleR_annot", size = 0.3)

dittoBarPlot(all_int, var = "SingleR_annot", group.by = "orig.ident")

DefaultAssay(all_int) <- "integrated"
all_int <- Seurat::FindNeighbors(all_int, dims = 1:25)


all_int <- Seurat::FindClusters(all_int, resolution = seq(0.1, 0.6, by=0.1))

dittoDimPlot(all_int, var = "integrated_snn_res.0.2")

dittoBarPlot(all_int, var = "integrated_snn_res.0.2", group.by = "orig.ident")

all_prob <- subset(all_int, subset = SingleR_annot == "Pro-B_cell_CD34+")
all_prob$type <- substring(all_prob$orig.ident, 1, nchar(all_prob$orig.ident) - 2)

DefaultAssay(all_prob) <- "RNA"

s.genes <- Seurat::cc.genes.updated.2019$s.genes
g2m.genes <- Seurat::cc.genes.updated.2019$g2m.genes

all_prob <- Seurat::CellCycleScoring(all_prob,
                                     s.features = s.genes,
                                     g2m.features = g2m.genes)

all_prob <- subset(all_prob, subset = Phase == "G1")
DefaultAssay(all_prob) <- "RNA"

all_prob <- Seurat::ScaleData(all_prob, assay = "RNA")
all_prob <- Seurat::RunPCA(all_prob, npcs = 30, assay = "RNA")
all_prob <- Seurat::RunUMAP(all_prob, reduction = "pca", dims = 1:30, assay = "RNA")
DimPlot(all_prob, group.by = "orig.ident")
DimPlot(all_prob, group.by = "type")
FeaturePlot(all_prob, features = "DNTT")

saveRDS(all_prob, "all_prob.rds")

counts <- Seurat::GetAssayData(all_prob, slot = "counts")
counts <- counts[rowSums(counts) != 0,]

dge <- edgeR::DGEList(counts = counts)
dge <- edgeR::calcNormFactors(dge)  

design <- model.matrix(~ 0 + type, data = all_prob@meta.data)
colnames(design)<-make.names(c("ETV6-RUNX1", "PBMMC"))

contrast.mat <- limma::makeContrasts(PBMMC - ETV6.RUNX1,
                                     levels = design)

vm <- limma::voom(dge, design = design, plot = TRUE)
fit <- limma::lmFit(vm, design = design)
fit.contrasts <- limma::contrasts.fit(fit, contrast.mat)
fit.contrasts <- limma::eBayes(fit.contrasts)

limma::topTable(fit.contrasts, number = 10, sort.by = "logFC")

Seurat::VlnPlot(all_prob, "CD52", split.by = "orig.ident")
Seurat::VlnPlot(all_prob, "IGLL1", split.by = "orig.ident")
Seurat::VlnPlot(all_prob, "CD74", split.by = "orig.ident")
Seurat::VlnPlot(all_prob, "DNTT", split.by = "orig.ident")

DefaultAssay(all_prob) <- "RNA"
tum_vs_norm <- FindMarkers(all_prob, 
                           ident.1 = "ETV6-RUNX1", 
                           ident.2 = "PBMMC", 
                           group.by = "type")

BiocManager::install("org.Hs.eg.db", update = FALSE)
library(org.Hs.eg.db)
AnnotationDbi::keytypes(org.Hs.eg.db)

tum_down <- subset(tum_vs_norm,
                   tum_vs_norm$avg_log2FC < -1 &
                     tum_vs_norm$p_val_adj < 0.05)
tum_down_genes <- rownames(tum_down)

tum_vs_norm_go <- clusterProfiler::enrichGO(tum_down_genes,
                                            "org.Hs.eg.db",
                                            keyType = "SYMBOL",
                                            ont = "BP",
                                            minGSSize = 50)

enr_go <- clusterProfiler::simplify(tum_vs_norm_go)
View(enr_go@result)


tum_up <- subset(tum_vs_norm,
                 tum_vs_norm$avg_log2FC > 1 &
                   tum_vs_norm$p_val_adj < 0.05)
tum_up_genes <- rownames(tum_up)

tum_vs_norm_go <- clusterProfiler::enrichGO(tum_up_genes,
                                            "org.Hs.eg.db",
                                            keyType = "SYMBOL",
                                            ont = "BP",
                                            minGSSize = 50)

enr_go <- clusterProfiler::simplify(tum_vs_norm_go)
View(enr_go@result)

install.packages("msigdbr")

gmt <- msigdbr::msigdbr(species = "human", category = "H")

tum_vs_norm_enrich <- clusterProfiler::enricher(gene = tum_down_genes,
                                                universe = rownames(all_prob),
                                                pAdjustMethod = "BH",
                                                pvalueCutoff  = 0.05,
                                                qvalueCutoff  = 0.05,
                                                TERM2GENE = gmt[,c("gs_name", "gene_symbol")])

View(tum_vs_norm_enrich@result)
VlnPlot(all_prob, features = c("TERF2", "IGLL1"), group.by = "type")

library(Seurat)
DimPlot(seu_int, group.by = "integrated_snn_res.0.3", repel = TRUE, label = TRUE)
seuB <- subset(seu_int, subset = integrated_snn_res.0.3 %in% c(11, 3, 7)) 
# seuB <- FindVariableFeatures(seuB, nfeatures = 200)
# seuB <- Seurat::ScaleData(seuB)
# seuB <- Seurat::RunPCA(seuB, npcs = 30)
# seuB <- Seurat::RunUMAP(seuB, reduction = "pca", dims = 1:30, n.neighbors = 10)
DimPlot(seuB)

library(monocle3)
feature_names <- as.data.frame(rownames(seu_int))
rownames(feature_names) <- rownames(seu_int)
colnames(feature_names) <- "gene_short_name"
seu_int_monocl <- monocle3::new_cell_data_set(seu_int@assays$RNA@counts,
                                              cell_metadata = seu_int@meta.data,
                                              gene_metadata = feature_names)

seu_int_monocl <- monocle3::preprocess_cds(seu_int_monocl)
monocle3::plot_pc_variance_explained(seu_int_monocl)


seu_int_monocl <- monocle3::reduce_dimension(seu_int_monocl, reduction_method = "UMAP")
monocle3::plot_cells(seu_int_monocl, 
                     color_cells_by = "integrated_snn_res.0.3", 
                     cell_size = 1, 
                     show_trajectory_graph = FALSE)
monocle3::plot_cells(seu_int_monocl, genes = "CD34", cell_size = 1) # to plot expression level of a gene
monocle3::plot_cells(seu_int_monocl, genes = "MS4A1", cell_size = 1) 
monocle3::plot_cells(seu_int_monocl, genes = "CD79A", 
                     show_trajectory_graph = FALSE, 
                     cell_size = 1)
monocle3::plot_cells(seu_int_monocl, genes = "CD8A", show_trajectory_graph = FALSE, cell_size = 1)

monocle3::plot_cells(seu_int_monocl, genes = c("CD79A", "CD34"),
                     show_trajectory_graph = FALSE, 
                     cell_size = 0.7, group_label_size = 5)


monocle3::plot_cells(seu_int_monocl, color_cells_by = "SingleR_annot", 
                     cell_size = 1, show_trajectory_graph = FALSE,
                     label_cell_groups = F)

p1 <- monocle3::plot_cells(seu_int_monocl, label_cell_groups = F)
p2 <- monocle3::plot_cells(seu_int_monocl, color_cells_by = "integrated_snn_res.0.3", 
                           label_cell_groups = F)
cowplot::plot_grid(p1, p2, ncol = 2)

seu_int_monocl <- monocle3::cluster_cells(seu_int_monocl, resolution=0.00025)

monocle3::plot_cells(seu_int_monocl, label_cell_groups = F)

seu_int_monocl <- monocle3::learn_graph(seu_int_monocl)
monocle3::plot_cells(seu_int_monocl)

seu_int_monocl<-monocle3::order_cells(seu_int_monocl)

monocle3::plot_cells(seu_int_monocl,
                     color_cells_by = "pseudotime",
                     label_cell_groups=F,
                     label_leaves=F,
                     label_branch_points=FALSE,
                     graph_label_size=1.5, cell_size = 1)

plot_genes_in_pseudotime(subset(seu_int_monocl, 
                                rowData(seu_int_monocl)$gene_short_name=="HMGB1"))

monocle3::plot_cells(seu_int_monocl, genes = "STMN1", show_trajectory_graph = FALSE, cell_size =  1) 

seuB <- choose_cells(seu_int_monocl)
plot_cells(seuB, show_trajectory_graph = FALSE, cell_size = 1)
pr_test <- graph_test(seuB, cores=4, neighbor_graph = "principal_graph")
pr_test <- pr_test[order(pr_test$morans_test_statistic, decreasing = TRUE),]

goi <- c("CD34", "MS4A1", "IGLL1", "IGLL5", 
         "MKI67", "CKS2")
plot_cells(seuB, label_cell_groups=FALSE, genes = goi,
           show_trajectory_graph=FALSE, cell_size = 1)

seuB@colData$monocle_cluster <- clusters(seuB)

plot_genes_in_pseudotime(subset(seuB, 
                                rowData(seu_int_monocl)$gene_short_name %in% goi),
                         min_expr=0.5, color_cells_by = "monocle_cluster")

plot_genes_in_pseudotime(subset(seuB, 
                                rowData(seuB)$gene_short_name %in% goi),
                         min_expr=0.5, color_cells_by = "monocle_cluster")


plot_genes_in_pseudotime(subset(seu_int_monocl, 
                                rowData(seu_int_monocl)$gene_short_name=="LTB"),
                         min_expr=0.5, color_cells_by = "integrated_snn_res.0.3")

BiocManager::install("org.Hs.eg.db")
bcell_pseudot <- clusterProfiler::enrichGO(pr_test$gene_short_name[1:200],
                                           "org.Hs.eg.db",
                                           keyType = "SYMBOL",
                                           ont = "BP",
                                           minGSSize = 50)
bcell_pseudot <- clusterProfiler::simplify(bcell_pseudot)

saveRDS(seu_int_monocl, "seu_int_monocl.rds")
