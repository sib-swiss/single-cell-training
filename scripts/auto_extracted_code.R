
####################################
#                                  #
#       Auto-extracted code        #
# DO NOT MANUALLY CHANGE THIS FILE #
#                                  #
####################################


## Code found in: index.md


## Code found in: precourse.md


## Code found in: course_schedule.md


## Code found in: day1/general_introduction.md


## Code found in: day1/setup.md


## Code found in: day1/introduction_scrnaseq.md


## Code found in: day1/analysis_tools_qc.md


## Code found in: day1/normalization_scaling.md
Seurat::GetAssayData(seu)[1:10,1:10]  

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

seu <- Seurat::SCTransform(seu)

DefaultAssay(seu) <- "RNA"

saveRDS(seu, "seu_day1.rds")

rm(list = ls())
gc()
.rs.restartR()

## Code found in: day2/dimensionality_reduction.md
seu <- readRDS("seu_day1.rds")

library(Seurat)

seu <- Seurat::RunPCA(seu)

Seurat::DimPlot(seu, reduction = "pca")

Seurat::FeaturePlot(seu, reduction = "pca", features = "percent.globin")

Seurat::FeaturePlot(seu, reduction = "pca", features = "HBA1")

Seurat::DimHeatmap(seu, dims = 1:12, cells = 500, balanced = TRUE)

Seurat::ElbowPlot(seu, ndims = 40)

seu <- Seurat::RunUMAP(seu, dims = 1:25)

Seurat::DimPlot(seu, reduction = "umap")

Seurat::FeaturePlot(seu, features = c("HBA1", "percent.globin", "IGKC", "percent.mito"))

seu <- Seurat::RunUMAP(seu, dims = 1:25, n.neighbors = 5)
Seurat::DimPlot(seu, reduction = "umap")

seu <- Seurat::RunUMAP(seu, dims = 1:5)
Seurat::DimPlot(seu, reduction = "umap")

seu <- Seurat::RunUMAP(seu, dims = 1:50) 
Seurat::DimPlot(seu, reduction = "umap")

seu <- Seurat::RunUMAP(seu, dims = 1:25)

## Code found in: day2/integration.md
Seurat::DimPlot(seu, reduction = "umap")

seu_list <- Seurat::SplitObject(seu, split.by = "orig.ident")

for (i in 1:length(seu_list)) {
seu_list[[i]] <- Seurat::NormalizeData(seu_list[[i]])
seu_list[[i]] <- Seurat::FindVariableFeatures(seu_list[[i]], selection.method = "vst", nfeatures = 2000,
verbose = FALSE)
}

seu_anchors <- Seurat::FindIntegrationAnchors(object.list = seu_list, dims = 1:25)

seu_int <- Seurat::IntegrateData(anchorset = seu_anchors, dims = 1:25)

Seurat::DefaultAssay(seu_int) <- "integrated"

seu_int <- Seurat::ScaleData(seu_int)
seu_int <- Seurat::RunPCA(seu_int, npcs = 30)
seu_int <- Seurat::RunUMAP(seu_int, reduction = "pca", dims = 1:25)

Seurat::DimPlot(seu_int, reduction = "umap")

saveRDS(seu_int, "seu_int_day2_part1.rds")

rm(list = ls())
gc()
.rs.restartR()

## Code found in: day2/clustering.md
seu_int <- readRDS("seu_int_day2_part1.rds")

seu_int <- Seurat::FindNeighbors(seu_int, dims = 1:25)

seu_int <- Seurat::FindClusters(seu_int, resolution = seq(0.1, 0.8, by=0.1))

head(seu_int@meta.data)

library(clustree)
clustree::clustree(seu_int@meta.data[,grep("integrated_snn_res", colnames(seu_int@meta.data))],
prefix = "integrated_snn_res.")

Seurat::DimPlot(seu_int, group.by = "integrated_snn_res.0.1")

Seurat::DimPlot(seu_int, group.by = "integrated_snn_res.0.3")

## Code found in: day2/cell_annotation.md
library(celldex)
library(SingleR)

seu_int <- Seurat::SetIdent(seu_int, value = seu_int$integrated_snn_res.0.3)

DefaultAssay(seu_int) <- "RNA"

Seurat::FeaturePlot(seu_int, "HBA1")

tcell_genes <- c("IL7R", "LTB", "TRAC", "CD3D")
monocyte_genes <- c("CD14", "CST3", "CD68", "CTSS")

Seurat::FeaturePlot(seu_int, tcell_genes, ncol=2)

Seurat::VlnPlot(seu_int,
features = tcell_genes,
ncol = 2)

Seurat::FeaturePlot(seu_int, monocyte_genes, ncol=2)

Seurat::VlnPlot(seu_int,
features = monocyte_genes,
ncol = 2)

seu_int <- Seurat::AddModuleScore(seu_int,
features = list(tcell_genes),
name = "tcell_genes")

Seurat::FeaturePlot(seu_int, "tcell_genes1")

Seurat::VlnPlot(seu_int,
"tcell_genes1")

s.genes <- Seurat::cc.genes.updated.2019$s.genes
g2m.genes <- Seurat::cc.genes.updated.2019$g2m.genes

seu_int <- Seurat::CellCycleScoring(seu_int,
s.features = s.genes,
g2m.features = g2m.genes)

Seurat::DimPlot(seu_int, group.by = "Phase")

ref <- celldex::NovershternHematopoieticData()
class(ref)
table(ref$label.main)

seu_int_SingleR <- SingleR::SingleR(test = Seurat::GetAssayData(seu_int, slot = "data"),
ref = ref,
labels = ref$label.main)

head(seu_int_SingleR)

SingleR::plotScoreHeatmap(seu_int_SingleR)

SingleR::plotDeltaDistribution(seu_int_SingleR)

singleR_labels <- seu_int_SingleR$labels
t <- table(singleR_labels)
other <- names(t)[t < 10]
singleR_labels[singleR_labels %in% other] <- "none"

seu_int$SingleR_annot <- singleR_labels

dittoSeq::dittoDimPlot(seu_int, "SingleR_annot", size = 0.7)

dittoSeq::dittoBarPlot(seu_int, var = "SingleR_annot", group.by = "orig.ident")

dittoSeq::dittoBarPlot(seu_int, 
var = "SingleR_annot", 
group.by = "integrated_snn_res.0.3")

saveRDS(seu_int, "seu_int_day2_part2.rds")

rm(list = ls())
gc()
.rs.restartR()

## Code found in: day3/differential_gene_expression.md
seu_int <- readRDS("seu_int_day2_part2.rds")

library(Seurat)
library(edgeR) # BiocManager::install("edgeR")
library(limma)

de_genes <- Seurat::FindAllMarkers(seu_int,  min.pct = 0.25,
only.pos = TRUE)

de_genes <- subset(de_genes, de_genes$p_val_adj<0.05)
write.csv(de_genes, "de_genes_FindAllMarkers.csv", row.names = F, quote = F)


library(dplyr)
top_specific_markers <- de_genes %>%
group_by(cluster) %>%
top_n(3, avg_log2FC)

dittoSeq::dittoDotPlot(seu_int, vars = unique(top_specific_markers$gene), 
group.by = "integrated_snn_res.0.3")

tcell_genes <- c("IL7R", "LTB", "TRAC", "CD3D")

de_genes[de_genes$gene %in% tcell_genes,]

seu_int <- Seurat::SetIdent(seu_int, value = "SingleR_annot")

deg_cd8_cd4 <- Seurat::FindMarkers(seu_int,
ident.1 = "CD8+ T cells",
ident.2 = "CD4+ T cells",
group.by = seu_int$SingleR_annot,
test.use = "wilcox")
deg_cd8_cd4 <- subset(deg_cd8_cd4, deg_cd8_cd4$p_val_adj<0.05)

View(deg_cd8_cd4)

deg_cd8_cd4[c("CD4", "CD8A", "CD8B"),]

Seurat::VlnPlot(seu_int, 
features = c("CD4", "CD8A", "CD8B"),
idents = c("CD8+ T cells", "CD4+ T cells"))

proB <- readRDS("course_data/proB.rds")

DimPlot(proB, group.by = "orig.ident")

table(proB@meta.data$type)
# ETV6-RUNX1      PBMMC 
#      2000       1021

head(proB@meta.data)


Seurat::DefaultAssay(proB) <- "RNA"
Seurat::Idents(proB) <- proB$orig.ident

Seurat::DimPlot(proB)

counts <- Seurat::GetAssayData(proB, slot = "counts")
counts <- counts[rowSums(counts) != 0,]
dim(counts)

dge <- edgeR::DGEList(counts = counts)
dge <- edgeR::calcNormFactors(dge)  

proB$patient.id<-gsub("ETV6-RUNX1", "ETV6_RUNX1", proB$orig.ident)
proB$patient.id<-sapply(strsplit(proB$patient.id, "-"), '[', 2)

design <- model.matrix(~ 0 + type + patient.id , 
data = proB@meta.data)

head(design)

# change column names to more simple group names: 
colnames(design)[c(1:2)] <- make.names(c("ETV6-RUNX1", "PBMMC"))


contrast.mat <- limma::makeContrasts(ETV6.RUNX1 - PBMMC,
levels = design)

vm <- limma::voom(dge, design = design, plot = TRUE)
fit <- limma::lmFit(vm, design = design)
fit.contrasts <- limma::contrasts.fit(fit, contrast.mat)
fit.contrasts <- limma::eBayes(fit.contrasts)

limma::topTable(fit.contrasts, number = 10, sort.by = "P")

limma_de <- limma::topTable(fit.contrasts, number = Inf, sort.by = "P")
length(which(limma_de$adj.P.Val<0.05))


Seurat::VlnPlot(proB, "S100A9", split.by = "type")
Seurat::VlnPlot(proB, "SOCS2", split.by = "type")

tum_vs_norm <- Seurat::FindMarkers(proB, 
ident.1 = "ETV6-RUNX1", 
ident.2 = "PBMMC", 
group.by = "type")
tum_vs_norm <- subset(tum_vs_norm, tum_vs_norm$p_val_adj<0.05)

dim(tum_vs_norm) 

merge_limma_FindMarkers <- merge(tum_vs_norm, 
limma_de, 
by="row.names",
all.x=T)
merge_limma_FindMarkers <- subset(merge_limma_FindMarkers,
merge_limma_FindMarkers$adj.P.Val<0.00001)

par(mar=c(4,4,4,4))
plot(merge_limma_FindMarkers$avg_log2FC,
merge_limma_FindMarkers$logFC,
xlab="log2FC Wilcoxon", ylab="log2FC limma",
pch=15, cex=0.5)
abline(a=0, b=1, col="red")

## Code found in: day3/enrichment_analysis.md
library(clusterProfiler)
library(enrichplot)

BiocManager::install("org.Hs.eg.db", update = FALSE)
library(org.Hs.eg.db)
AnnotationDbi::keytypes(org.Hs.eg.db)

tum_down <- subset(tum_vs_norm,
tum_vs_norm$avg_log2FC < -1 &
tum_vs_norm$p_val_adj < 0.05)
tum_down_genes <- rownames(tum_down)

?enrichGO
tum_vs_norm_go <- clusterProfiler::enrichGO(tum_down_genes,
"org.Hs.eg.db",
keyType = "SYMBOL",
ont = "BP",
minGSSize = 50)

View(tum_vs_norm_go@result)

enr_go <- clusterProfiler::simplify(tum_vs_norm_go)
View(enr_go@result)

enrichplot::emapplot(enrichplot::pairwise_termsim(enr_go),
showCategory = 30, cex_label_category = 0.5)

gmt <- msigdbr::msigdbr(species = "human", category = "H")

tum_vs_norm_enrich <- clusterProfiler::enricher(gene = tum_down_genes,
universe = rownames(proB),
pAdjustMethod = "BH",
pvalueCutoff  = 0.05,
qvalueCutoff  = 0.05,
TERM2GENE = gmt[,c("gs_name", "gene_symbol")])

View(tum_vs_norm_enrich@result[which(tum_vs_norm_enrich@result$p.adjust<0.05),])

rm(list = ls())
gc()
.rs.restartR()

## Code found in: day3/advanced_analyses.md


## Code found in: day3/trajectory_analysis.md
BiocManager::install("scater")

library(SingleCellExperiment)
library(scater)
library(slingshot)
library(ggplot2)
library(ggbeeswarm)

deng_SCE <- readRDS("course_data/deng-reads.rds")

deng_SCE$cell_type2 <- factor(deng_SCE$cell_type2,
levels = c("zy",
"early2cell",
"mid2cell",
"late2cell",
"4cell",
"8cell",
"16cell",
"earlyblast",
"midblast",
"lateblast"))

deng_SCE <- scater::runPCA(deng_SCE, ncomponents = 50)

pca <- SingleCellExperiment::reducedDim(deng_SCE, "PCA")

head(pca)

deng_SCE$PC1 <- pca[, 1]
deng_SCE$PC2 <- pca[, 2]

ggplot(as.data.frame(colData(deng_SCE)), aes(x = PC1, y = PC2, color = cell_type2)) +
geom_point(size=2, shape=20) +
theme_classic() +
xlab("PC1") + ylab("PC2") + ggtitle("PC biplot")

deng_SCE$pseudotime_PC1 <- rank(deng_SCE$PC1)  # rank cells by their PC1 score

ggplot(as.data.frame(colData(deng_SCE)), aes(x = pseudotime_PC1, y = cell_type2,
colour = cell_type2)) +
ggbeeswarm::geom_quasirandom(groupOnX = FALSE) +
theme_classic() +
xlab("PC1") + ylab("Timepoint") +
ggtitle("Cells ordered by first principal component")

sce <- slingshot::slingshot(deng_SCE, reducedDim = 'PCA')

PCAplot_slingshot <- function(sce, draw_lines = TRUE, variable = NULL, legend = FALSE, ...){
# set palette for factorial variables
palf <- colorRampPalette(RColorBrewer::brewer.pal(8, "Set2"))
# set palette for numeric variables
paln <- colorRampPalette(RColorBrewer::brewer.pal(9, "Blues"))
# extract pca from SingleCellExperiment object
pca <- SingleCellExperiment::reducedDims(sce)$PCA

if(is.null(variable)){
col <- "black"
}
if(is.character(variable)){
variable <- as.factor(variable)
}
if(is.factor(variable)){
colpal <- palf(length(levels(variable)))
colors <- colpal[variable]
}
if(is.numeric(variable)){
colpal <- paln(50)
colors <- colpal[cut(variable,breaks=50)]
}

# draw the plot
plot(pca, bg = colors, pch = 21)
# draw lines
if(draw_lines){
lines(slingshot::SlingshotDataSet(sce), lwd = 2, ... )
}
# add legend
if(legend & is.factor(variable)){
legend("bottomright", pt.bg = colpal,legend = levels(variable),pch=21)

}
}

PCAplot_slingshot(sce, variable = sce$slingPseudotime_1, draw_lines = TRUE)

ggplot(as.data.frame(colData(deng_SCE)), aes(x = sce$slingPseudotime_1,
y = cell_type2,
colour = cell_type2)) +
ggbeeswarm::geom_quasirandom(groupOnX = FALSE) +
theme_classic() +
xlab("Slingshot pseudotime") + ylab("Timepoint") +
ggtitle("Cells ordered by Slingshot pseudotime")

gcdata <- Seurat::CreateSeuratObject(counts = SingleCellExperiment::counts(deng_SCE),
project = "slingshot")

gcdata <- Seurat::NormalizeData(object = gcdata,
normalization.method = "LogNormalize",
scale.factor = 10000)

gcdata <- Seurat::FindVariableFeatures(object = gcdata,
mean.function = ExpMean,
dispersion.function = LogVMR)

gcdata <- Seurat::ScaleData(object = gcdata,
do.center = T,
do.scale = F)

gcdata <- Seurat::RunPCA(object = gcdata,
pc.genes = gcdata@var.genes)

gcdata <- Seurat::FindNeighbors(gcdata,
reduction = "pca",
dims = 1:5)

# clustering with resolution of 0.6
gcdata <- Seurat::FindClusters(object = gcdata,
resolution = 0.6)

deng_SCE$Seurat_clusters <- as.character(Idents(gcdata))  # go from factor to character

sce <- slingshot::slingshot(deng_SCE,
clusterLabels = 'Seurat_clusters',
reducedDim = 'PCA',
start.clus = "2")

SlingshotDataSet(sce)

PCAplot_slingshot(sce, variable = sce$slingPseudotime_2)

ggplot(data.frame(cell_type2 = deng_SCE$cell_type2,
slingPseudotime_1 = sce$slingPseudotime_1),
aes(x = slingPseudotime_1, y = cell_type2,
colour = cell_type2)) +
ggbeeswarm::geom_quasirandom(groupOnX = FALSE) +
theme_classic() +
xlab("Slingshot pseudotime") + ylab("Timepoint") +
ggtitle("Cells ordered by Slingshot pseudotime")

ggplot(data.frame(cell_type2 = deng_SCE$cell_type2,
slingPseudotime_2 = sce$slingPseudotime_2),
aes(x = slingPseudotime_2, y = cell_type2,
colour = cell_type2)) +
ggbeeswarm::geom_quasirandom(groupOnX = FALSE) +
theme_classic() +
xlab("Slingshot pseudotime") + ylab("Timepoint") +
ggtitle("Cells ordered by Slingshot pseudotime")

PCAplot_slingshot(sce,
variable = deng_SCE$Seurat_clusters,
type = 'lineages',
col = 'black',
legend = TRUE)

PCAplot_slingshot(sce,
variable = deng_SCE$cell_type2,
type = 'lineages',
col = 'black',
legend = TRUE)

sce <- slingshot::slingshot(deng_SCE,
clusterLabels = 'Seurat_clusters',
reducedDim = 'PCA',
end.clus = c("0", "3", "5")) ## check which would be the best according to bio

rm(list = ls())
gc()
.rs.restartR()

seu_int <- readRDS("seu_int_day2_part2.rds")

library(monocle3)

# get matrix and filter for minimum number of cells and features (the latter is a fix for backward compatibility)
mat_tmp <- seu_int@assays$RNA@counts
seu_tmp <- Seurat::CreateSeuratObject(mat_tmp, min.cells = 3,
min.features = 100)

feature_names <- as.data.frame(rownames(seu_tmp))
rownames(feature_names) <- rownames(seu_tmp)
colnames(feature_names) <- "gene_short_name"

seu_int_monocl <- monocle3::new_cell_data_set(seu_tmp@assays$RNA@counts,
cell_metadata = seu_int@meta.data,
gene_metadata = feature_names)

?preprocess_cds

seu_int_monocl <- monocle3::preprocess_cds(seu_int_monocl)

monocle3::plot_pc_variance_explained(seu_int_monocl)

seu_int_monocl <- monocle3::reduce_dimension(seu_int_monocl, reduction_method = "UMAP")

monocle3::plot_cells(seu_int_monocl, 
color_cells_by = "integrated_snn_res.0.3", 
cell_size = 1, 
show_trajectory_graph = FALSE)

monocle3::plot_cells(seu_int_monocl, genes = "CD79A", 
show_trajectory_graph = FALSE, 
cell_size = 1)


seu_int_monocl <- monocle3::cluster_cells(seu_int_monocl, resolution=0.00025)
monocle3::plot_cells(seu_int_monocl, label_cell_groups = F)

seu_int_monocl <- monocle3::learn_graph(seu_int_monocl)
monocle3::plot_cells(seu_int_monocl)

monocle3::plot_cells(seu_int_monocl, genes = c("CD79A", "CD34"),
show_trajectory_graph = FALSE, 
cell_size = 0.7, group_label_size = 4)

seu_int_monocl<-monocle3::order_cells(seu_int_monocl)#

monocle3::plot_cells(seu_int_monocl,
color_cells_by = "pseudotime",
label_cell_groups=F,
label_leaves=F,
label_branch_points=FALSE,
graph_label_size=1.5, cell_size = 1)

seuB <- choose_cells(seu_int_monocl)

plot_cells(seuB, show_trajectory_graph = FALSE, cell_size = 1)

pr_test <- graph_test(seuB, 
cores=4, 
neighbor_graph = "principal_graph")
# order by test statistic
pr_test <- pr_test[order(pr_test$morans_test_statistic, 
decreasing = TRUE),]
View(pr_test)

goi <- c("CD34", "MS4A1", "IGLL1", "IGLL5", 
"MKI67", "CKS2")
plot_cells(seuB, label_cell_groups=FALSE, genes = goi,
show_trajectory_graph=FALSE, cell_size = 1)

seuB@colData$monocle_cluster <- clusters(seuB)

plot_genes_in_pseudotime(subset(seuB, 
rowData(seuB)$gene_short_name %in% goi),
min_expr=0.5, color_cells_by = "monocle_cluster")
