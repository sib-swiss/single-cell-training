# Auto-extracted code from markdown files

## Code found in: index.md


## Code found in: precourse.md


## Code found in: course_schedule.md


## Code found in: day1/general_introduction.md


## Code found in: day1/setup.md


## Code found in: day1/introduction_scrnaseq.md
library(Seurat)

gbm.data <- Seurat::Read10X(data.dir = "data/gbm_dataset/filtered_feature_bc_matrix/")

gbm.data[c("PECAM1", "CD8A", "TSPAN1"), 1:30]

gbm <- Seurat::CreateSeuratObject(counts = gbm.data,
project = "gbm",
min.cells = 3,
min.features = 100)

hist(gbm$nCount_RNA)

hist(gbm@meta.data$nCount_RNA)

Seurat::FeatureScatter(gbm, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

gbm[["percent.mt"]] <- Seurat::PercentageFeatureSet(gbm, pattern = "^MT-")

head(gbm@meta.data)

Seurat::VlnPlot(gbm, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

Seurat::FeatureScatter(gbm, feature1 = "nCount_RNA", feature2 = "percent.mt")

gbm <- subset(gbm,
subset = nFeature_RNA > 500 & nFeature_RNA < 7500 & percent.mt < 20)

## Code found in: day1/normalization_scaling.md
Seurat::GetAssay(gbm)[1:10,1:10]  

gbm <- Seurat::NormalizeData(gbm,
normalization.method = "LogNormalize",
scale.factor = 10000)

gbm <- Seurat::FindVariableFeatures(gbm,
selection.method = "vst",
nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(Seurat::VariableFeatures(gbm), 10)
top10

vf_plot <- Seurat::VariableFeaturePlot(gbm)
Seurat::LabelPoints(plot = vf_plot,
points = top10, repel = TRUE)

gbm <- Seurat::ScaleData(gbm,
features = rownames(gbm))

## Code found in: day1/quality_control.md
Seurat::cc.genes.updated.2019

s.genes <- Seurat::cc.genes.updated.2019$s.genes
g2m.genes <- Seurat::cc.genes.updated.2019$g2m.genes

gbm <- Seurat::CellCycleScoring(gbm,
s.features = s.genes,
g2m.features = g2m.genes)

head(gbm)
table(gbm$Phase)
#   G1  G2M    S
# 2887  711 1493

Seurat::RidgePlot(gbm, features = c("PCNA", "MKI67"),
group.by = "orig.ident",
ncol = 2)

library(scater)
library(SingleCellExperiment)

cts <- Seurat::GetAssayData(gbm, slot = "counts")

gbm_sce <- SingleCellExperiment::SingleCellExperiment(
assays = list(counts = cts),
colData = gbm@meta.data,
rowData = rownames(gbm)
)
class(gbm)
class(gbm_sce)
gbm_sce

dissoc_genes <- readLines("data/gbm_dataset/dissocation_genes.txt")
ribo_genes <- rownames(gbm)[grep(pattern = "^RP[S|L]", rownames(gbm), perl = T)]
mito_genes <- rownames(gbm)[grep(pattern = "^MT-", rownames(gbm))]

gbm_sce <- scuttle::addPerCellQC(gbm_sce,
subsets=list(mito_genes=which(rownames(gbm_sce) %in% mito_genes),
dissoc_genes=which(rownames(gbm_sce) %in% dissoc_genes),
ribo_genes=which(rownames(gbm_sce) %in% ribo_genes)))

SingleCellExperiment::colData(gbm_sce)

scater::plotColData(gbm_sce, x = "sum", y="detected")
scater::plotColData(gbm_sce, x = "detected", y="subsets_mito_genes_percent")
scater::plotColData(gbm_sce, x = "detected", y="subsets_dissoc_genes_percent")
scater::plotColData(gbm_sce, x = "subsets_mito_genes_percent", y="subsets_ribo_genes_percent")

scater::plotHighestExprs(gbm_sce, exprs_values = "counts", n = 30)

gbm_sce <- scater::logNormCounts(gbm_sce)  # alternative to Seurat's normalization here using scater

vars <- scater::getVarianceExplained(gbm_sce,
variables = "Phase")
head(vars)

scater::plotExplanatoryVariables(vars)

saveRDS(gbm, "gbm_day1.rds")

rm(list = ls())
gc()
.rs.restartR()

## Code found in: day2/dimensionality_reduction.md
gbm <- readRDS("gbm_day1.rds")

library(Seurat)
library(clustree)

gbm <- Seurat::RunPCA(gbm)

Seurat::DimPlot(gbm, reduction = "pca")

Seurat::DimPlot(gbm, reduction = "pca", group.by = "Phase")

Seurat::DimHeatmap(gbm, dims = 1:12, cells = 500, balanced = TRUE)

Seurat::ElbowPlot(gbm, ndims = 40)

gbm <- Seurat::RunUMAP(gbm, dims = 1:25)

Seurat::DimPlot(gbm, reduction = "umap")

Seurat::DimPlot(gbm, reduction = "umap", group.by = "Phase")

## Code found in: day2/clustering.md
gbm <- Seurat::FindNeighbors(gbm, dims = 1:25)

gbm <- Seurat::FindClusters(gbm, resolution = seq(0.1, 0.8, by=0.1))

head(gbm@meta.data)

library(clustree)
clustree::clustree(gbm@meta.data[,grep("RNA_snn_res", colnames(gbm@meta.data))],
prefix = "RNA_snn_res.")

Seurat::DimPlot(gbm, group.by = "RNA_snn_res.0.1")

Seurat::DimPlot(gbm, group.by = "RNA_snn_res.0.2")

saveRDS(gbm, "gbm_day2_part1.rds")

rm(list = ls())
gc()
.rs.restartR()

## Code found in: day2/cell_annotation.md
gbm <- readRDS("gbm_day2_part1.rds")

library(celldex)
library(SingleR)

gbm <- Seurat::SetIdent(gbm, value = gbm$RNA_snn_res.0.2)

Seurat::FeaturePlot(gbm, "PMP2")

immune_genes<-c("GZMA", "CD3E", "CD3D")
microglia_genes<-c("CCL4", "CCL3", "P2RY12", "C1QB", "CSF1ER", "CY3CR1")

Seurat::FeaturePlot(gbm, immune_genes, ncol=2)

Seurat::VlnPlot(gbm,
features = immune_genes,
ncol = 2)

Seurat::FeaturePlot(gbm, microglia_genes, ncol=2)

Seurat::VlnPlot(gbm,
features = microglia_genes,
ncol = 2)

gbm <- Seurat::AddModuleScore(gbm,
features = list(immune_genes),
name = "immune_genes")

Seurat::FeaturePlot(gbm, "immune_genes1")

Seurat::VlnPlot(gbm,
"immune_genes1")

hpca.se <- celldex::HumanPrimaryCellAtlasData()
class(hpca.se)
table(hpca.se$label.main)

gbm_SingleR <- SingleR::SingleR(test = Seurat::GetAssayData(gbm, slot = "data"),
ref = hpca.se,
labels = hpca.se$label.main)

head(gbm_SingleR)

gbm$SingleR_annot <- gbm_SingleR$labels

Seurat::DimPlot(gbm, group.by = "SingleR_annot", label = T, repel = T)

mean_scores <- tapply(gbm$immune_genes1, gbm$SingleR_annot, mean)
mean_scores[order(mean_scores, decreasing = TRUE)[1:6]]

saveRDS(gbm, "gbm_day2_part2.rds")

rm(list = ls())
gc()
.rs.restartR()

## Code found in: day2/integration.md
pancreas.data <- readRDS(file = "data/pancreas_dataset/pancreas_expression_matrix.rds")
metadata <- readRDS(file = "data/pancreas_dataset/pancreas_metadata.rds")

pancreas <- Seurat::CreateSeuratObject(pancreas.data, meta.data = metadata)

pancreas <- Seurat::NormalizeData(pancreas)
pancreas <- Seurat::FindVariableFeatures(pancreas, selection.method = "vst", nfeatures = 2000)
pancreas <- Seurat::ScaleData(pancreas)
pancreas <- Seurat::RunPCA(pancreas, npcs = 30)
pancreas <- Seurat::RunUMAP(pancreas, reduction = "pca", dims = 1:30)

Seurat::DimPlot(pancreas, reduction = "umap", group.by = "tech")

Seurat::DimPlot(pancreas, reduction = "umap", group.by = "celltype")

pancreas.list <- Seurat::SplitObject(pancreas, split.by = "tech")

for (i in 1:length(pancreas.list)) {
pancreas.list[[i]] <- Seurat::NormalizeData(pancreas.list[[i]])
pancreas.list[[i]] <- Seurat::FindVariableFeatures(pancreas.list[[i]], selection.method = "vst", nfeatures = 2000,
verbose = FALSE)
}

pancreas.anchors <- Seurat::FindIntegrationAnchors(object.list = pancreas.list, dims = 1:30)

pancreas.integrated <- Seurat::IntegrateData(anchorset = pancreas.anchors, dims = 1:30)

Seurat::DefaultAssay(pancreas.integrated) <- "integrated"

pancreas.integrated <- Seurat::ScaleData(pancreas.integrated)
pancreas.integrated <- Seurat::RunPCA(pancreas.integrated, npcs = 30)
pancreas.integrated <- Seurat::RunUMAP(pancreas.integrated, reduction = "pca", dims = 1:30)

Seurat::DimPlot(pancreas.integrated, reduction = "umap", group.by = "tech")
Seurat::DimPlot(pancreas.integrated, reduction = "umap", group.by = "celltype", label = TRUE, repel = TRUE)

saveRDS(pancreas.integrated, "pancreas.integrated.rds")

rm(list = ls())
gc()
.rs.restartR()

## Code found in: day3/differential_gene_expression.md
gbm <- readRDS("gbm_day2_part2.rds")

library(Seurat)
library(edgeR)
library(limma)

de_genes <- Seurat::FindAllMarkers(gbm,  min.pct = 0.25)

de_genes <- subset(de_genes, de_genes$p_val_adj < 0.05)
View(de_genes)

immune_genes <- c("GZMA", "CD3E", "CD3D")

de_genes[de_genes$gene %in% immune_genes,]

gbm <- Seurat::SetIdent(gbm, value = "SingleR_annot")

DEG_astro_vs_macro <- Seurat::FindMarkers(gbm,
ident.1 = "Astrocyte",
ident.2 = "Macrophage",
group.by = gbm$SingleR_annot,
test.use = "wilcox")

top_order <- order(DEG_astro_vs_macro$p_val_adj)
DEG_astro_vs_macro[top_order[1:10],]

Seurat::VlnPlot(gbm, features = "SLC2A5")

pancreas.integrated <- readRDS("pancreas.integrated.rds")

Seurat::DefaultAssay(pancreas.integrated) <- "RNA"
Seurat::Idents(pancreas.integrated) <- pancreas.integrated$celltype

Seurat::DimPlot(pancreas.integrated)

pancreas.dg <- subset(pancreas.integrated, idents = c("delta", "gamma"))

counts <- Seurat::GetAssayData(pancreas.dg, slot = "counts")
counts <- counts[rowSums(counts) != 0,]

dge <- edgeR::DGEList(counts = counts)
dge <- edgeR::calcNormFactors(dge)  

design <- model.matrix(~ 0 + celltype + tech, data = pancreas.dg@meta.data)
colnames(design)<-c("delta", "gamma", "celseq2", "fluidigmc1", "smartseq2")

contrast.mat <- limma::makeContrasts(delta - gamma,
levels = design)

vm <- limma::voom(dge, design = design, plot = TRUE)
fit <- limma::lmFit(vm, design = design)
fit.contrasts <- limma::contrasts.fit(fit, contrast.mat)
fit.contrasts <- limma::eBayes(fit.contrasts)

limma::topTable(fit.contrasts, number = 10, sort.by = "P")

Seurat::VlnPlot(pancreas.dg, "PPY", split.by = "tech")
Seurat::VlnPlot(pancreas.dg, "RBP4", split.by = "tech")

## Code found in: day3/enrichment_analysis.md
library(clusterProfiler)
library(enrichplot)

BiocManager::install("org.Hs.eg.db", update = FALSE)
library(org.Hs.eg.db)
AnnotationDbi::keytypes(org.Hs.eg.db)

AC_up_DEG <- subset(DEG_astro_vs_macro,
DEG_astro_vs_macro$avg_log2FC > 0 &
DEG_astro_vs_macro$p_val_adj < 0.05)
AC_up_genes <- rownames(AC_up_DEG)

AC_MAC_GO <- clusterProfiler::enrichGO(AC_up_genes, # vector of up regulated genes
"org.Hs.eg.db", # orgdb= package that contains gene label types correspondances
keyType = "SYMBOL", # indicate that genes are labeled using symbols
ont = "BP", # which of the GO categories to test, here the "Biological Processes"
minGSSize = 50) # exclude gene sets that contain less than 50 genes

View(AC_MAC_GO@result)

enrichplot::emapplot(enrichplot::pairwise_termsim(AC_MAC_GO),
showCategory = 30, cex_label_category = 0.5)

gmt <- clusterProfiler::read.gmt("data/gbm_dataset/h.all.v7.2.symbols.xls")
head(gmt)

AC_MAC_enrich <- clusterProfiler::enricher(gene = AC_up_genes,
universe = rownames(gbm),
pAdjustMethod = "BH",
pvalueCutoff  = 0.05,
qvalueCutoff  = 0.05,
TERM2GENE = gmt)

View(AC_MAC_enrich@result)

myc_target_genes <- gmt$gene[gmt$term=="HALLMARK_MYC_TARGETS_V1"]

gbm <- Seurat::AddModuleScore(gbm,
features = list(myc_target_genes=myc_target_genes),
name = "myc_target_genes")

Seurat::VlnPlot(gbm, "myc_target_genes1", group.by = "SingleR_annot")

saveRDS(gbm, "gbm_day3.rds")

rm(list = ls())
gc()
.rs.restartR()

## Code found in: day3/trajectory_analysis.md
library(SingleCellExperiment)
library(scater)
library(slingshot)
library(ggplot2)
library(ggbeeswarm)

deng_SCE <- readRDS("data/deng_dataset/deng-reads.rds")

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

deng_SCE <- slingshot::slingshot(deng_SCE,
clusterLabels = 'Seurat_clusters',
reducedDim = 'PCA',
start.clus = "2")

head(colData(deng_SCE))

SlingshotDataSet(deng_SCE)

PCAplot_slingshot(deng_SCE, variable = deng_SCE$slingPseudotime_2)

ggplot(data.frame(cell_type2 = deng_SCE$cell_type2,
slingPseudotime_1 = deng_SCE$slingPseudotime_1),
aes(x = slingPseudotime_1, y = cell_type2,
colour = cell_type2)) +
ggbeeswarm::geom_quasirandom(groupOnX = FALSE) +
theme_classic() +
xlab("Slingshot pseudotime") + ylab("Timepoint") +
ggtitle("Cells ordered by Slingshot pseudotime")

ggplot(data.frame(cell_type2 = deng_SCE$cell_type2,
slingPseudotime_2 = deng_SCE$slingPseudotime_2),
aes(x = slingPseudotime_2, y = cell_type2,
colour = cell_type2)) +
ggbeeswarm::geom_quasirandom(groupOnX = FALSE) +
theme_classic() +
xlab("Slingshot pseudotime") + ylab("Timepoint") +
ggtitle("Cells ordered by Slingshot pseudotime")

PCAplot_slingshot(deng_SCE,
variable = deng_SCE$Seurat_clusters,
type = 'lineages',
col = 'black',
legend = TRUE)

PCAplot_slingshot(deng_SCE,
variable = deng_SCE$cell_type2,
type = 'lineages',
col = 'black',
legend = TRUE)

rm(list = ls())
gc()
.rs.restartR()

gbm <- readRDS("gbm_day3.rds")

library(monocle3)

feature_names <- as.data.frame(rownames(gbm))
rownames(feature_names) <- rownames(gbm)
colnames(feature_names) <- "gene_short_name"
gbm_monocl <- monocle3::new_cell_data_set(gbm@assays$RNA@counts,
cell_metadata = gbm@meta.data,
gene_metadata = feature_names)

?preprocess_cds

gbm_monocl <- monocle3::preprocess_cds(gbm_monocl)

monocle3::plot_pc_variance_explained(gbm_monocl)

gbm_monocl <- monocle3::reduce_dimension(gbm_monocl, reduction_method = "UMAP")

monocle3::plot_cells(gbm_monocl, color_cells_by = "RNA_snn_res.0.2")
monocle3::plot_cells(gbm_monocl, genes = "PMP2") # to plot expression level of a gene

gbm_monocl <- monocle3::cluster_cells(gbm_monocl, resolution=0.00025)
p1 <- monocle3::plot_cells(gbm_monocl, label_cell_groups = F)
p2 <- monocle3::plot_cells(gbm_monocl, color_cells_by = "RNA_snn_res.0.2", label_cell_groups = F)
cowplot::plot_grid(p1, p2, ncol = 2) # Are there differences?

gbm_monocl <- monocle3::learn_graph(gbm_monocl)
monocle3::plot_cells(gbm_monocl)
monocle3::plot_cells(gbm_monocl, color_cells_by = "RNA_snn_res.0.2")

gbm_monocl@clusters$UMAP$clusters <- colData(gbm_monocl)$RNA_snn_res.0.2
names(gbm_monocl@clusters$UMAP$clusters) <- rownames(colData(gbm_monocl))
gbm_monocl <- monocle3::learn_graph(gbm_monocl)
monocle3::plot_cells(gbm_monocl, label_cell_groups = F)

gbm_monocl<-monocle3::order_cells(gbm_monocl)#

monocle3::plot_cells(gbm_monocl,
color_cells_by = "pseudotime",
label_cell_groups=F,
label_leaves=F,
label_branch_points=FALSE,
graph_label_size=1.5, cell_size = 1)

plot_genes_in_pseudotime(subset(gbm_monocl, rowData(gbm_monocl)$gene_short_name=="PMP2"))
