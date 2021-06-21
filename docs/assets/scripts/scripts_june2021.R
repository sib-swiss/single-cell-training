# scRNAseq SIB training June 2021

### Day 1
library(Seurat)
library(scater)
library(SingleCellExperiment)
library(clustree)
library(celldex)
library(SingleR)
library(edgeR)
library(limma)
library(clusterProfiler)
library(enrichplot)
# BiocManager::install("org.Hs.eg.db", update = FALSE)
library(org.Hs.eg.db)
library(scater)
library(slingshot)
library(ggplot2)
library(ggbeeswarm)


# Import 10x data set
gbm.data <- Seurat::Read10X(data.dir = "data/data/gbm_dataset/filtered_feature_bc_matrix/")
gbm.data[c("PECAM1", "CD8A", "TSPAN1"), 1:30]
dim(gbm.data)
# 1] 36601  5604

# create seurat object:
gbm <- Seurat::CreateSeuratObject(counts = gbm.data,
                                  project = "gbm",
                                  min.cells = 3,
                                  min.features = 100)
gbm
# An object of class Seurat 
# 24363 features across 5553 samples within 1 assay 
# Active assay: RNA (24363 features, 0 variable features)
str(gbm)
head(gbm@meta.data)

hist(gbm@meta.data$nCount_RNA)
# scatter plot between 2 features:
Seurat::FeatureScatter(gbm, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

# percent mitochondrial gene expression
gbm[["percent.mt"]] <- Seurat::PercentageFeatureSet(gbm, pattern = "^MT-")
Seurat::VlnPlot(gbm, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

Seurat::FeatureScatter(gbm, feature1 = "nCount_RNA", feature2 = "percent.mt")

# Filtering of cells
gbm <- subset(gbm,
              subset = nFeature_RNA > 500 & nFeature_RNA < 7500 & percent.mt < 20)

# Normalization:
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

# scaling:
gbm <- Seurat::ScaleData(gbm,
                         features = rownames(gbm))

# cell cycle scoring using Seurat built-in cell cycle markers:
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

# Ridge plots of some cell cycle markers:
Seurat::RidgePlot(gbm, features = c("PCNA", "MKI67"),
                    group.by = "orig.ident",
                    ncol = 2)

# QC with scran and scater:
# extract counts to create sce:
cts <- Seurat::GetAssayData(gbm, slot = "counts")

gbm_sce <- SingleCellExperiment::SingleCellExperiment(
  assays = list(counts = cts),
  colData = gbm@meta.data,
  rowData = rownames(gbm)
)
class(gbm)
class(gbm_sce)
gbm_sce

# Import dissociation-related genes, ribosomal protein genes and mitochondrial genes:
dissoc_genes <- readLines("data/data/gbm_dataset/dissocation_genes.txt")
ribo_genes <- rownames(gbm)[grep(pattern = "^RP[S|L]", rownames(gbm), perl = T)]
mito_genes <- rownames(gbm)[grep(pattern = "^MT-", rownames(gbm))]

# calculate QC metrics using scran including the "QC"-related genes:
gbm_sce <- scuttle::addPerCellQC(gbm_sce,
                                 subsets=list(mito_genes=which(rownames(gbm_sce) %in% mito_genes),
                                              dissoc_genes=which(rownames(gbm_sce) %in% dissoc_genes),
                                              ribo_genes=which(rownames(gbm_sce) %in% ribo_genes)))

# check the meta data of the sce object:
SingleCellExperiment::colData(gbm_sce)

# create scater plots of QC metrics:
scater::plotColData(gbm_sce, x = "sum", y="detected")
scater::plotColData(gbm_sce, x = "detected", y="subsets_mito_genes_percent")
scater::plotColData(gbm_sce, x = "detected", y="subsets_dissoc_genes_percent")
scater::plotColData(gbm_sce, x = "subsets_mito_genes_percent", y="subsets_ribo_genes_percent")

# plot the top 30 most highly expressed genes:
scater::plotHighestExprs(gbm_sce, exprs_values = "counts", n = 30)

# alternative to Seurat's normalization here using scater
gbm_sce <- scater::logNormCounts(gbm_sce)

# calculate, for each gene, the variance it explains for several variables and plot
vars <- scater::getVarianceExplained(gbm_sce,
                                     variables = c("Phase", "nFeature_RNA"))
head(vars)
dim(vars )
scater::plotExplanatoryVariables(vars)

# saveRDS(gbm, "data/gbm_day1.rds")



#### Day 2  



# import data processed on day1:
# gbm <- readRDS("gbm_day1.rds")

# PCA:
gbm <- Seurat::RunPCA(gbm)
Seurat::DimPlot(gbm, reduction = "pca")
Seurat::DimPlot(gbm, reduction = "pca", group.by = "Phase")

# Create a PCA plot and color the cells according to gene expression:
FeaturePlot(gbm, reduction = "pca", features = "CLU")
FeaturePlot(gbm, reduction = "pca", features = "TYROBP")

# Heatmap of the genes correlating with the top 12 PCs:
Seurat::DimHeatmap(gbm, dims = 1:12, cells = 500, balanced = TRUE)
Seurat::ElbowPlot(gbm, ndims = 40)

# UMAP:
gbm <- Seurat::RunUMAP(gbm, dims = 1:25)

Seurat::DimPlot(gbm, reduction = "umap")
Seurat::DimPlot(gbm, reduction = "umap", group.by = "Phase")
Seurat::FeaturePlot(gbm, "nFeature_RNA")
Seurat::FeaturePlot(gbm, "percent.mt")
Seurat::FeaturePlot(gbm, "nCount_RNA")


# clustering:
gbm <- Seurat::FindNeighbors(gbm, dims = 1:25)
gbm <- Seurat::FindClusters(gbm, resolution = seq(0.1, 0.8, by=0.1))
head(gbm@meta.data)

clustree::clustree(gbm@meta.data[,grep("RNA_snn_res", colnames(gbm@meta.data))],
                   prefix = "RNA_snn_res.")

Seurat::DimPlot(gbm, group.by = "RNA_snn_res.0.1")
Seurat::DimPlot(gbm, group.by = "RNA_snn_res.0.1")


# cell type annotation
library(celldex)
library(SingleR)

# annotate cells clustering at res 0.2:
DimPlot(gbm, group.by = "RNA_snn_res.0.2")
gbm <- Seurat::SetIdent(gbm, value = gbm$RNA_snn_res.0.2)

DimPlot(gbm)
Seurat::FeaturePlot(gbm, "PMP2") # P2 protein is a constituent of the 
# myelin in nervous system

# Add module score:
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

?AddModuleScore

gbm <- Seurat::AddModuleScore(gbm,
                              features = list(immune_genes),
                              name = "immune_genes")
head(gbm@meta.data)
Seurat::FeaturePlot(gbm, "immune_genes1")
Seurat::VlnPlot(gbm,
                "immune_genes1")

# neuron markers downloaded from CellMarker
neurons<-read.csv("data/CellMarker_brain_neuron.csv", header = T)
dim(neurons)

head(neurons)

neurons<-subset(neurons, neurons$Cell.Marker %in% rownames(gbm))

gbm <- Seurat::AddModuleScore(gbm,
                              features = list(neurons$Cell.Marker),
                              name = "neuron_genes")
Seurat::FeaturePlot(gbm, "neuron_genes1",
                    cols = c("yellow", "darkblue")) + 
  labs(title="Neuron signature")
Seurat::VlnPlot(gbm, "neuron_genes1")

p0<-Seurat::DimPlot(gbm, group.by = "Phase")
p1<-Seurat::DimPlot(gbm, group.by = "RNA_snn_res.0.2")
p2<-Seurat::VlnPlot(gbm, "neuron_genes1", pt.size = 0)
cowplot::plot_grid(p0, p1, p2, nrow = 1)

# scatter plot of immune gene score vs neuron gene score:
FeatureScatter(gbm, "immune_genes1", "neuron_genes1")


# Annotate using singleR and the Human primary cell atlas data as a reference:
hpca.se <- celldex::HumanPrimaryCellAtlasData()
class(hpca.se)
head(colData(hpca.se))
table(hpca.se$label.main)

?SingleR
gbm_SingleR <- SingleR::SingleR(test = Seurat::GetAssayData(gbm, slot = "data"),
                                ref = hpca.se,
                                # clusters = gbm$RNA_snn_res.0.2,
                                labels = hpca.se$label.main)

head(gbm_SingleR)
# view single R scores for each cell for each reference cell type:
View(gbm_SingleR$scores)
table(gbm_SingleR$labels)

# add the cell type annotation to the gbm@meta.data:
gbm$SingleR_annot <- gbm_SingleR$labels
Seurat::DimPlot(gbm, group.by = "SingleR_annot", label = T, repel = T)

# calculate the average immune score per cell type:
mean_scores <- tapply(gbm$immune_genes1, gbm$SingleR_annot, mean)
mean_scores[order(mean_scores, decreasing = TRUE)]# [1:6]]

# check the distribution of cell types per cluster:
table(gbm$SingleR_annot, gbm$RNA_snn_res.0.2)

FeaturePlot(gbm, c("GFAP", "SOX9", "ATP13A4", "CBS")) # astrocyte markers according to CellMarker
head(neurons)
FeaturePlot(gbm, c("PDGFD", "RYR3", "SDK2", "STK32B")) # neuron markers according to CellMarker

# SingleR QC diagnostics:
# plot the distribution of  distribution of correlation scores per
# annotated cell, the cells labeled in orange are considered as bad quality annotation, 
# so will be labeled as NA in the "pruned.labels" column
plotDeltaDistribution(gbm_SingleR)

# generate a heatmap of all scores for each annotated cell. This may help in checking 
# if any group of cells has high scores for several reference cell types.
plotScoreHeatmap(gbm_SingleR)

# annotate the each cluster instead of each cell:
gbm_SingleR_cl <- SingleR::SingleR(test = Seurat::GetAssayData(gbm, slot = "data"),
                                ref = hpca.se,
                                clusters = gbm$RNA_snn_res.0.2,
                                labels = hpca.se$label.main)
gbm_SingleR_cl$scores
# check the score per cluster per cell type:
View(gbm_SingleR_cl$scores)
gbm_SingleR_cl$labels
# cluster 2 had an NA, so it's annotation was of bad quality:
gbm_SingleR_cl$pruned.labels

# saveRDS(gbm, "gbm_day2_part2.rds")


# Integration with pancreas data:
pancreas.data <- readRDS(file = "data/data/pancreas_dataset/pancreas_expression_matrix.rds")
metadata <- readRDS(file = "data/data/pancreas_dataset/pancreas_metadata.rds")

# normalize and generate UMAP
pancreas <- Seurat::CreateSeuratObject(pancreas.data, meta.data = metadata)
pancreas <- Seurat::NormalizeData(pancreas)
pancreas <- Seurat::FindVariableFeatures(pancreas, selection.method = "vst", nfeatures = 2000)
pancreas <- Seurat::ScaleData(pancreas)
pancreas <- Seurat::RunPCA(pancreas, npcs = 30)
pancreas <- Seurat::RunUMAP(pancreas, reduction = "pca", dims = 1:30)

Seurat::DimPlot(pancreas, reduction = "umap", group.by = "tech")
Seurat::DimPlot(pancreas, reduction = "umap", group.by = "celltype")

# # split according to scRNAseq technology:
pancreas.list <- Seurat::SplitObject(pancreas, split.by = "tech")
for (i in 1:length(pancreas.list)) {
  pancreas.list[[i]] <- Seurat::NormalizeData(pancreas.list[[i]])
  pancreas.list[[i]] <- Seurat::FindVariableFeatures(pancreas.list[[i]], selection.method = "vst", nfeatures = 2000,
                                                     verbose = FALSE)
}

# integrate:
pancreas.anchors <- Seurat::FindIntegrationAnchors(object.list = pancreas.list, dims = 1:30)

pancreas.integrated <- Seurat::IntegrateData(anchorset = pancreas.anchors, dims = 1:30)

# rerun UMAP and clustering:
Seurat::DefaultAssay(pancreas.integrated) <- "integrated"

pancreas.integrated <- Seurat::ScaleData(pancreas.integrated)
pancreas.integrated <- Seurat::RunPCA(pancreas.integrated, npcs = 30)
pancreas.integrated <- Seurat::RunUMAP(pancreas.integrated, reduction = "pca", dims = 1:30)

Seurat::DimPlot(pancreas.integrated, reduction = "umap", group.by = "tech")
Seurat::DimPlot(pancreas.integrated, reduction = "umap", group.by = "celltype", 
                label = TRUE, repel=T)

saveRDS(pancreas.integrated, "data/pancreas.integrated.rds")


### Day 3
# determine cluster markers: i.e. 1 cluster of cells against all others:
de_genes <- Seurat::FindAllMarkers(gbm,  min.pct = 0.25)
de_genes <- subset(de_genes, de_genes$p_val_adj < 0.05)
View(de_genes)

# check if the immune genes are part of the markers of any cluster:
immune_genes <- c("GZMA", "CD3E", "CD3D")
de_genes[de_genes$gene %in% immune_genes,]

# pairwise DEG between 2 cell types:
gbm <- Seurat::SetIdent(gbm, value = "SingleR_annot")
DEG_astro_vs_macro <- Seurat::FindMarkers(gbm,
                                          ident.1 = "Astrocyte",
                                          ident.2 = "Macrophage",
                                          group.by = gbm$SingleR_annot,
                                          test.use = "wilcox")
top_order <- order(DEG_astro_vs_macro$p_val_adj)
DEG_astro_vs_macro[top_order[1:10],]

Seurat::VlnPlot(gbm, features = "SLC2A5")

# DEG using edgeR and limma in the case of the pancreas data set, using tech as a covariate:
pancreas.integrated <- readRDS("pancreas.integrated.rds")
Seurat::DefaultAssay(pancreas.integrated) <- "RNA"
Seurat::Idents(pancreas.integrated) <- pancreas.integrated$celltype
Seurat::DimPlot(pancreas.integrated)

pancreas.dg <- subset(pancreas.integrated, idents = c("delta", "gamma"))
counts <- Seurat::GetAssayData(pancreas.dg, slot = "counts")
counts <- counts[rowSums(counts) != 0,]
dge <- edgeR::DGEList(counts = counts)
dge <- edgeR::calcNormFactors(dge)  

# technology (batch) used as covariate:
design <- stats::model.matrix(~ 0 + celltype + tech, data = pancreas.dg@meta.data)
# could also check for interaction:
# design <- model.matrix(~ 0 + gender + gender:genotype)
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

# compare limma to FindMarkers output
# save topTable output of limma to a data frame:
deg_limma_delta_vs_gamma<-limma::topTable(fit.contrasts, number = Inf, sort.by = "P")
deg_limma_delta_vs_gamma$gene<-rownames(deg_limma_delta_vs_gamma)
# count how many genes are DE:
length(which(deg_limma_delta_vs_gamma$adj.P.Val<=0.05)) #3118

# Run FindMarkers of Seurat to compare the same cell types as done with limma, 
# but this does not allow to include the technology as a covariate:
Seurat::DefaultAssay(pancreas.integrated) <- "RNA"
Seurat::Idents(pancreas.integrated) <- pancreas.integrated$celltype
Seurat::DimPlot(pancreas.integrated)

# compare wilcoxon to limma
wilcox_delta_gamma<-FindMarkers(pancreas.integrated,
                                ident.1 = "delta",
                                ident.2 = "gamma")
wilcox_delta_gamma<-subset(wilcox_delta_gamma, wilcox_delta_gamma$p_val_adj<=0.05)
wilcox_delta_gamma$gene<-rownames(wilcox_delta_gamma)
dim(wilcox_delta_gamma)
# 1133 6

# compare FindMarkers to limma + covariates by merging the 2 DE results:
compare_deg<-merge(wilcox_delta_gamma, deg_limma_delta_vs_gamma, all=T)
compare_deg$avg_log2FC<-ifelse(is.na(compare_deg$p_val_adj), 0, compare_deg$avg_log2FC)
compare_deg$logFC<-ifelse(compare_deg$adj.P.Val>0.05, 0, compare_deg$logFC)
par(mar=c(5,5,5,5))
plot(compare_deg$avg_log2FC, compare_deg$logFC,
     xlab="Log2FC FindMarkers (wilcoxon)",
     ylab="Log2FC limma with covariates",
     pch=15, cex=0.5)
abline(a=0, b = 1, col="red")

# saveRDS(de_genes, "data/de_genes.rds")

# saveRDS(DEG_astro_vs_macro, "data/DEG_astro_vs_macro.rds")

# Over-representation analysis:
library(clusterProfiler)
AnnotationDbi::keytypes(org.Hs.eg.db)
AC_up_DEG <- subset(DEG_astro_vs_macro,
                    DEG_astro_vs_macro$avg_log2FC > 0 &
                      DEG_astro_vs_macro$p_val_adj < 0.05)
AC_up_genes <- rownames(AC_up_DEG)

# Over-rep analysis of GO terms:
AC_MAC_GO <- clusterProfiler::enrichGO(AC_up_genes, # vector of up regulated genes
                                       "org.Hs.eg.db", # orgdb= package that contains gene label types correspondances
                                       keyType = "SYMBOL", # indicate that genes are labeled using symbols
                                       ont = "BP", # which of the GO categories to test, here the "Biological Processes"
                                       minGSSize = 50) # exclude gene sets that contain less than 50 genes
View(AC_MAC_GO@result)

# genes contained within each GO term as included in the output as a list:
AC_MAC_GO@geneSets
AC_MAC_GO@geneSets$`GO:0000002`

# network plot of top 30 significant gene sets:
enrichplot::emapplot(enrichplot::pairwise_termsim(AC_MAC_GO),
                     showCategory = 30, cex_label_category = 0.5)

# export results to a text file to visualize the results online with Revigo:
write.table(AC_MAC_GO@result, "data/AC_MAC_GO.txt",
            sep="\t", row.names = F, quote = F)


# hallmark ORA:
gmt <- clusterProfiler::read.gmt("data/data/gbm_dataset/h.all.v7.2.symbols.xls")
head(gmt)

AC_MAC_enrich <- clusterProfiler::enricher(gene = AC_up_genes,
                                           universe = rownames(gbm),
                                           pAdjustMethod = "BH",
                                           pvalueCutoff  = 0.05,
                                           qvalueCutoff  = 0.05,
                                           TERM2GENE = gmt)
View(AC_MAC_enrich@result)

# extract the gene list for the myc target gene set:
myc_target_genes <- gmt$gene[gmt$term=="HALLMARK_MYC_TARGETS_V1"]

# extract only the genes that are up-regulated and part of the myc target gene set:
myc_target_genes_up<-unlist(strsplit(AC_MAC_enrich@result$geneID[1], "\\/"))

# calculate module score for all genes included in the myc target gene set:
gbm <- Seurat::AddModuleScore(gbm,
                              features = list(myc_target_genes=myc_target_genes),
                              name = "myc_target_genes")
Seurat::VlnPlot(gbm, "myc_target_genes1", group.by = "SingleR_annot",
                pt.size = 0)

# calculate the module score only for genes that are part of the myc target gene set AND up-regulated:
gbm <- Seurat::AddModuleScore(gbm,
                              features = list(myc_target_genes=myc_target_genes_up),
                              name = "myc_target_genes_up")

Seurat::VlnPlot(gbm, "myc_target_genes_up1", pt.size = 0)
# plot cell type annotation and module score side by side:
p1<-DimPlot(gbm, group.by = "SingleR_annot")
p2<-FeaturePlot(gbm, "myc_target_genes_up1", col=c("grey98", "darkblue"))
cowplot::plot_grid(p1, p2, nrow = 1)


# trajectory analysis using slingshot
library(SingleCellExperiment)
library(scater)
library(slingshot)
library(ggplot2)
library(ggbeeswarm)

# import the SingleCellExperiment object:
deng_SCE <- readRDS("/export/scratch/twyss/SIB_scRNAseq_course/June2021/data/data/deng_dataset/deng-reads.rds")
# re-organize the cell type factor to a chronological order instead of alphabetical order:
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
# PCA
deng_SCE <- scater::runPCA(deng_SCE, ncomponents = 50)
pca <- SingleCellExperiment::reducedDim(deng_SCE, "PCA")

head(pca)
deng_SCE$PC1 <- pca[, 1]
deng_SCE$PC2 <- pca[, 2]

# plot PCA
ggplot(as.data.frame(colData(deng_SCE)), aes(x = PC1, y = PC2, color = cell_type2)) +
  geom_point(size=2, shape=20) +
  theme_classic() +
  xlab("PC1") + ylab("PC2") + ggtitle("PC biplot")

# rank cells according to their PC1 score
deng_SCE$pseudotime_PC1 <- rank(deng_SCE$PC1)  

# plot cells with the PC1 score on the x-axis and the cell type on the y axis:
ggplot(as.data.frame(colData(deng_SCE)), aes(x = pseudotime_PC1, y = cell_type2,
                                             colour = cell_type2)) +
  ggbeeswarm::geom_quasirandom(groupOnX = FALSE) +
  theme_classic() +
  xlab("PC1") + ylab("Timepoint") +
  ggtitle("Cells ordered by first principal component")

# calculate slingshot trajectory without providing a cell label
sce <- slingshot::slingshot(deng_SCE, reducedDim = 'PCA')

# Geert's function to create a plot using a slingshot object:
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

# Create PCA coloring the cells according to pseudotime and adding the 
# trajectory on top:
PCAplot_slingshot(sce, variable = sce$slingPseudotime_1, draw_lines = TRUE)

# plot the cells with the pseudotime on the x axis and the cell type on the y axis:
ggplot(as.data.frame(colData(deng_SCE)), aes(x = sce$slingPseudotime_1,
                                             y = cell_type2,
                                             colour = cell_type2)) +
  ggbeeswarm::geom_quasirandom(groupOnX = FALSE) +
  theme_classic() +
  xlab("Slingshot pseudotime") + ylab("Timepoint") +
  ggtitle("Cells ordered by Slingshot pseudotime")

# use Seurat to perform clustering :
gcdata <- Seurat::CreateSeuratObject(counts = SingleCellExperiment::counts(deng_SCE),
                                     project = "slingshot")

gcdata <- Seurat::NormalizeData(object = gcdata,
                                normalization.method = "LogNormalize",
                                scale.factor = 10000)

gcdata <- Seurat::FindVariableFeatures(object = gcdata)# ,
                                       # mean.function = ExpMean,
                                       # dispersion.function = LogVMR)

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
gcdata$cell_type<-sapply(strsplit(rownames(gcdata@meta.data), "\\."), "[",1)
Seurat::DimPlot(gcdata)
Seurat::DimPlot(gcdata, group.by = "cell_type")

# add cluster ID of each cell to singleCellExperiment object:
head(colData(deng_SCE))
deng_SCE$Seurat_clusters <- as.character(gcdata$seurat_clusters)  # go from factor to character
head(colData(deng_SCE))

# run slingshot including information about cell cluster ID:
deng_SCE  <- slingshot::slingshot(deng_SCE,
                                  clusterLabels = "Seurat_clusters",
                                  reducedDim = "PCA",
                                  start.clus = "2")

head(colData(deng_SCE))
# now slingshot calculated 2 trajectories:
SlingshotDataSet(deng_SCE)

# create a PCA plot coloring the cells according to pseudotime and adding the 2 trajectories on top:
PCAplot_slingshot(deng_SCE, variable = deng_SCE$slingPseudotime_2)

# plot the cells with pseudotime 1 on the x-axis and cell type on the y-axis:
ggplot(data.frame(cell_type2 = deng_SCE$cell_type2,
                  slingPseudotime_1 = deng_SCE$slingPseudotime_1),
       aes(x = slingPseudotime_1, y = cell_type2,
           colour = cell_type2)) +
  ggbeeswarm::geom_quasirandom(groupOnX = FALSE) +
  theme_classic() +
  xlab("Slingshot pseudotime") + ylab("Timepoint") +
  ggtitle("Cells ordered by Slingshot pseudotime")

# plot the cells with pseudotime 2 on the x-axis and cell type on the y-axis:
ggplot(data.frame(cell_type2 = deng_SCE$cell_type2,
                  slingPseudotime_2 = deng_SCE$slingPseudotime_2),
       aes(x = slingPseudotime_2, y = cell_type2,
           colour = cell_type2)) +
  ggbeeswarm::geom_quasirandom(groupOnX = FALSE) +
  theme_classic() +
  xlab("Slingshot pseudotime") + ylab("Timepoint") +
  ggtitle("Cells ordered by Slingshot pseudotime")

# PCA plot with cells colored according to cluster ID and lineage added on top:
PCAplot_slingshot(deng_SCE,
                  variable = deng_SCE$Seurat_clusters,
                  type = 'lineages',
                  col = 'black',
                  legend = TRUE)

# PCA plot with cells colored according to cell type and lineage added on top:
PCAplot_slingshot(deng_SCE,
                  variable = deng_SCE$cell_type2,
                  type = 'lineages',
                  col = 'black',
                  legend = TRUE)

# Monocle3:
library(monocle3)

DimPlot(gbm, group.by = "SingleR_annot")

# create a cell_data_set type of object specific for monocle3
feature_names <- as.data.frame(rownames(gbm))
rownames(feature_names) <- rownames(gbm)
colnames(feature_names) <- "gene_short_name"
gbm_monocl <- monocle3::new_cell_data_set(gbm@assays$RNA@counts,
                                          cell_metadata = gbm@meta.data,
                                          gene_metadata = feature_names)
# normalize and perform PCA
gbm_monocl <- monocle3::preprocess_cds(gbm_monocl)
# elbow plot:
monocle3::plot_pc_variance_explained(gbm_monocl)
# UMAP:
gbm_monocl <- monocle3::reduce_dimension(gbm_monocl, reduction_method = "UMAP")
# plot UMAP coloring the cells according to Seurat clustering:
monocle3::plot_cells(gbm_monocl, color_cells_by = "RNA_snn_res.0.2")
# to UMAP coloring the cells according to expression level of a gene
monocle3::plot_cells(gbm_monocl, genes = "PMP2") 

# cluster cells using monocl algorithm and plot UMAP with monocl and seurat clusters side by side:
gbm_monocl <- monocle3::cluster_cells(gbm_monocl, resolution=0.00025)
p1 <- monocle3::plot_cells(gbm_monocl, label_cell_groups = F)
p2 <- monocle3::plot_cells(gbm_monocl, color_cells_by = "RNA_snn_res.0.2", label_cell_groups = F)
cowplot::plot_grid(p1, p2, ncol = 2) # Are there differences?

# learn trajectory using monocl clustering ID and plot
gbm_monocl <- monocle3::learn_graph(gbm_monocl)
monocle3::plot_cells(gbm_monocl) + labs(title="Monocle3 clustering")
monocle3::plot_cells(gbm_monocl, color_cells_by = "RNA_snn_res.0.2") + labs(title="Seurat clustering res.0.2")

# replace cluster ids by Seurat cluster IDs:
gbm_monocl@clusters$UMAP$clusters <- colData(gbm_monocl)$RNA_snn_res.0.2
names(gbm_monocl@clusters$UMAP$clusters) <- rownames(colData(gbm_monocl))
gbm_monocl <- monocle3::learn_graph(gbm_monocl)
monocle3::plot_cells(gbm_monocl, label_cell_groups = F)
monocle3::plot_cells(gbm_monocl, label_cell_groups = F, color_cells_by = "SingleR_annot")

# determine which is the root (i.e. starting point) of each trajectory to 
# calculate pseudotime for each cell:
gbm_monocl<-monocle3::order_cells(gbm_monocl)#

# plot UMAP coloring cells according to  pseudotime 
monocle3::plot_cells(gbm_monocl,
                     color_cells_by = "pseudotime",
                     label_cell_groups=F,
                     label_leaves=F,
                     label_branch_points=FALSE,
                     graph_label_size=1.5, cell_size = 1)

# plot gene expression vs pseudotime:
plot_genes_in_pseudotime(subset(gbm_monocl, rowData(gbm_monocl)$gene_short_name=="GFAP"))
