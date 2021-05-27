library(Seurat)

gbm.data <- Seurat::Read10X(data.dir = "data/filtered_feature_bc_matrix/")

gbm <- Seurat::CreateSeuratObject(counts = gbm.data,
                                  project = "gbm",
                                  min.cells = 3,
                                  min.features = 100)

gbm[["percent.mt"]] <- Seurat::PercentageFeatureSet(gbm, pattern = "^MT-")
Seurat::VlnPlot(gbm, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

Seurat::FeatureScatter(gbm, feature1 = "nCount_RNA", feature2 = "percent.mt")

gbm <- subset(gbm,
              subset = nFeature_RNA > 500 & nFeature_RNA < 7500 & percent.mt < 20)

Seurat::GetAssay(gbm)[1:10,1:10]

gbm <- Seurat::NormalizeData(gbm,
                             normalization.method = "LogNormalize",
                             scale.factor = 10000)


gbm <- FindVariableFeatures(gbm, 
                            selection.method = "vst", 
                            nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(gbm), 10)



vf_plot <- Seurat::VariableFeaturePlot(gbm)
Seurat::LabelPoints(plot = vf_plot,
                    points = top10, repel = TRUE)

gbm <- Seurat::ScaleData(gbm,
                         features = rownames(gbm))

cc.genes.updated.2019

s.genes <- Seurat::cc.genes.updated.2019$s.genes
g2m.genes <- Seurat::cc.genes.updated.2019$g2m.genes

gbm <- Seurat::CellCycleScoring(gbm, s.features = s.genes, g2m.features = g2m.genes)

Seurat::RidgePlot(gbm, features = c("PCNA", "MKI67"),
                  group.by = "orig.ident",
                  ncol = 2)

cts <- Seurat::GetAssayData(gbm, slot = "counts")

gbm_sce <- SingleCellExperiment::SingleCellExperiment(
  assays = list(counts = cts), 
  colData = gbm@meta.data, 
  rowData = rownames(gbm)
)
gbm_sce

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

writeLines(dissoc_genes, "data/dissocation_genes.txt")

dissoc_genes <- readLines("data/dissocation_genes.txt")

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

gbm_sce <- scater::logNormCounts(gbm_sce) 

vars <- scater::getVarianceExplained(gbm_sce,
                                     variables = "Phase"
                                     #, variables=c("tissue", "detected", "sex", "age")
)
head(vars)

scater::plotExplanatoryVariables(vars) 

gbm <- RunPCA(gbm) 
DimPlot(gbm, reduction = "pca", group.by = "Phase")

Seurat::DimHeatmap(gbm, dims = 1:12, cells = 500, balanced = TRUE) 

gbm <- Seurat::RunUMAP(gbm, dims = 1:25)

Seurat::DimPlot(gbm, reduction = "umap")
Seurat::DimPlot(gbm, reduction = "umap", group.by = "Phase")

gbm <- Seurat::FindNeighbors(gbm, dims = 1:25)

gbm <- Seurat::FindClusters(gbm, resolution = seq(0.1, 0.8, by=0.1))

library(clustree)
clustree::clustree(gbm@meta.data[,grep("RNA_snn_res", colnames(gbm@meta.data))],
                   prefix = "RNA_snn_res.")
Seurat::DimPlot(gbm, group.by = "RNA_snn_res.0.1")
png("UMAP_res.0.2.png")
Seurat::DimPlot(gbm, group.by = "RNA_snn_res.0.2")
dev.off()
Seurat::DimPlot(gbm, group.by = "RNA_snn_res.0.3")

FeaturePlot(gbm, "PMP2")

endothelial_genes<-c("RGS5", "VWF", "CDH5", "PTPRB", "CD34", "ABCB1")
astrocyte_genes<-c("PMP2", "AQP4", "SLC1A2", "GJA1", "GJB6", "SLC4A4")
neural_genes<-c("VIP", "RELN", "GAD2", "GAD1", "SCG2", "SYNPR")
immune_genes<-c("GZMA", "CD3E", "CD3D")
microglia_genes<-c("CCL4", "CCL3", "P2RY12", "C1QB", "CSF1ER", "CY3CR1")

FeaturePlot(gbm, immune_genes, 
            ncol=2)
Seurat::VlnPlot(gbm, features = immune_genes, ncol = 2, group.by = "RNA_snn_res.0.2")
png("UMAP_res.0.2_labels.png")
DimPlot(gbm, group.by = "RNA_snn_res.0.2", label = TRUE, repel = TRUE)
dev.off()

png("featureplots_microglia.png", width = 800, height = 800)
FeaturePlot(gbm, microglia_genes, ncol = 2)
dev.off()

png("violinplots_microglia.png")
Seurat::VlnPlot(gbm, features = microglia_genes, ncol = 2, group.by = "RNA_snn_res.0.2")
dev.off()

gbm <- Seurat::AddModuleScore(gbm, features = list(immune_genes), name = "immune_genes")

png("UMAP_immune_genes_modules.png")
Seurat::FeaturePlot(gbm, "immune_genes1")
dev.off()

png("violinplot_immune_genes_modules.png")
Seurat::VlnPlot(gbm, "immune_genes1") 
dev.off()

hpca.se <- celldex::HumanPrimaryCellAtlasData()
table(hpca.se$label.main) 

gbm_SingleR <- SingleR::SingleR(test = gbm@assays$RNA@data, # log2 count data
                                ref = hpca.se, # reference data also as log2 count data
                                labels = hpca.se$label.main #,
                                # clusters = gbm$RNA_snn_res.0.2
)

gbm_SingleR <- SingleR::SingleR(test = Seurat::GetAssayData(gbm, slot = "data"),
                                ref = hpca.se,
                                labels = hpca.se$label.main) 

gbm$SingleR_annot<-gbm_SingleR$labels

Seurat::DimPlot(gbm, group.by = "SingleR_annot", label = T, repel = T)


mean_scores <- tapply(gbm$immune_genes1, gbm$SingleR_annot, mean)
mean_scores[order(mean_scores, decreasing = TRUE)]

# de_genes_res.0.2 <- Seurat::FindAllMarkers(gbm,  min.pct = 0.25)
# de_genes_res.0.2 <- subset(de_genes_res.0.2, de_genes_res.0.2$p_val_adj < 0.05)
# head(de_genes_res.0.2[order(de_genes_res.0.2$p_val_adj),])

de_genes <- read.csv("data/de_genes_gbm_res.0.2.csv", row.names = 1)
de_genes_cl6 <- de_genes[de_genes$cluster == 6,]
de_genes_cl6[immune_genes,]

gbm <- Seurat::SetIdent(gbm, value = "SingleR_annot")

DEG_astro_vs_macro <- Seurat::FindMarkers(gbm,
                                          ident.1 = "Astrocyte", 
                                          ident.2 = "Macrophage",
                                          group.by = gbm$SingleR_annot,
                                          test.use = "wilcox")

top_order <- order(DEG_astro_vs_macro$p_val_adj)
DEG_astro_vs_macro[top_order[1:10],]

png("violinplot_SLC2A5.png", width = 600)
Seurat::VlnPlot(gbm, features = "SLC2A5")
dev.off()

pancreas.integrated <- readRDS("pancreas.integrated.rds")
Seurat::DefaultAssay(pancreas.integrated) <- "RNA"
Seurat::Idents(pancreas.integrated) <- pancreas.integrated$celltype 

Seurat::DimPlot(pancreas.integrated)

pancreas.dg <- subset(pancreas.integrated, idents = c("delta", "gamma"))

table(pancreas.dg$tech, pancreas.dg$celltype)

counts <- GetAssayData(pancreas.dg, slot = "counts")

counts<-counts[rowSums(counts) != 0,]

dge <- edgeR::DGEList(counts = counts)
dge <- edgeR::calcNormFactors(dge) 

design <- model.matrix(~ 0 + celltype + tech, data = pancreas.dg@meta.data)
colnames(design)<-c("delta", "gamma", "celseq2", "fluidigmc1", "smartseq2")

contrast.mat <- limma::makeContrasts(delta - gamma,
                                     levels = design)

# Using limma, fit the linear model, apply the contrasts and calculate
# t-statistics for each gene:
vm <- limma::voom(dge, design = design, plot = TRUE)
fit <- limma::lmFit(vm, design = design)
fit.contrasts <- limma::contrasts.fit(fit, contrast.mat)
fit.contrasts <- limma::eBayes(fit.contrasts)

# Show the top differentially expressed genes:
limma::topTable(fit.contrasts, number = 10, sort.by = "P")

VlnPlot(pancreas.dg, "PPY", split.by = "tech")
VlnPlot(pancreas.dg, "RBP4", split.by = "tech")


BiocManager::install("org.Hs.eg.db")
library(org.Hs.eg.db)
AnnotationDbi::keytypes(org.Hs.eg.db) 

AC_up_DEG <- subset(DEG_astro_vs_macro,
                    DEG_astro_vs_macro$avg_log2FC > 0 & DEG_astro_vs_macro$p_val_adj < 0.05)
AC_up_genes <- rownames(AC_up_DEG)

AC_MAC_GO<-clusterProfiler::enrichGO(AC_up_genes, # vector of up regulated genes
                                     "org.Hs.eg.db", # orgdb= package that contains gene label types correspondances
                                     keyType = "SYMBOL", # indicate that genes are labeled using symbols
                                     ont = "BP", # which of the GO categories to test, here the "Biological Processes"
                                     minGSSize = 50) # exclude gene sets that contain less than 50 genes

enrichplot::emapplot(enrichplot::pairwise_termsim(AC_MAC_GO), 
                     showCategory = 30)


# Other gene set collections can be tested eg the hallmark collection
# Repository of gene set collections: MSigDB

h <- clusterProfiler::read.gmt("data/h.all.v7.2.symbols.xls")

head(h)

AC_MAC_enrich <- clusterProfiler::enricher(gene = AC_up_genes,
                                           universe = rownames(gbm),
                                           pAdjustMethod = "BH",
                                           pvalueCutoff  = 0.05,
                                           qvalueCutoff  = 0.05,
                                           TERM2GENE = h)

View(AC_MAC_enrich@result)

enrichplot::dotplot(AC_MAC_enrich)

myc_target_genes<-subset(h, h$term=="HALLMARK_MYC_TARGETS_V1")$gene

gbm<-Seurat::AddModuleScore(gbm, 
                            features = list(myc_target_genes=myc_target_genes), 
                            name = "myc_target_genes")

head(gbm@meta.data)

# We focused on astrocytes vs macrophages, we see that the score for
# the MYC target gene set is higher in astrocytes than in macrophages
png("violinplot_MYC.png", width = 600)
VlnPlot(gbm, "myc_target_genes1", group.by = "SingleR_annot")
dev.off()
FeaturePlot(gbm, features = "myc_target_genes1")
DimPlot(gbm, group.by = "SingleR_annot")

deng_SCE <- readRDS("data/deng-reads.rds")

class(deng_SCE)

structure(deng_SCE)
# How many mouse cells are at each stage?
table(deng_SCE$cell_type1)
table(deng_SCE$cell_type2)

deng_SCE$cell_type2 <- factor(deng_SCE$cell_type2,
                              levels = c("zy", "early2cell", "mid2cell", "late2cell",
                                         "4cell", "8cell", "16cell", "earlyblast", "midblast",
                                         "lateblast"))



deng_SCE <- BiocSingular::runPCA(deng_SCE, ncomponents = 50)

# Read the Slingshot documentation (?slingshot) and then run Slingshot below.
# Given your understanding of the algorithm and the documentation, what is one
# major set of parameters we omitted here when running Slingshot?
sce <- slingshot::slingshot(deng_SCE, reducedDim = 'PCA')  # no clusters
# Check how the slingshot object looks like
slingshot::SlingshotDataSet(sce)

PCAplot_slingshot <- function(sce, draw_lines = TRUE, variable = NULL, legend = FALSE, ...){
  palf <- colorRampPalette(RColorBrewer::brewer.pal(8, "Set2"))
  paln <- colorRampPalette(RColorBrewer::brewer.pal(9, "Blues"))
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
  plot(pca, bg = colors, pch = 21)
  if(draw_lines){
    slingshot::lines(slingshot::SlingshotDataSet(sce), lwd = 2, ... )
  }
  if(legend & is.factor(variable)){
    legend("bottomright", pt.bg = colpal,legend = levels(variable),pch=21)
    
  }
}

PCAplot_slingshot(sce, variable = sce$slingPseudotime_1)

ggplot(as.data.frame(colData(deng_SCE)), aes(x = sce$slingPseudotime_1, y = cell_type2,
                                             colour = cell_type2)) +
  ggbeeswarm::geom_quasirandom(groupOnX = FALSE) +
  theme_classic() +
  xlab("Slingshot pseudotime") + ylab("Timepoint") +
  ggtitle("Cells ordered by Slingshot pseudotime")

gcdata <- Seurat::CreateSeuratObject(counts = counts(deng_SCE), project = "slingshot")
gcdata <- Seurat::NormalizeData(object = gcdata, normalization.method = "LogNormalize",
                        scale.factor = 10000)
gcdata <- Seurat::FindVariableFeatures(object = gcdata, mean.function = ExpMean, dispersion.function = LogVMR)
topgenes <- Seurat::VariableFeatures(gcdata)

gcdata <- Seurat::ScaleData(object = gcdata, do.center = T, do.scale = F)
gcdata <- Seurat::RunPCA(object = gcdata, pc.genes = gcdata@var.genes, do.print = TRUE, pcs.print = 1:5,
                 genes.print = 5)
gcdata <- Seurat::FindNeighbors(gcdata, reduction="pca", dims = 1:5)

gcdata <- Seurat::FindClusters(object = gcdata,
                       resolution = 0.6)

# Add clustering information from Seurat to the deng_SCE object
# Then run Slingshot using these cluster assignments.

deng_SCE$Seurat_clusters <- as.character(Idents(gcdata))  # go from factor to character
deng_SCE <- slingshot::slingshot(deng_SCE, clusterLabels = 'Seurat_clusters', reducedDim = 'PCA')

head(colData(deng_SCE))


PCAplot_slingshot(deng_SCE, variable = deng_SCE$slingPseudotime_2)

ggplot(as.data.frame(colData(deng_SCE)), aes(x = slingPseudotime_2, y = cell_type2,
                                             colour = cell_type2)) +
  ggbeeswarm::geom_quasirandom(groupOnX = FALSE) +
  theme_classic() +
  xlab("Slingshot pseudotime") + ylab("Timepoint") +
  ggtitle("Cells ordered by Slingshot pseudotime")

PCAplot_slingshot(deng_SCE, variable = deng_SCE$Seurat_clusters, type = 'lineages', col = 'black', legend = TRUE)
PCAplot_slingshot(deng_SCE, variable = deng_SCE$cell_type2, type = 'lineages', col = 'black', legend = TRUE)

### Monocle3 trajectory analysis:

# Monocle3 requires an object of type cell_data_set, so we need to
# create it by extracting the counts and metadata from the Seurat object:

gbm<-readRDS("gbm.rds") # or "data/gbm.rds" for our version

feature_names<-as.data.frame(rownames(gbm))
rownames(feature_names)<-rownames(gbm)
colnames(feature_names)<-"gene_short_name"
library(monocle3)
gbm_monocl<-monocle3::new_cell_data_set(gbm@assays$RNA@counts,
                              cell_metadata = gbm@meta.data,
                              gene_metadata = feature_names)





# We pre-process the newly created object. What does it involve? Check:
?preprocess_cds
gbm_monocl <- monocle3::preprocess_cds(gbm_monocl, num_dim = 50)
monocle3::plot_pc_variance_explained(gbm_monocl)

# Perform UMAP using the implemention in the Monocle3 package and its default parameters:
gbm_monocl<- monocle3::reduce_dimension(gbm_monocl, reduction_method = "UMAP")

# Plot the monocle3 UMAP coloring cells according to the cluster ID ran with Seurat:
monocle3::plot_cells(gbm_monocl, color_cells_by = "RNA_snn_res.0.2")
monocle3::plot_cells(gbm_monocl, genes = "PMP2") # to plot expression level of a gene

# cluster cells using monocle3's clustering function:
gbm_monocl <- monocle3::cluster_cells(gbm_monocl, resolution=0.00025)
p1<-monocle3::plot_cells(gbm_monocl, label_cell_groups = F)
p2<-monocle3::plot_cells(gbm_monocl, color_cells_by = "RNA_snn_res.0.2", label_cell_groups = F)
cowplot::plot_grid(p1, p2, ncol = 2) # Are there differences?

# learn graph (i.e. identify trajectory) using monocle3 umap and clustering:
gbm_monocl<-monocle3::learn_graph(gbm_monocl)
monocle3::plot_cells(gbm_monocl)
monocle3::plot_cells(gbm_monocl, color_cells_by = "RNA_snn_res.0.2")

# replace monocle cluster id by seurat cluster id if we want to keep the same information:
gbm_monocl@clusters$UMAP$clusters<-colData(gbm_monocl)$RNA_snn_res.0.2
names(gbm_monocl@clusters$UMAP$clusters)<-rownames(colData(gbm_monocl))
gbm_monocl<-monocle3::learn_graph(gbm_monocl)
monocle3::plot_cells(gbm_monocl, label_cell_groups = F)

# Select the "initial" cells to calculate pseudotime.
# A pop up window will open and you need to click on the "initial" cells
# (one node per trajectory), then click "Done"
gbm_monocl<-monocle3::order_cells(gbm_monocl)#
monocle3::plot_cells(gbm_monocl,
           color_cells_by = "pseudotime",
           label_cell_groups=F,
           label_leaves=F,
           label_branch_points=FALSE,
           graph_label_size=1.5, cell_size = 1)

# Plot a gene's expression vs pseudotime
monocle3::plot_genes_in_pseudotime(subset(gbm_monocl, rowData(gbm_monocl)$gene_short_name=="PMP2"))

# Record R and package versions!
sessionInfo()

# PLEASE: empty out environment, run garbage collection then restart
rm(list = ls())
gc()
quit()
