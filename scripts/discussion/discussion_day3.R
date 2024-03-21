# --------  Day 3

library(Seurat)
# BiocManager::install("edgeR")
library(edgeR)
library(limma)
library(dittoSeq)
library(dplyr)

### Find marker genes ("differential" gene expression)
# setwd("/export/scratch/twyss/SIB_scRNAseq_course/July2023/data/")
setwd("/home/rstudio/single_cell_course/")

seu <- readRDS("seu_day2_part2.rds")

# Use wilcoxon test from Seurat (make sure the default cell identity
# is set to what you want it to be)
Idents(seu) # res 0.3 with clusters 0-12
# make sure the Default assay is set to RNA:
?FindAllMarkers # check Value section for type of output

# this takes a while! Already run:
de_genes <- Seurat::FindAllMarkers(seu,  min.pct = 0.25,
                                    only.pos = TRUE)
head(de_genes) 

# sometimes p>0.05 are included, remove non-significant genes
range(de_genes$p_val_adj) # [1] 0 1
de_genes<-subset(de_genes, de_genes$p_val_adj<0.05)

# write to csv file if you want to store (and to avoid loosing time recalculating)
write.csv(de_genes, "de_genes_FindAllMarkers.csv", row.names = F, quote = F)

# create dotplot with top 3 markers per cluster:
library(dplyr) 
# use of %>% pipe of magrittr package
# What the function does is to pass the left hand side of the operator 
# to the first argument of the right hand side of the operator. 
top_specific_markers <- de_genes %>% 
  dplyr::group_by(cluster) %>%
  dplyr::top_n(3, avg_log2FC) # select top 3 rows by value

View(top_specific_markers)
dittoSeq::dittoDotPlot(seu, vars = unique(top_specific_markers$gene), 
                       group.by = "RNA_snn_res.0.3")
dittoSeq::dittoDotPlot(seu, vars = unique(top_specific_markers$gene), 
                       group.by = "SingleR_annot")

# check if T cell genes are within DE genes
tcell_genes <- c("IL7R", "LTB", "TRAC", "CD3D")
de_genes[de_genes$gene %in% tcell_genes,]

# marker detection between CD4+ and CD8+ T cells
seu <- Seurat::SetIdent(seu, value = "SingleR_annot")

deg_cd8_cd4 <- Seurat::FindMarkers(seu,
                                   ident.1 = "CD8+ T cells",
                                   ident.2 = "CD4+ T cells",
                                   group.by = seu$SingleR_annot,
                                   test.use = "wilcox")
deg_cd8_cd4<-subset(deg_cd8_cd4, deg_cd8_cd4$p_val_adj<0.05)

View(deg_cd8_cd4)
# to have info on the output, see "Value" section in help:
?FindMarkers
deg_cd8_cd4[c("CD4", "CD8A", "CD8B"),]

# How many of the DE genes are part of the variable features?
length(which(rownames(deg_cd8_cd4) %in% VariableFeatures(seu)))
# 54

# Plotting the T cell genes only in T cells:
Seurat::VlnPlot(seu, 
                features = c("CD4", "CD8A", "CD8B"),
                idents = c("CD8+ T cells", "CD4+ T cells"),
                group.by = "SingleR_annot")
# note on surface markers, eg CD4, which is often low at the mRNA level
# CITE-seq for surface marker identification at the protein level

# Pseudo-bulk DGE analysis with limma, summing counts per cell:
##New script for generating pseudobulk
##This script follows the vignette on this page 
# http://bioconductor.org/books/3.17/OSCA.multisample/multi-sample-comparisons.html
# taking the proB data 
proB <- readRDS("course_data/proB.rds")

DimPlot(proB, group.by = "orig.ident")

table(proB@meta.data$type)
# ETV6-RUNX1      PBMMC 
#      2000       1021

DimPlot(proB, group.by = "type")

head(proB@meta.data)

#taking the proB data 
Seurat::DefaultAssay(proB) <- "RNA"
Seurat::Idents(proB) <- proB$orig.ident

## add the patient id also for paired DGE
proB$patient.id<-gsub("ETV6-RUNX1", "ETV6_RUNX1", proB$orig.ident)
proB$patient.id<-sapply(strsplit(proB$patient.id, "-"), '[', 2)

## Here we do perform pseudo-bulk:
##first a mandatory column of sample needs to be added to the meta data that is the grouping factor, should be the samples
proB$sample <- factor(proB$orig.ident)

# aggergate the cells per sampple
bulk <- Seurat::AggregateExpression(proB, group.by = "sample",
                                    return.seurat = TRUE,
                                    assay = "RNA")
bulk@meta.data

# create a metadata data frame based on the aggregated cells
meta_data <- unique(proB@meta.data[, c("orig.ident",
                                       "sample", "type",
                                       "patient.id")])
rownames(meta_data) <- meta_data$orig.ident
bulk@meta.data <- meta_data[colnames(bulk), ]

##have a look at the counts
counts <- Seurat::GetAssayData(bulk, layer = "counts") |> as.matrix()

head(counts)

#have a look at the colData of our new object summed, can you see type and 
#patient.id are there
head(bulk@meta.data)

#As in the standard limma analysis generate a DGE object
library(edgeR)
library(limma)

#As in the standard limma analysis generate a DGE object

y <- edgeR::DGEList(counts, samples = bulk@meta.data)

##filter lowly expressed (recommanded for limma)
keep <- edgeR::filterByExpr(y, group = bulk$type)
y <- y[keep,]

##see how many genes were kept 
summary(keep)
#    Mode   FALSE    TRUE 
# logical   11086   10017 

## Create the design matrix and include the patient ID (or scRNAseq technology, etc) as a covariate:
design <- model.matrix(~0 + y$samples$type + y$samples$patient.id)

# Have a look
design

# change column/rownames names to more simple group names: 
colnames(design) <- make.names(c("ETV6-RUNX1", "PBMMC","patient2","patient3"))
rownames(design)<-rownames(bulk@meta.data)

# Create contrasts, i.e. specify which groups we want to compare, here we want
# to find genes differentially expressed between cluster 1 and cluster 2.
contrast.mat <- limma::makeContrasts(ETV6.RUNX1 - PBMMC,
                                     levels = design)

dge <- edgeR::calcNormFactors(y)  

# Run limma
# What is voom doing?
# Counts are transformed to log2 counts per million reads (CPM), where “per million reads” is defined based on the 
# normalization factors we calculated with calcNormFactors().
# A linear model is fitted to the log2 CPM for each gene, and the residuals are calculated.
# A smoothed curve is fitted to the sqrt(residual standard deviation) by average expression (see red line in plot mean-variance trend plot).
# The smoothed curve is used to obtain weights for each gene and sample that are passed into limma along with the log2 CPMs.

?voom
vm <- limma::voom(dge, design = design, plot = TRUE)
?lmFit
fit <- limma::lmFit(vm, design = design)
fit.contrasts <- limma::contrasts.fit(fit, contrast.mat)
?eBayes
fit.contrasts <- limma::eBayes(fit.contrasts)

# Show the top differentially expressed genes:
limma::topTable(fit.contrasts, number = 10, sort.by = "P")
limma_de <- limma::topTable(fit.contrasts, number = Inf, sort.by = "P")
length(which(limma_de$adj.P.Val<0.05))

# Already run: 
tum_vs_norm <- Seurat::FindMarkers(proB, 
                                   ident.1 = "ETV6-RUNX1", 
                                   ident.2 = "PBMMC", 
                                   group.by = "type")
tum_vs_norm <- subset(tum_vs_norm, tum_vs_norm$p_val_adj<0.05)

merge_limma_FindMarkers <- merge(tum_vs_norm, limma_de, by="row.names",
                                 all.x=T)

par(mar=c(4,4,4,4))
plot(merge_limma_FindMarkers$avg_log2FC,
     merge_limma_FindMarkers$logFC,
     xlab="log2FC Wilcoxon", ylab="log2FC limma",
     pch=15, cex=0.5)
abline(a=0, b=1, col="red")


# ---- Enrichment analysis
library(clusterProfiler)
library(enrichplot)
# Need a gene annotation package for Human:
# BiocManager::install("org.Hs.eg.db", update = FALSE)
library(org.Hs.eg.db)
# check gene label types allowed, we can use SYMBOL:
AnnotationDbi::keytypes(org.Hs.eg.db)

# select genes down-regulated in tumor:
tum_down  <- subset(limma_de,
                    limma_de$logFC < -1 
                    & limma_de$adj.P.Val <  0.05)
tum_down_genes <- rownames(tum_down)
length(tum_down_genes) # 958

# over-representation analysis (Fisher test) for down-reg genes:
?enrichGO
# Takes a while, already run:
tum_vs_norm_go <- clusterProfiler::enrichGO(gene = tum_down_genes,
                                            OrgDb =  "org.Hs.eg.db",
                                            keyType = "SYMBOL",
                                            ont = "BP",
                                            universe = rownames(limma_de), 
                                            minGSSize = 50)
View(tum_vs_norm_go@result)
class(tum_vs_norm_go@geneSets)
head(names(tum_vs_norm_go@geneSets))
# Check for key word in Description of GO terms:
tum_vs_norm_go@result[grep("cell cycle", tum_vs_norm_go@result$Description), c(1,2,6)]

# remove redundant gene sets (already run):
enr_go <- clusterProfiler::simplify(tum_vs_norm_go)
View(enr_go@result)

# check if significant gene sets share overlapping genes:
enrichplot::emapplot(enrichplot::pairwise_termsim(enr_go),
                     showCategory = 30, cex_label_category = 0.5)

# Test for hallmark gene set over-representation:
?msigdbr
gmt <- msigdbr::msigdbr(species = "Homo sapiens", category = "H")
?clusterProfiler::read.gmt

tum_vs_norm_enrich <- clusterProfiler::enricher(gene = tum_down_genes,
                                                universe = rownames(limma_de),
                                                pAdjustMethod = "BH",
                                                TERM2GENE = gmt[,c("gs_name", "gene_symbol")])
View(tum_vs_norm_enrich@result[which(tum_vs_norm_enrich@result$p.adjust<0.05),])

# Bonus code:
# GSEA example using gseGO(), with t-statistic of limma output:
# create ranked list of t-statistics:
gene.list<-limma_de$t
names(gene.list)<-rownames(limma_de)
gene.list<-sort(gene.list, decreasing = T)
head(gene.list)
#      RPS4Y2       SDC2       CTGF AP005530.2      GNG11   HLA-DQA1 
#    15.39674   15.33434   14.89301   14.27130   13.70183   13.39202 

set.seed(1234)
# Already run:
ego_1<-gseGO(geneList = gene.list,
           ont="BP",
           OrgDb = "org.Hs.eg.db",
           keyType = "SYMBOL",
           minGSSize = 60,
           eps=1e-60,
           seed=T)
ego <- clusterProfiler::simplify(ego_1)
head(ego@result[,c(2:7)])
#                                 Description setSize enrichmentScore       NES       pvalue     p.adjust
# GO:0000280                 nuclear division     292      -0.5139084 -2.583266 3.544230e-23 4.377997e-20
# GO:0007059           chromosome segregation     267      -0.5233834 -2.608855 6.808704e-23 4.377997e-20
# GO:0140014         mitotic nuclear division     225      -0.5404579 -2.661770 2.584499e-21 8.309163e-19
# GO:0098813   nuclear chromosome segregation     208      -0.5369762 -2.606928 2.340434e-19 4.299712e-17
# GO:0051301                    cell division     420      -0.4367175 -2.300783 3.118075e-19 5.012305e-17
# GO:0010564 regulation of cell cycle process     443      -0.4125029 -2.183002 2.049784e-17 2.928913e-15

# barcode plot of the top GO term:
gseaplot(ego, geneSetID = "GO:0000280", title="GO:0000280: nuclear division")

# Bonus code over-representation of KEGG collection
# Needs convertion of gene symbols to only 3 allowed keyType, eg we will use ncbi-geneid=entrezID
?enrichKEGG
# Use bitr to convert the gene symbols in tum_down_genes to entrezID
tum_down_genes

keytypes(org.Hs.eg.db)
# convert from= "ENSEMBL" to "SYMBOL" and "ENTREZID"
gene_convert <- clusterProfiler::bitr(as.character(tum_down_genes), 
                     fromType="SYMBOL", 
                     toType=c("ENTREZID"), OrgDb="org.Hs.eg.db")
gene_convert_universe<-clusterProfiler::bitr(as.character(rownames(proB)), 
                                             fromType="SYMBOL", 
                                             toType=c("ENTREZID"), OrgDb="org.Hs.eg.db")
# Check the format of the data frame obtained after conversion:                     
head(gene_convert)
dim(gene_convert)

tum_vs_norm_kegg<-clusterProfiler::enrichKEGG(gene = gene_convert$ENTREZID,
                                              organism = "hsa",
                                              keyType = "ncbi-geneid",
                                              universe = gene_convert_universe$ENTREZID)
View(tum_vs_norm_kegg@result)


##### Trajectory analysis
library(Seurat)
library(SingleCellExperiment)
library(scater)
library(slingshot)
library(ggplot2)
library(ggbeeswarm)
library(batchelor)
# sudo apt-get update
# sudo apt-get install libcairo2-dev libxt-dev
# devtools::install_github('cole-trapnell-lab/monocle3')
# devtools::install_github('cole-trapnell-lab/monocle3', ref="develop")
library(monocle3)

# ---- Trajectory with slingshot:

# Download the dataset from github within Terminal tab:
#  cd course_data/
#  wget https://github.com/hemberg-lab/nrg-paper-figures/blob/master/deng-reads.rds?raw=true
# mv deng-reads.rds\?raw\=true deng-reads.rds

# Import SingleCellExperiment object:
deng_SCE <- readRDS("course_data/deng-reads.rds")
class(deng_SCE)

# Change levels as developmental order:
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

# Run PCA using scater:
deng_SCE <- scater::runPCA(deng_SCE, ncomponents = 50)

# Use the reducedDim function to access the PCA and store the results.
pca <- SingleCellExperiment::reducedDim(deng_SCE, "PCA")

# Describe how the PCA is stored in a matrix. Why does it have this structure?
head(pca)

# Add 2 first principal components to SCE object:
deng_SCE$PC1 <- pca[, 1]
deng_SCE$PC2 <- pca[, 2]

# To view the meta data as a data frame:
View(as.data.frame(colData(deng_SCE)))

# PCA plot with cells colored according to developmental stage:
ggplot(as.data.frame(colData(deng_SCE)), aes(x = PC1, y = PC2, color = cell_type2)) +
  geom_point(size=2, shape=20) +
  theme_classic() +
  xlab("PC1") + ylab("PC2") + ggtitle("PC biplot")

# Plot PC1 vs cell_type2.
deng_SCE$pseudotime_PC1 <- base::rank(deng_SCE$PC1)  # rank cells by their PC1 score

# Create a jitter plot of developmental stage against PC1:
ggplot(as.data.frame(colData(deng_SCE)), aes(x = pseudotime_PC1, y = cell_type2,
                                             colour = cell_type2)) +
  ggbeeswarm::geom_quasirandom(groupOnX = FALSE) +
  theme_classic() +
  xlab("PC1") + ylab("Timepoint") +
  ggtitle("Cells ordered by first principal component")

# Run Slingshot 
?slingshot::slingshot
sce <- slingshot::slingshot(deng_SCE, reducedDim = 'PCA')
# Check how many curves were determined:
SlingshotDataSet(sce)

# custom function to plot the PCA based on a slingshot object, to color the cells
# either according to a factor or to a numeric variable, with and without the slingshot
# lines
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

# Plot PC1 versus PC2 with slingshot pseudotime added as line and coloring the cells:
PCAplot_slingshot(sce, variable = sce$slingPseudotime_1, draw_lines = TRUE)

# Plot developmental stage against pseudotime:
  ggplot(as.data.frame(colData(deng_SCE)), aes(x = sce$slingPseudotime_1,
                                               y = cell_type2,
                                               colour = cell_type2)) +
  ggbeeswarm::geom_quasirandom(groupOnX = FALSE) +
  theme_classic() +
  xlab("Slingshot pseudotime") + ylab("Timepoint") +
  ggtitle("Cells ordered by Slingshot pseudotime")

  
# Slingshot will use the centroid of the clusters   
# Generate unsupervised clustering of the cells using seurat: 
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
View(gcdata@meta.data) # zygote = cluster 2
  
# Now we can add these clusters to the SCE object and re-run slingshot providing cluster labels :
deng_SCE$Seurat_clusters <- as.character(Idents(gcdata))  # go from factor to character

sce <- slingshot::slingshot(deng_SCE,
                              clusterLabels = 'Seurat_clusters',
                              reducedDim = 'PCA',
                              start.clus = "2") # zygote
  
# Check how the slingshot object has evolved, it now contains 2 curves:
SlingshotDataSet(sce)
  
# Plot PC1 versus PC2 colored by slingshot pseudotime:
PCAplot_slingshot(sce, variable = sce$slingPseudotime_2)
  
# Plot Slingshot pseudotime 1 or 2 vs cell stage.
ggplot(data.frame(cell_type2 = deng_SCE$cell_type2,
                    slingPseudotime_1 = sce$slingPseudotime_1),
         aes(x = slingPseudotime_1, y = cell_type2,
             colour = cell_type2)) +
    ggbeeswarm::geom_quasirandom(groupOnX = FALSE) +
    theme_classic() +
    xlab("Slingshot pseudotime") + ylab("Timepoint") +
    ggtitle("Cells ordered by Slingshot pseudotime 1")
  
ggplot(data.frame(cell_type2 = deng_SCE$cell_type2,
                    slingPseudotime_2 = sce$slingPseudotime_2),
         aes(x = slingPseudotime_2, y = cell_type2,
             colour = cell_type2)) +
    ggbeeswarm::geom_quasirandom(groupOnX = FALSE) +
    theme_classic() +
    xlab("Slingshot pseudotime") + ylab("Timepoint") +
    ggtitle("Cells ordered by Slingshot pseudotime 2")
  
#  Particularly the later stages, separation seems to improve. Since we have included the Seurat clustering, we can plot the PCA, with colors according to these clusters:
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

# Slingshot analysis with end clusters indicated:
sce <- slingshot::slingshot(deng_SCE,
                            clusterLabels = 'Seurat_clusters',
                            reducedDim = 'PCA',
                            end.clus = c("0", "3", "5")) ## late blastula clusters, 
                              # check which would be the best cluster according to biology knowledge

# PCA plot with curves ending at 3 clusters:
PCAplot_slingshot(sce,
                  variable = deng_SCE$Seurat_clusters,
                  type = 'lineages',
                  col = 'black',
                  legend = F)



# ---- Trajectory with Monocle3 
# Generate a monocle3 object (with class cell_data_set) from our Seurat object:

# get matrix and filter for minimum number of cells and 
# features (the latter is a fix for backward compatibility)
# create gene metadata data.frame
feature_names <- as.data.frame(rownames(seu))
rownames(feature_names) <- rownames(seu)
colnames(feature_names) <- "gene_short_name"

# initiate monocle object from seurat count table 
seu_monocl <- monocle3::new_cell_data_set(Seurat::GetAssayData(seu,
                                                               layer = "counts"),
                                          cell_metadata = seu@meta.data,
                                          gene_metadata = feature_names)

# # Preprocess the dataset:
?preprocess_cds
# already ran: (takes a while)
seu_monocl <- monocle3::preprocess_cds(seu_monocl)

# elbow plot:
monocle3::plot_pc_variance_explained(seu_monocl)

# Perform UMAP using the implementation in the monocle3 package and its default parameters:
seu_monocl <- monocle3::reduce_dimension(seu_monocl, reduction_method = "UMAP")

# Plot the monocle3 UMAP coloring cells according to the cluster ID ran with Seurat:
monocle3::plot_cells(seu_monocl, 
                     color_cells_by = "RNA_snn_res.0.3", 
                     cell_size = 1, 
                     show_trajectory_graph = FALSE)
# Slightly different UMAP between Seurat and Monocle3:
p1<-DimPlot(seu, group.by = "RNA_snn_res.0.3", label = T)
p2<-monocle3::plot_cells(seu_monocl, 
                         color_cells_by = "RNA_snn_res.0.3", 
                         cell_size = 0.7, 
                         show_trajectory_graph = FALSE,
                         label_cell_groups = T, 
                         group_label_size=6)
cowplot::plot_grid(p1,p2, ncol = 2)

# plot B cell marker:
monocle3::plot_cells(seu_monocl, genes = "CD79A", 
                     show_trajectory_graph = FALSE, 
                     cell_size = 1, norm_method = "log")
# Cluster cells using monocle3‘s clustering function:
?cluster_cells
seu_monocl <- monocle3::cluster_cells(seu_monocl, resolution=0.00025)
# monocle3::partitions(seu_monocl)
monocle3::plot_cells(seu_monocl, label_cell_groups = F)


# learn graph (i.e. identify trajectory) using monocle3 UMAP and clustering:
seu_monocl <- monocle3::learn_graph(seu_monocl)
monocle3::plot_cells(seu_monocl)

# Find the CD34+ B-cell cluster in the monocle UMAP. This cluster has a high expression of CD79A and expresses CD34.
monocle3::plot_cells(seu_monocl, genes = c("CD79A", "CD34"),
                     show_trajectory_graph = FALSE, 
                     cell_size = 0.7, group_label_size = 4)

# Select the “initial” cells in the B-cell cluster to calculate pseudotime. 
# The initial cells in this case are the CD34+ B-cells we have just identified. 
# A pop up window will open and you need to click on the “initial” cells (one node per trajectory), then click “Done”.
# seu_monocl <- monocle3::learn_graph(seu_monocl)

seu_monocl<-monocle3::order_cells(seu_monocl)

monocle3::plot_cells(seu_monocl,
                     color_cells_by = "pseudotime",
                     label_cell_groups=F,
                     label_leaves=F,
                     label_branch_points=FALSE,
                     graph_label_size=1.5, cell_size = 1)

# In order to find genes which expression is affected by pseudotime, 
# we first have to isolate the B-cell cluster. Therefore, extract 
# all cells in the B-cell cluster with the interactive choose_cells function:
seuB <- choose_cells(seu_monocl)
class(seuB)

# Check whether you have selected the right cells:
plot_cells(seuB, show_trajectory_graph = T, cell_size = 1)
plot_cells(seuB)

# Now we can use the cells in this trajectory to test which genes are affected by the trajectory:
?graph_test # Moran's I test, to test for spatial autocorrelation
# already ran:
pr_test <- graph_test(seuB, 
                      cores=4, 
                      neighbor_graph = "principal_graph")
# order by test statistic
pr_test <- pr_test[order(pr_test$morans_test_statistic, 
                         decreasing = TRUE),]
View(pr_test)
# There are some interesting genes in there, for example related to cell cycling (MKI67, CKS2), 
# related to B-cell development 
# (CD34, MS4A1) and immunoglobulins (IGLL1 and IGLL5). We can plot those in the UMAP:
goi <- c("CD34", "MS4A1", "IGLL1", "IGLL5", 
          "MKI67", "CKS2")
plot_cells(seuB, label_cell_groups=FALSE, genes = goi,
           show_trajectory_graph=FALSE, cell_size = 1)

# But also against pseudotime:
?monocle3::clusters
table(clusters(seuB))
# 1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16 
# 0 770   0   0   0   0 464   0   0   0   0   0   0   0   0   0 
seuB@colData$monocle_cluster <- clusters(seuB)
plot_genes_in_pseudotime(subset(seuB, 
                                rowData(seuB)$gene_short_name %in% goi),
                         min_expr=0.5, color_cells_by = "monocle_cluster")
