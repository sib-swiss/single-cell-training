library(Seurat)
library(stringr)
library(ggplot2)
library(Matrix)
library(dittoSeq)

# download from:
# https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE202212&format=file


h5_files <- list.files("data/project1/GSE202212_RAW/", full.names = TRUE)

seu_list <- lapply(h5_files,
  function(x) {
    sample_id <- basename(x) |> word(2, sep = "_")
    
    Read10X_h5(x, use.names = TRUE, unique.features = TRUE) |> 
      CreateSeuratObject(project = sample_id) 
  }
)

names(seu_list) <- basename(h5_files) |> word(2, sep = "_")

seu <- merge(seu_list[[1]], y = seu_list[2:length(seu_list)],
             add.cell.ids = names(seu_list),
             project = "zebra") |>
  JoinLayers()

Seurat::FeatureScatter(seu, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

Seurat::VlnPlot(seu, features = c("nCount_RNA",
                                  "nFeature_RNA"))

genes <- rownames(seu)
genes[grepl("^hb[^(p)]", genes)]

# mitochondrial geseu# mitochondrial genes
  seu <- Seurat::PercentageFeatureSet(seu, 
                                      pattern = "^mt-", 
                                      col.name = "percent.mito")

# ribosomal genes
seu <- Seurat::PercentageFeatureSet(seu, 
                                    pattern = "^rp[sl]",
                                    col.name = "percent.ribo")

seu <- Seurat::PercentageFeatureSet(seu,
                                    pattern = "^hb[^(p)]",
                                    col.name = "percent.globin")


Seurat::VlnPlot(seu, features = c("percent.mito",
                                  "percent.ribo",
                                  "percent.globin"))

Seurat::FeatureScatter(seu, 
                       feature1 = "percent.globin", 
                       feature2 = "percent.ribo")

most_expressed_boxplot <- function(object, ngenes = 20){
  
  # matrix of raw counts
  cts <- Seurat::GetAssayData(object, assay = "RNA", layer = "counts")
  
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
                percent.mito < 25)

Seurat::VlnPlot(seu, features = c("nFeature_RNA",
                                  "percent.mito"))


# Don't run it yet! Read the exercise first
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

seu <- Seurat::ScaleData(seu)

seu <- Seurat::RunPCA(seu)

Seurat::DimPlot(seu, reduction = "pca")

Seurat::FeaturePlot(seu, reduction = "pca", features = "percent.ribo")

Seurat::DimHeatmap(seu, dims = 1:12, cells = 500, balanced = TRUE)

Seurat::ElbowPlot(seu, ndims = 40)

seu <- Seurat::RunUMAP(seu, dims = 1:25)

Seurat::DimPlot(seu, reduction = "umap", group.by = "orig.ident")

Seurat::FeaturePlot(seu, "islr2")
Seurat::FeaturePlot(seu, "GFPx")

mg_cell_markers <- c("fabp7a", "rlbp1a", "slc1a2b")

Seurat::FeaturePlot(seu, mg_cell_markers[3])

Seurat::FeaturePlot(seu, "ppdpfa")
Seurat::FeaturePlot(seu, "pde6a")

Seurat::FeaturePlot(seu, "meig1")

Seurat::FeatureScatter(seu, "ppdpfa", "pde6a")

seu <- Seurat::FindNeighbors(seu, dims = 1:25, reduction = "pca")

seu <- Seurat::FindClusters(seu, resolution = seq(0.1, 0.8, by=0.1))

Seurat::DimPlot(seu, group.by = "RNA_snn_res.0.8")

dittoDimPlot(seu, var = "RNA_snn_res.0.7")

Seurat::VlnPlot(seu, features = "GFPx", group.by = "RNA_snn_res.0.3")

Seurat::VlnPlot(seu, features = "pde6a", group.by = "RNA_snn_res.0.7")

