pancreas.data <- readRDS(file = "data/pancreas_v3_files/pancreas_expression_matrix.rds")
metadata <- readRDS(file = "data/pancreas_v3_files/pancreas_metadata.rds")

pancreas <- Seurat::CreateSeuratObject(pancreas.data, meta.data = metadata)


pancreas <- Seurat::NormalizeData(pancreas)
pancreas <- Seurat::FindVariableFeatures(pancreas, selection.method = "vst", nfeatures = 2000)
pancreas <- Seurat::ScaleData(pancreas)
pancreas <- Seurat::RunPCA(pancreas, npcs = 30)
pancreas <- Seurat::RunUMAP(pancreas, reduction = "pca", dims = 1:30)

Seurat::DimPlot(pancreas, reduction = "umap", group.by = "tech")
Seurat::DimPlot(pancreas, reduction = "umap", group.by = "celltype", label = TRUE, repel = TRUE)
png("umap_celltype.png", width = 600)
Seurat::DimPlot(pancreas, reduction = "umap", group.by = "celltype")
dev.off()

pancreas.list <- Seurat::SplitObject(pancreas, split.by = "tech")

for (i in 1:length(pancreas.list)) {
  pancreas.list[[i]] <- Seurat::NormalizeData(pancreas.list[[i]], verbose = FALSE)
  pancreas.list[[i]] <- Seurat::FindVariableFeatures(pancreas.list[[i]], selection.method = "vst", nfeatures = 2000, 
                                                     verbose = FALSE)
}

pancreas.anchors <- Seurat::FindIntegrationAnchors(object.list = pancreas.list, dims = 1:30)

pancreas.integrated <- Seurat::IntegrateData(anchorset = pancreas.anchors, dims = 1:30)

Seurat::DefaultAssay(pancreas.integrated) <- "integrated"

pancreas.integrated <- Seurat::ScaleData(pancreas.integrated, verbose = FALSE)
pancreas.integrated <- Seurat::RunPCA(pancreas.integrated, npcs = 30, verbose = FALSE)
pancreas.integrated <- Seurat::RunUMAP(pancreas.integrated, reduction = "pca", dims = 1:30)
Seurat::DimPlot(pancreas.integrated, reduction = "umap", group.by = "tech")
Seurat::DimPlot(pancreas.integrated, reduction = "umap", group.by = "celltype", label = TRUE, repel = TRUE) +
  NoLegend()
plot_grid(p1, p2)

png("umap_integrated_tech.png", width = 600)
Seurat::DimPlot(pancreas.integrated, reduction = "umap", group.by = "tech")
dev.off()
png("umap_integrated_celltype.png", width = 600)
Seurat::DimPlot(pancreas.integrated, reduction = "umap", group.by = "celltype", label = TRUE, repel = TRUE)
dev.off()

saveRDS(pancreas.integrated, "pancreas.integrated.rds")
