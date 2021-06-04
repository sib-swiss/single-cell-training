library(Seurat)

gbm.data <- Seurat::Read10X(data.dir = "data/gbm_dataset/filtered_feature_bc_matrix/")

gbm <- Seurat::CreateSeuratObject(counts = gbm.data,
                                  project = "gbm",
                                  min.cells = 3,
                                  min.features = 100)

gbm <- Seurat::NormalizeData(gbm,
                             normalization.method = "LogNormalize",
                             scale.factor = 10000)


gbm <- FindVariableFeatures(gbm,
                            selection.method = "vst",
                            nfeatures = 2000)

gbm <- Seurat::ScaleData(gbm,
                         features = rownames(gbm))

gbm <- RunPCA(gbm)

bm <- Seurat::RunUMAP(gbm, dims = 1:25)

gbm <- Seurat::FindNeighbors(gbm, dims = 1:25)

gbm <- Seurat::FindClusters(gbm, resolution = 0.2)

library(monocle3)

feature_names <- as.data.frame(rownames(gbm))
rownames(feature_names) <- rownames(gbm)
colnames(feature_names) <- "gene_short_name"
gbm_monocl <- monocle3::new_cell_data_set(gbm@assays$RNA@counts,
                                          cell_metadata = gbm@meta.data,
                                          gene_metadata = feature_names)

gbm_monocl <- monocle3::preprocess_cds(gbm_monocl, num_dim = 50)
monocle3::plot_pc_variance_explained(gbm_monocl)

gbm_monocl <- monocle3::reduce_dimension(gbm_monocl, reduction_method = "UMAP")

plot_cells(gbm_monocl, color_cells_by = "RNA_snn_res.0.2")
plot_cells(gbm_monocl, genes = "PMP2") # to plot expression level of a gene

gbm_monocl <- monocle3::cluster_cells(gbm_monocl, resolution=0.00025)
p1 <- monocle3::plot_cells(gbm_monocl, label_cell_groups = F)
p2 <- monocle3::plot_cells(gbm_monocl, color_cells_by = "RNA_snn_res.0.2", label_cell_groups = F)
cowplot::plot_grid(p1, p2, ncol = 2) # Are there differences?

gbm_monocl <- monocle3::learn_graph(gbm_monocl)
monocle3::plot_cells(gbm_monocl)
monocle3::plot_cells(gbm_monocl, color_cells_by = "RNA_snn_res.0.2")

gbm_monocl@clusters$UMAP$clusters <- colData(gbm_monocl)$RNA_snn_res.0.2
names(gbm_monocl@clusters$UMAP$clusters) <-rownames(colData(gbm_monocl))
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
