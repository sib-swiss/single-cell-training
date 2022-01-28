## Learning outcomes

**After having completed this chapter you will be able to:**

- Calculate and visualize quality measures based on:
    - mitochondrial genes
    - ribosomal genes
    - hemoglobin genes
    - relative gene expression
- Interpret the above quality measures per cell.
- Perform cell filtering based on quality thresholds

## Material

[:fontawesome-solid-file-pdf: Download the presentation](../assets/pdf/quality_control.pdf){: .md-button }

## Exercises

### Visualizing QC per cell and gene

While generating the `Seurat` object, there were already some quality measures calculated for each cell, namely the total UMI counts per cell (`nCount_RNA`) and the total number of detected features per cell (`nFeature_RNA`). We can plot those in a violin plot and evaluate their distribution per sample:

```R
Seurat::VlnPlot(seu, features = c("nCount_RNA",
                                  "nFeature_RNA"))
```

You can see that there is quite a wide range for both. A cell with low number of detected features or counts might not give you a lot of information, while a high number of detected features/counts might point to doublets.

Single cells have often undergone sampling and/or dissociation and/or sorting. Therefore, there are often cells and genes in your dataset that cause variation due to technical reasons. In the following steps, we will visualize those and make decisions on whether or not to remove cells or genes with extreme values.

We will start with calculating the percentage of counts coming from transcript types:

- **Mitochondrial genes**: If a cell membrane is damaged, it looses free RNA quicker compared to mitochondrial RNA, because the latter is part of the mitochondrion. A high relative amount of mitochondrial counts can therefore point to damaged cells ([Lun et al. 2016](https://f1000research.com/articles/5-2122)). 
- **Ribosomal genes**: Although rRNA transcripts are depleted due to poly-A enrichment, they are often abundant. They do not point to specific issues, but it can be good to have a look at their relative abundance. Sometimes, they might even have a biological relevance (e.g. [Caron et al. 2020](https://www.nature.com/articles/s41598-020-64929-x)).
- **Hemoglobin genes**: these transcripts are very abundant in erythrocytes. Depending on your application, you can expect 'contamination' of erythrocytes and select against it. 

In order to have an idea about the relative counts of these type of genes in our dataset, we can calculate their expression as relative counts in each cell. We do that by selecting genes based on patterns (e.g. `^MT-` matches with all gene names starting with `MT`, i.e. mitochondrial genes):

```R
# mitochondrial genes
seu <- Seurat::PercentageFeatureSet(seu, 
                                    pattern = "^MT-", 
                                    col.name = "percent.mito")

# ribosomal genes
seu <- Seurat::PercentageFeatureSet(seu, 
                                    pattern = "^RP[SL]",
                                    col.name = "percent.ribo")

# hemoglobin genes (but not HBP)
seu <- Seurat::PercentageFeatureSet(seu,
                                    pattern = "^HB[^(P)]",
                                    col.name = "percent.globin")
```

**Exercise:** Run the commands and check out the metadata data frame at `sc@meta.data`. What has changed?

??? done "Answer"
    If we type `head(sc@meta.data)` it returns:

    ```
                               orig.ident nCount_RNA nFeature_RNA percent.mito percent.ribo percent.globin
    PBMMC-1_AAACCTGCAGACGCAA-1    PBMMC-1       2401          909     2.540608     28.65473      0.1665973
    PBMMC-1_AAACCTGTCATCACCC-1    PBMMC-1       3532          760     5.181200     55.03964      0.1981880
    PBMMC-1_AAAGATGCATAAAGGT-1    PBMMC-1       3972         1215     4.934542     30.43807      0.3776435
    PBMMC-1_AAAGCAAAGCAGCGTA-1    PBMMC-1       3569          894     3.250210     55.02942      0.3642477
    PBMMC-1_AAAGCAACAATAACGA-1    PBMMC-1       2982          730     3.688799     54.49363      0.1006036
    PBMMC-1_AAAGCAACATCAGTCA-1    PBMMC-1      22284         3108     3.181655     23.40693     36.9682283
    ```
    So, the function `PercentageFeatureSet` adds a column to `meta.data`, specifying the percentage of counts for the specified gene sets. 


Now we can plot the distribution of these percentages in a violin plot:

```R
Seurat::VlnPlot(seu, features = c("percent.mito",
                                  "percent.ribo",
                                  "percent.globin"))
```

You can see that `PBMMC-2` is quite different from the two others, it has a group of cells with very low ribosomal counts and very high globin counts. Maybe these two percentages are negatively correlated? Let's have a look, by plotting the two percentages against each other:

```R
Seurat::FeatureScatter(seu, 
                       feature1 = "percent.globin", 
                       feature2 = "percent.ribo")
```

**Exercise:** Are they correlated? What kind of cells might have a high abundance of hemoglobin transcripts and low ribosomal transcripts? 

??? done "Answer"
    Yes there is a negative correlation:

    <figure>
    <img src="../../assets/images/globin_vs_ribo.png" width="500"/>
    </figure>

    Erythrocytes (red blood cells) have a high abundance of hemoglobin transcripts and low abundance of ribosomal transcripts. As they don't have a nucleus, we don't expect them in this set of Bone Marrow **Mononuclear** Cells (BMMCs). 

We can also evaluate the relative expression of other genes in our dataset, for example, the ones that are most highly expressed. Some very highly expressed genes might point to a technical cause, and we might consider to remove them. Below you will find a simple function to generate a boxplot of relative counts per gene per cell. Load it into your environment and run it on our `seu` object:

```R
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
```

As for most 10X based poly-A enriched single cell datasets, we find a relatively high expression of MALAT1. Many researchers choose to remove it, but it can have biological relevance (e.g. [Shaat et al. 2021](https://www.nature.com/articles/s41420-020-00383-y)). 

### Cell filtering

Based on the QC process we went through we can come to the following conclusions:

- There no cells with very high mitochondrial gene counts.
- There are some cells with a hemoglobin and low ribosomal counts, and these are probably erythrocytes.
- There are some cells with a very low and very high number of features. These might point to non-informative cells and doublets respectively. 
- The 'usual suspect' MALAT1 sometimes makes up a large part of the counts per cell. As we do not see any other suggestions of dying/stressed cells, we leave it in. 

In this case, a sensible decision would be to do mild filtering on the number of features per cell and mitochondrial counts. We can leave the possible erythrocytes in for now, and see where they end up later during the dimensionality reduction. 

In the M&M of the [publication](https://www.nature.com/articles/s41598-020-64929-x#Sec7), the authors describe that they have used a threshold of < 8% mitochondrial counts and > 200 features per cell. To filter against possible doublets, here, we also filter out cells with > 5000 detected features/cell. Filtering `Seurat` objects can be done with the `subset` method for class `SeuratObject`:

```R
seu <- subset(seu, subset = nFeature_RNA > 200 & 
                nFeature_RNA < 5000 &
                percent.mito < 8)
```

To evaluate this did the trick we can visualize those parameters again in a violin plot:

```R
VlnPlot(seu, features = c("nFeature_RNA",
                          "percent.mito"))
```

### Save the dataset and clear environment

Now, save the dataset so you can use it tomorrow:

```R
saveRDS(seu, "seu_day1.rds")
```

Clear your environment:

```R
rm(list = ls())
gc()
.rs.restartR()
```
