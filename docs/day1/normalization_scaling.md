## Learning outcomes

**After having completed this chapter you will be able to:**

- Describe and perform standard procedures for normalization and scaling with the package `Seurat`
- Select the most variable genes from a `Seurat` object for downstream analyses

## Material

- [Seurat vignette](https://satijalab.org/seurat/articles/pbmc3k_tutorial.html)

## Exercises

> :material-zodiac-cancer: This chapter uses the `gbm` dataset

### Normalization

!!! bug
    Recommendations are unclear here. Should we use the `NormalizeData` -> `FindVariableFeatures` -> `ScaleData` workflow or `SCTransform`??

After removing unwanted cells from the dataset, the next step is to normalize the data.
By default, Seurat employs a global-scaling normalization method `"LogNormalize"` that normalizes the feature expression measurements for each cell by the total expression, multiplies this by a scale factor (10,000 by default), and log-transforms the result.
Normalized values are stored in the "RNA" slot of the gbm object.

**Exercise:** Have a look at the assay data before and after running `NormalizeData()`. Did it change?

!!! hint
    You can extract assay data with the function `Seurat::GetAssay`

??? done "Answer"
    You can check out some assay data with:

    ```R
    Seurat::GetAssay(gbm)[1:10,1:10]  
    ```
    Returning:

    === "Before normalization"

        ```
        10 x 10 sparse Matrix of class "dgCMatrix"
          [[ suppressing 10 column names ‘AAACCCAAGGCGATAC-1’, ‘AAACCCACAAGTCCCG-1’, ‘AAACCCACAGATGCGA-1’ ... ]]

        AL627309.1 . . . . . . . . . .
        AL627309.5 . . . . . . . . . .
        AP006222.2 . . . . . . . . . .
        LINC01409  . . . 1 . . . . . .
        FAM87B     . . . . . . . . . .
        LINC01128  . . 1 . . . . 1 . 1
        LINC00115  . . . . . . . . . .
        FAM41C     . . . . . . . . . .
        AL645608.6 . . . . . . . . . .
        AL645608.2 . . . . . . . . . .
        ```

    === "After normalization"

        ```
        10 x 10 sparse Matrix of class "dgCMatrix"
           [[ suppressing 10 column names ‘AAACCCAAGGCGATAC-1’, ‘AAACCCACAAGTCCCG-1’, ‘AAACCCACAGATGCGA-1’ ... ]]

        AL627309.1 . . .         .         . . . .         . .        
        AL627309.5 . . .         .         . . . .         . .        
        AP006222.2 . . .         .         . . . .         . .        
        LINC01409  . . .         0.7438965 . . . .         . .        
        FAM87B     . . .         .         . . . .         . .        
        LINC01128  . . 0.7991683 .         . . . 0.5091777 . 0.3826447
        LINC00115  . . .         .         . . . .         . .        
        FAM41C     . . .         .         . . . .         . .        
        AL645608.6 . . .         .         . . . .         . .        
        AL645608.2 . . .         .         . . . .         . .  
        ```

```R
gbm <- Seurat::NormalizeData(gbm,
                     normalization.method = "LogNormalize",
                     scale.factor = 10000)
```

!!! note "Updating `gbm`"
    As you might have noticed, this function takes the object `gbm` as input, and it returns it to an object named `gbm`. We can do this because the output of such calculations are added to the object, without loosing information.

### Variable features

We next calculate a subset of features that exhibit high cell-to-cell variation in the
dataset (i.e, they are highly expressed in some cells, and lowly expressed in others).
Focusing on these genes in downstream analysis helps to highlight biological signal
in single-cell datasets.
The procedure in Seurat models the mean-variance relationship inherent in single-cell
data, and is implemented in the `FindVariableFeatures()` function.
By default, 2,000 genes (features) per dataset are returned and these will be used in
downstream analysis, like PCA.

```R
gbm <- Seurat::FindVariableFeatures(gbm,
                            selection.method = "vst",
                            nfeatures = 2000)
```

Let's have a look at the 10 most variable genes:

```R
# Identify the 10 most highly variable genes
top10 <- head(Seurat::VariableFeatures(gbm), 10)
top10
```

We can plot them in a nicely labeled scatterplot:

```R
vf_plot <- Seurat::VariableFeaturePlot(gbm)
Seurat::LabelPoints(plot = vf_plot,
            points = top10, repel = TRUE)
```

### Scaling

Next, we apply scaling, a linear transformation that is a standard pre-processing
step prior to dimensional reduction techniques like PCA. The `ScaleData()` function

1. shifts the expression of each gene, so that the mean expression across cells is 0
2. scales the expression of each gene, so that the variance across cells is 1

This step gives equal weight in downstream analyses, so that highly-expressed genes do not dominate. The results of this are stored in `gbm$RNA@scale.data`

```R
gbm <- Seurat::ScaleData(gbm,
                 features = rownames(gbm))
```

Here's a function to plot mean versus variance in a count matrix (credits go to Wandrille Duchemin):

```R
plot_mean_variance <- function(count_mat){
  # get log10 of mean per gene
  log_mean <- log10(rowMeans(count_mat))
  # get log10 of variance per gene
  log_var <- log10(apply(count_mat, 1, var))
  # get density colors
  cols <- densCols(log_mean, log_var)

  plot(log_mean, log_var, col = cols, pch = 19, cex = 0.2,
        xlab = "log10(mean)", ylab = "log10(variation)")
  # draw the line y=x
  abline(0,1, col = "red")
}
```

**Exercise:** Use the above function to plot the mean variance of:

- raw counts (stored in the slot `"counts"`)
- normalized data (stored in the slot `"data"`)
- scaled data (stored in the slot `"scale.data"`)

!!! hint "Retrieving count data"
    You can retrieve count data from a specific slot with the function `GetAssayData`, e.g. for the raw counts this would be:
    ```R
    Seurat::GetAssayData(gbm, slot = "counts")
    ```

??? done "Answer"
    === "Raw counts"

        Running

        ```R
        plot_mean_variance(Seurat::GetAssayData(gbm, slot =  "counts"))
        ```

        Returns:

        <figure>
          <img src="../../assets/images/plot_mean_var_counts.png" width="400"/>
        </figure>

    === "Normalized data"

        Running

        ```R
        plot_mean_variance(Seurat::GetAssayData(gbm, slot =  "data"))
        ```

        Returns:

        <figure>
          <img src="../../assets/images/plot_mean_var_data.png" width="400"/>
        </figure>

    === "Scaled and normalized"

        Running

        ```R
        plot_mean_variance(Seurat::GetAssayData(gbm, slot =  "scale.data"))
        ```

        As you remember, scaling results in a variance of 1 for each gene, so the `log10` of the variance would be 0. In addition, the mean is scaled to zero, so a `log10` of 0 would return `-Inf` or a very negative number (in case it's very close to 0). It therefore returns:

        <figure>
          <img src="../../assets/images/plot_mean_var_scale.data.png" width="400"/>
        </figure>
