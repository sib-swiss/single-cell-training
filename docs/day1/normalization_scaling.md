## Learning outcomes

**After having completed this chapter you will be able to:**

- Describe and perform standard procedures for normalization and scaling with the package `Seurat`
- Select the most variable genes from a `Seurat` object for downstream analyses

## Material

- [Seurat vignette](https://satijalab.org/seurat/articles/pbmc3k_tutorial.html)

## Exercises

> :fontawesome-solid-ribbon: This chapter uses the `gbm` dataset

### Normalization

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

!!! note "The use of `Seurat::SCTransform`"
    The functions `NormalizeData`, `VariableFeatures` and `ScaleData` can be replaced by the function `SCTransform`. The latter uses a more sophisticated way to perform the normalization and scaling, and is [argued to perform better](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1874-1). However, it is slower, and a bit less transparent compared to using the three separate functions. Therefore, we chose not to use `SCTransform` for the exercises.
