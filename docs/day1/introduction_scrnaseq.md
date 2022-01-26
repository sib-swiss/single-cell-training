## Learning outcomes

**After having completed this chapter you will be able to:**

- Explain what kind of information single cell RNA-seq can give you to answer a biological question
- Describe essential considerations during the design of a single cell RNA-seq experiment
- Describe the pros and cons of different single cell sequencing methods
- Load single cell data into R
- Explain the basic structure of a `Seurat` object and extract count data and metadata
- Perform a basic quality control by:
    - Evaluating the percentage of UMIs originating from mitochondrial genes
    - Detecting doublets


## Material

[:fontawesome-solid-file-pdf: Download the presentation](../assets/pdf/introduction_scRNAseq.pdf){: .md-button }

- Single cell introductory [video on iBiology](https://www.youtube.com/watch?v=k9VFNLLQP8c)
- Seurat [website](https://satijalab.org/seurat/)
- [Paper](https://doi.org/10.3389/fcell.2018.00108) on experimental considerations
- [Paper](https://doi.org/10.1093/bib/bby007) on experimental design
- [SMART-seq3 protocol](https://www.protocols.io/view/smart-seq3-protocol-bcq4ivyw) at protocols.io
- `cellranger` [system requirements](https://support.10xgenomics.com/single-cell-gene-expression/software/overview/system-requirements) and [installation](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/installation)
- [Review](https://www.nature.com/articles/s41596-020-00409-w?proof=t) by Tallulah Andrews
- [Paper](https://www.biorxiv.org/content/10.1101/749473v3.full) on correlation between mRNA and protein level in single cells

### Running `cellranger count`

Have a look in the directory `course_data/reads` and `reference`. In the `reads` directory you will find reads on one sample: `ETV6-RUNX1_1`. In the analysis part of the course we will work with six samples, but due to time and computational limitations we will run `cellranger count` on one of the samples, and only reads originating from chromsome 21 and 22. 

The only thing you need to run `cellranger count` are the sequence reads and a reference. Here, we have prepared a reference only with chromosome 21 and 22, but in 'real life' you would of course get the full reference genome of your species. The reference has a specific format. You can download precomputed human and mouse references from the [10X website](https://support.10xgenomics.com/single-cell-gene-expression/software/downloads/latest). If your species of interest is not one of those, you will have to generate it yourself. For that, have a look [here](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/tutorial_mr).

Have a look at the documentation of [`cellranger count`](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/count) (scroll down to *Command-line argument reference*), fill out the missing arguments (at `FIXME`) in the script below, and run it:

```sh
cellranger count \
--id=FIXME \
--sample=FIXME \
--transcriptome=FIXME \
--fastqs=FIXME \
--localcores=4 
```

!!! note "This will take a while.."
    Once started, the process will need approximately 15 minutes to finish. Have a coffee and/or have a look at the other exercises. 

!!! hint "Running a bash command with Rstudio"
    You can run a bash script or command using the terminal tab in Rstudio server: 
    <figure>
    <img src="../../assets/images/select_terminal_tab.png" width="300"/>
    </figure>

??? done "Answer"
    ```sh
    cellranger count \
    --id=ETV6-RUNX1_1 \
    --sample=ETV6-RUNX1_1 \
    --transcriptome=cellranger_index \
    --fastqs=course_data/reads \
    --localcores=4
    ```

Have a look out the output directory (i.e. `~/ETV6-RUNX1_1/outs`). The analysis report (`web_summary.html`) is usually a good place to start. 

!!! tip "Open html files in Rstudio server"
    You can use the file browser in the bottom right (tab "Files") to open html files:

    <figure>
    <img src="../../assets/images/open_html.png" width="300"/>
    </figure>



### Loading scRNAseq data

The next step after the generation of the count matrices with `cellranger count`, is the data analysis. The `R` package `Seurat` is currently the most popular software to do this, and therefore we will focus on it in this course. To start working with `Seurat` you can load it into your environment like this:

```R
library(Seurat)
```

!!! tip "Tip: make an R script"
    You could type and copy-paste the commands of these exercises directly in the console. However, that makes it hard to track what you have done. In addition, it can be nice to add comments to your code, so you can read back why you have made certain choices. In order to do that, do not write commands in the console, but write them in a script, and send them to the console with ++ctrl+enter++ (Windows) or ++cmd+enter++ (MacOS).

First, we will load a file specifying the different samples, and create an object specifying the location of the count data:

```R
sample_info <- read.csv("sample_info_course.csv")
datadirs <- file.path("course_data", "count_matrices", sample_info$SampleName,
                      "outs", "filtered_feature_bc_matrix")
names(datadirs) <- gsub("_", "-", sample_info$SampleName)
datadirs
```

The object `datadirs` is a named vector specifying the paths of the count directories for each sample:

```
                                                                  PBMMC-1 
     "course_data/count_matrices/PBMMC_1/outs/filtered_feature_bc_matrix" 
                                                                  PBMMC-2 
     "course_data/count_matrices/PBMMC_2/outs/filtered_feature_bc_matrix" 
                                                                  PBMMC-3 
     "course_data/count_matrices/PBMMC_3/outs/filtered_feature_bc_matrix" 
```

To run through a typical `Seurat` analysis, we will use the files that are in the directory `data/filtered_feature_bc_matrix`. This directory is part of the output generated by [`cellranger`](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/what-is-cell-ranger). To load in this data into R and generate a sparse matrix, run the following command:

```R
sparse_matrix <- Seurat::Read10X(data.dir = datadirs)
```

Basically, this imports a raw count matrix. Have a look at the counts of the first 30 cells of three genes by running:

```R
sparse_matrix[c("PECAM1", "CD8A", "TSPAN1"), 1:30]
```

To generate a `Seurat` object, we will run `CreateSeuratObject`:

```R
seu <- Seurat::CreateSeuratObject(counts = sparse_matrix,
                                  project = "pbmmc")
```

!!! note "Function notation with `::`"
    Here, we define the function together with its associated package. We do that by the syntax `package::function`. Of course, you can also call `library(package)`, and only type the function name. Since we use many different packages in this course, it can be confusing which function comes from which package. Therefore, we chose to always associate the package with the called function.

**Exercise:** check what's in the `seu` object, by typing `seu` in the R console. How many features are in there? And how many cells?

??? done "Answer"
    Typing `seu` should return:

    ```
    An object of class Seurat 
    18673 features across 6946 samples within 1 assay 
    Active assay: RNA (18673 features, 0 variable features)
    ```

    This means that we have 18673 genes (features) in there, and 6946 cells (samples)

### The `Seurat` object

The `seu` object we have created has the class `Seurat`. The object contains multi-level slots and lists. Each `Seurat` object contains exactly the same slots that are specified in the image below (get slot descriptions by typing `?SeuratObject::Seurat`). You can get the information inside a slot with `@`, in the same way as you would use the `$` for lists (e.g. `seu@meta.data` will return the a data frame with information on each cell). Slots can be filled with other R objects, like lists, vectors, data frames or any other class. Here's an overview of all slots that are in a `Seurat` object:

<figure>
  <img src="../../assets/images/seurat_object.png" width="500"/>
</figure>

In addition to the original count table, the `Seurat` object can therefore store a lot of information that is generated during your analysis, like results of a normalization (`@assays$RNA@data`) a PCA or UMAP (`@reductions`) and the clustering (`@graphs`). It also tracks all the commands that have been used to generate the object in its current state (`@commands`). Therefore, while going through the analysis steps, the same object gets more and more of its slots filled. Because most `Seurat` functions return the input object + adjusted slots, we can use this syntax:

```R
seurat_object <- Seurat::function(seurat_object)
```

!!! note "Getting specific information out of the `Seurat` object"
        In order to get specific data you can use the `@` and `$` symbols to browse through the objects. However, `Seurat` comes with a lot of convenience functions, that are easier to use. So, e.g. to get the raw count matrix, you could type `seurat_object@assays$RNA@counts`, however, this is equivalent to `GetAssayData(object = seurat_object, slot = "counts")`. More information on these convenience functions [here](https://satijalab.org/seurat/articles/essential_commands.html).

**Exercise:**

A. Have a look at the `seu` object by running `View(seu)`. What is in there? What is stored in `@active.ident`? 

B. Have a look at the `data.frame` stored at `seu@meta.data` what kind of information is in there?

??? done "Answer"
    === "Answer A"

        There are many slots as described in the above figure. The slot `@active.ident` contains data specifying the samples, e.g. `table(seu@active.ident)` returns:

        ```
        PBMMC-1 PBMMC-2 PBMMC-3 
        1612    3105    2229 
        ```

        Which are the number of cells per sample. 

    === "Answer B"
        Running `head(seu@meta.data)` returns:

        ```
                                   orig.ident nCount_RNA nFeature_RNA
        PBMMC-1_AAACCTGCAGACGCAA-1    PBMMC-1       2401          909
        PBMMC-1_AAACCTGTCATCACCC-1    PBMMC-1       3532          760
        PBMMC-1_AAAGATGCATAAAGGT-1    PBMMC-1       3972         1215
        PBMMC-1_AAAGCAAAGCAGCGTA-1    PBMMC-1       3569          894
        PBMMC-1_AAAGCAACAATAACGA-1    PBMMC-1       2982          730
        PBMMC-1_AAAGCAACATCAGTCA-1    PBMMC-1      22284         3108
        ```

        Giving you the names of three columns and a row for each cell:

        * `orig_ident`: the original identity (origin) of a cell.
        * `nCount_RNA`: the number of reads assigned to a cell.
        * `nFeature_RNA`: the number of expressed features (genes) per cell.

Luckily, usually you do not have to dive into this structure to retrieve information. For example, information in the slot `@meta.data` can be retrieved and set by using `$` or `[[]]`.

!!! note
    There is a subtle difference here between `$` and `[[]]`. While `$` returns a vector of the column in `@meta.data`, `[[]]` returns a `data.frame`.


 **Exercise:** Generate a histogram of the column `nCount_RNA` at `seu@meta.data`, with the base function `hist`.

??? done "Answer"
    ```R
    hist(seu$nCount_RNA)
    ```

    or

    ```R
    hist(seu@meta.data$nCount_RNA)
    ```

There are also built-in functions to plot data from `Seurat` object, for example `FeatureScatter`. This function enables you easily draw a scatterplot from a `Seurat` object:

```R
Seurat::FeatureScatter(seu, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
```

### Basic Quality Control

A high number of UMIs originating from mitochondrial genes can point to dying cells. In order to calculate the percentage of UMIs coming from mitochondrial genes for each cell, we use the function `PercentageFeatureSet`. In our count matrix, the names of mitochondrial genes all start with `MT-`, so we can use that pattern to search these genes:

```R
seu[["percent.mt"]] <- Seurat::PercentageFeatureSet(seu, pattern = "^MT-")
```

!!! note "Finding mitochondrial genes for mouse"
    For mouse, mitochondrial genes start with `mt-`, so the function call with look like this:

    ```
    seu[["percent.mt"]] <- Seurat::PercentageFeatureSet(seu, pattern = "^mt-")
    ```

Now check out whether the column was added to the `meta.data` slot:

```R
head(seu@meta.data)
```

Let's have a look at the distribution of the three columns stored in `@meta.data`:

```R
Seurat::VlnPlot(seu, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
```

There is often a relationship between number of UMIs and percentage of mitochondrial genes. We can have a look at their relationship:

```R
Seurat::FeatureScatter(seu, feature1 = "nCount_RNA", feature2 = "percent.mt")
```

Based on what we know now, it would be sensible to keep cells in which less than 15% of the UMIs come from mitochondrial genes. We can filter like this:

```R  
seu <- subset(seu,
              subset = percent.mt < 15)
```

**Exercise:** How many cells did we filter out?

??? done "Answer"
    Just running `seu` returns:
    ```
    An object of class Seurat 
    18673 features across 6939 samples within 1 assay 
    Active assay: RNA (18673 features, 0 variable features)
    ```

    Meaning that 6946 - 6939 = 7 cells were filtered out.
