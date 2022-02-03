## Material

- [MSigDB](http://www.gsea-msigdb.org/gsea/msigdb/index.jsp)
- `clusterProfiler` [vignette](https://bioconductor.org/packages/release/bioc/vignettes/clusterProfiler/inst/doc/clusterProfiler.html)
- [Revigo](http://revigo.irb.hr/)
- [Signaling Pathway Impact Analysis (SPIA)](https://bioconductor.org/packages/release/bioc/html/SPIA.html)
- Original [paper](https://www.pnas.org/content/102/43/15545) on GSEA
- [STRING](https://string-db.org/) for protein-protein interactions
- [GO figure!](https://gitlab.com/evogenlab/GO-Figure) for plotting GO terms and the associated [paper](https://www.frontiersin.org/articles/10.3389/fbinf.2021.638255/full)

## Exercises

Load the following packages:

```R
library(clusterProfiler)
library(enrichplot)
```

If the `FindMarkers` or `FindAllMarkers` functions were used,
we have a table containing only the significant genes,
but we don't have any information for the non-significant
genes. Therefore, we can use the over-representation analysis
which is a threshold-based method.
Using our list of significant genes, we can test
if any gene set is over-represented in our data or not using a test
similar to a Fisher test to compare differences in proportions.

The `clusterProfiler` package provides functions for over-representation
analysis of Gene Ontology gene sets (among other functions) or KEGG gene sets.

Genes can be labeled using different types of labels, eg
symbol, ensembl ID, Entrez ID. To list the allowed
label types use:

```R
BiocManager::install("org.Hs.eg.db", update = FALSE)
library(org.Hs.eg.db)
AnnotationDbi::keytypes(org.Hs.eg.db)
```

Let's select a set of genes that are upregulated in the astrocytes compared to the macrophages:

```R
tum_down <- subset(tum_vs_norm,
                   tum_vs_norm$avg_log2FC < -1 &
                     tum_vs_norm$p_val_adj < 0.05)
tum_down_genes <- rownames(tum_down)
```

We can do a gene ontology term enrichment analysis based on this set of genes:

```R
tum_vs_norm_go <- clusterProfiler::enrichGO(tum_down_genes,
                                            "org.Hs.eg.db",
                                            keyType = "SYMBOL",
                                            ont = "BP",
                                            minGSSize = 50)
```

The results are stored in the `@result` slot:

```R
enr_go <- clusterProfiler::simplify(tum_vs_norm_go)
View(enr_go@result)
```

We can quite easily generate an enrichment map with the `enrichplot` package:

```R
enrichplot::emapplot(enrichplot::pairwise_termsim(enr_go),
                     showCategory = 30, cex_label_category = 0.5)
```

In stead of testing for gene ontology terms, we can also test for other gene set collections. For example the hallmark collection from [MSigDB](http://www.gsea-msigdb.org/gsea/msigdb/index.jsp):

```R
gmt <- msigdbr::msigdbr(species = "human", category = "H")
```

We can use the function `enricher` to test for enrichment of any set of genes. But we would have to test it against a "universe", i.e. the background genes:

```R
tum_vs_norm_enrich <- clusterProfiler::enricher(gene = tum_down_genes,
                                                universe = rownames(all_prob),
                                                pAdjustMethod = "BH",
                                                pvalueCutoff  = 0.05,
                                                qvalueCutoff  = 0.05,
                                                TERM2GENE = gmt[,c("gs_name", "gene_symbol")])
```

The most signifcantly enriched group of genes is `HALLMARK_G2M_CHECKPOINT`:

```R
View(tum_vs_norm_enrich@result)
```

### Save the dataset and clear environment

Now, save the dataset so you can use it later today:

```R
saveRDS(seu_int, "seu_int_day3.rds")
```

Clear your environment:

```R
rm(list = ls())
gc()
.rs.restartR()
```
