##New script for generating pseudobulk
##This script follows the vignette on this page 
#http://bioconductor.org/books/3.14/OSCA.multisample/multi-sample-comparisons.html
#here 6-11 lines are the same as in the course
#taking the proB data 
Seurat::DefaultAssay(proB) <- "RNA"
Seurat::Idents(proB) <- proB$orig.ident

## add the patient id also for paired DGE
proB$patient.id<-gsub("ETV6-RUNX1", "ETV6_RUNX1", proB$orig.ident)
proB$patient.id<-sapply(strsplit(proB$patient.id, "-"), '[', 2)

##here it is new
##first a mandatory column of sample needs to be added to the meta data that is the grouping factor, should be the samples
proB$sample <- factor(proB$orig.ident)

##first an sce object is needed
sce_proB <- as.SingleCellExperiment(proB)

#The needed package
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("scuttle")

library(scuttle)

##aggregateAcrossCells here it is only aggregated by sample, one could imagine
##to aggregate by sample and by celltype for instance
summed <- aggregateAcrossCells(sce_proB, 
                               id=colData(sce_proB)[,c("sample")])

##have a look at the counts
counts(summed)[1:3,]

#have a look at the colData of our new object summed, can you see type and 
#patient.id are there
head(colData(summed))

#As in the standard limma analysis generate a DGE object

y <- DGEList(counts(summed), samples=colData(summed)$sample)
 
##filter lowly expressed (recommanded for limma)
keep <- filterByExpr(y, group=summed$type)
y <- y[keep,]
  
##see how many genes were kept 
summary(keep)
  
## Create the design matrix and include the technology as a covariate:
design <- model.matrix(~0 + summed$type + summed$patient.id)
  
# Have a look
design

# change column/rownames names to more simple group names: 
 colnames(design) <- make.names(c("ETV6-RUNX1", "PBMMC","patient2","patient3"))
 rownames(design)<-colData(summed)$sample

# Create contrasts, i.e. specify which groups we want to compare, here we want
# to find genes differentially expressed between cluster 1 and cluster 2.
  contrast.mat <- limma::makeContrasts(ETV6.RUNX1 - PBMMC,
                                       levels = design)
  
  
#FROM HERE IT IS THE SAME AS BEFORE
dge <- edgeR::calcNormFactors(y)  

#Do limma
  vm <- limma::voom(dge, design = design, plot = TRUE)
  fit <- limma::lmFit(vm, design = design)
  fit.contrasts <- limma::contrasts.fit(fit, contrast.mat)
  fit.contrasts <- limma::eBayes(fit.contrasts)
  
  # Show the top differentially expressed genes:
  limma::topTable(fit.contrasts, number = 10, sort.by = "P")
  limma_de <- limma::topTable(fit.contrasts, number = Inf, sort.by = "P")
  length(which(limma_de$adj.P.Val<0.05))
 