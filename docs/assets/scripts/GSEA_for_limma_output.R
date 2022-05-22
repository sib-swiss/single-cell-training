# GSEA example using gseGO(), with t-statistic of limma output:
# create ranked list of t-statistics:
gene.list<-limma_de$t
names(gene.list)<-rownames(limma_de)
gene.list<-sort(gene.list, decreasing = T)
head(gene.list)
#    SOCS2     CD74   RPS4Y2 HLA-DPA1    HLA-C     DNTT 
# 44.04190 43.54079 42.65128 38.54829 37.13147 36.50141 

set.seed(1234)
ego<-gseGO(geneList = gene.list,
           ont="BP",
           OrgDb = "org.Hs.eg.db",
           keyType = "SYMBOL",
           minGSSize = 60,
           eps=1e-60,
           seed=T)
ego <- clusterProfiler::simplify(ego)
head(ego@result[,c(2:7)])
#                                                          Description setSize enrichmentScore       NES       pvalue     p.adjust
# GO:0007059                                    chromosome segregation     297      -0.4562146 -3.164213 1.880259e-35 2.918162e-32
# GO:0000070                      mitotic sister chromatid segregation     160      -0.5447013 -3.552161 1.310400e-32 1.016870e-29
# GO:0000280                                          nuclear division     368      -0.4013811 -2.834243 6.824621e-31 2.647953e-28
# GO:0098813                            nuclear chromosome segregation     237      -0.4602798 -3.132855 1.113668e-29 2.469162e-27
# GO:0044770                               cell cycle phase transition     436      -0.3175543 -2.271288 5.563948e-19 1.079406e-16
# GO:1902850 microtubule cytoskeleton organization involved in mitosis     137      -0.4678129 -2.986132 3.639766e-18 6.276573e-16

# barcode plot of the top GO term:
gseaplot(ego, geneSetID = "GO:0007059", title="GO:0007059, chromosome segregation")
