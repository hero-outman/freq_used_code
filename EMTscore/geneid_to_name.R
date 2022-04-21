library(biomaRt)
library(clusterProfiler)
library('org.Hs.eg.db')
library(AnnotationDbi)

setwd('/Users/jplab/Desktop/DAILY_CODE_DATA/2022-4/code/4-6_EMTScore/counts/')
deseq2_normed_counts <- read.table("./mRNA_B14.titrition.tsv", header=TRUE, sep ="\t")

deseq2_normed_counts$external_gene_name <- mapIds(org.Hs.eg.db,
                                             keys=deseq2_normed_counts$geneid,
                                             column="SYMBOL",
                                             keytype="ENSEMBL",
                                             multiVals="first")


write.table(deseq2_normed_counts, file="mRNA_B14.titrition.tsv",quote=F, sep="\t", row.names=T, col.names=T)
