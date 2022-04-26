library(DESeq2)
library(DEGreport)
library(ggplot2)
library(pheatmap)
library(amap)
library("RColorBrewer")
library(clusterProfiler)
library('org.Hs.eg.db')
library(AnnotationDbi)
library(ggpubr)
library(ggthemes)

setwd('/Users/jplab/Desktop/DAILY_CODE_DATA/2022-3/data/3-30_noTGF_shFOXA1_deseq2_norm/')

######################## read data and set coldata ########################
# sample name and order
# sampleNames <- c("nc_B14.rep1", "nc_B14.rep2", "nc_DMSO.rep1", "nc_DMSO.rep2", "nc_TGFb.rep1", "nc_TGFb.rep2", "shFOXA1_B14.rep1", "shFOXA1_B14.rep2", "shFOXA1_DMSO.rep1", "shFOXA1_DMSO.rep2", "shFOXA1_TGFb.rep1", "shFOXA1_TGFb.rep2")

# load gc normed counts, low expr was removed already
read_counts <- read.table("/Users/jplab/Desktop/DAILY_CODE_DATA/2022-3/data/3-30_noTGF_shFOXA1_deseq2_norm/counts_noVersion.tsv", header=TRUE, sep ="\t")

sampleNames <- colnames(read_counts)[2:13]
geneid <- read_counts[1]

rownames(read_counts) <- read_counts[,1] # first column(gene_id) for row names
read_counts <- read_counts[,c(2,3,4,5,6,7,8,9,10,11,12,13)]


# 去掉在多于两个样本中counreat<1的值，如果样本数多，这个数值可以适当增加;排除极低表达基因的干扰
# read_counts <- read_counts[rowSums(read_counts)>2,]

# set sample name
names(read_counts) <- sampleNames
# check if have NA elements
na_ele <- read_counts[rowSums(is.na(read_counts)) > 0, ] 
# take off NA 
read_counts <- na.omit(read_counts)


# knockdown or not
knockdown <- rep(c('nc','sh'),times=c(6,6) )

# drug added condition
# treatment <- rep(c('B14','B14','DMSO','DMSO','TGFb','TGFb'),times=2)
treatment <- c('shFOXA1_B14', 'shFOXA1_B14', 'shFOXA1_DMSO', 'shFOXA1_DMSO', 'shFOXA1_TGFb', 'shFOXA1_TGFb', 'shNC_B14', 'shNC_B14', 'shNC_DMSO', 'shNC_DMSO', 'shNC_TGFb', 'shNC_TGFb')

# coldata
colData <- data.frame(name=sampleNames, 
                      # knockdown=factor(knockdown),
                      treatment=factor(treatment))

# re-level conditions: set control as ref in each condition
colData$treatment <- relevel(colData$treatment,ref="TGFb")

# colData$treatment <- relevel(colData$treatment,ref="DMSO")

colData$knockdown <- relevel(colData$knockdown,ref="nc")

# coldata rownames
rownames(colData) <- sampleNames


######################## DESeq2 for DEG ########################
# deseq2 will use coldata's conditon as dds's condtion; when re-level above, dds condition will also re-leveled
dds <- DESeqDataSetFromMatrix(
  countData = read_counts, 
  colData=colData, 
  design = ~ treatment
  # design= ~ knockdown + treatment + knockdown:treatment
  )
dds <- DESeq(dds)

# norm by size factor
normalized_counts <- counts(dds, normalized=TRUE)
head(normalized_counts)




# 根据基因在不同的样本中表达变化的差异程度mad值对数据排序，差异越大的基因排位越前
# 作用是把表达变化大的放到前面
# apply ：apply函数，对矩阵或list中每一行或列进行处理
# 1 代表行
# mad ，估算离散程度的方法
normalized_counts_mad <- apply(normalized_counts, 1, mad)
normalized_counts <- normalized_counts[order(normalized_counts_mad, decreasing=T), ]

# 标准化后的数据输出; 表头有偏移，手动拖一下
write.table(normalized_counts, file="counts_deseq2Normed.tsv",quote=F, sep="\t", row.names=T, col.names=T)

knockdown_sh_B14_vs_TGF <- list(c('treatment_B14_vs_TGFb','knockdownsh.treatmentB14'))
res_knockdown_sh_B14_vs_TGF<- results(dds, contrast=knockdown_sh_B14_vs_TGF)
res_knockdown_sh_B14_vs_TGF$symbol <- mapIds(org.Hs.eg.db,
                                            keys=rownames(res_knockdown_sh_B14_vs_TGF),
                                            column="SYMBOL",
                                            keytype="ENSEMBL",
                                            multiVals="first")
res_knockdown_sh_B14_vs_TGF$entrez <- mapIds(org.Hs.eg.db,
                                            keys=rownames(res_knockdown_sh_B14_vs_TGF),
                                            column="ENTREZID",
                                            keytype="ENSEMBL",
                                            multiVals="first")
head(res_knockdown_sh_B14_vs_TGF,5)
res_knockdown_sh_B14_vs_TGF <- merge(as.data.frame(res_knockdown_sh_B14_vs_TGF), as.data.frame(normalized_counts),by="row.names",sort=FALSE)

res_knockdown_sh_B14_vs_TGF_compare <- "sh_B14_vs_TGF"
res_knockdown_sh_B14_vs_TGF$compare <- res_knockdown_sh_B14_vs_TGF_compare
res_knockdown_sh_B14_vs_TGF$change <- as.factor(
  ifelse(
    res_knockdown_sh_B14_vs_TGF$padj<=0.05 & abs(res_knockdown_sh_B14_vs_TGF$log2FoldChange)>=1,
    ifelse(res_knockdown_sh_B14_vs_TGF$log2FoldChange>=1, "Up", "Down"),
    "NoDiff"
  )
)

write.table(res_knockdown_sh_B14_vs_TGF, file="sh_B14_vs_TGF.tsv",quote=F, sep="\t", row.names=T, col.names=NA)


# use diff contrast to got diff results
# use resultsNames(dds) to show contrast
treatment_B14_sh_vs_nc <- list(c('knockdown_sh_vs_nc','knockdownsh.treatmentB14'))

# use contrast to get the result table
res_treatment_B14_sh_vs_nc <- results(dds, contrast=treatment_B14_sh_vs_nc)

# use AnnotationDbi convert ENSEMBL id to gene symbol and ENTREZID
res_treatment_B14_sh_vs_nc$symbol <- mapIds(org.Hs.eg.db,
                     keys=rownames(res_treatment_B14_sh_vs_nc),
                     column="SYMBOL",
                     keytype="ENSEMBL",
                     multiVals="first")
res_treatment_B14_sh_vs_nc$entrez <- mapIds(org.Hs.eg.db,
                     keys=rownames(res_treatment_B14_sh_vs_nc),
                     column="ENTREZID",
                     keytype="ENSEMBL",
                     multiVals="first")
head(res_treatment_B14_sh_vs_nc,5)




treatment_DMSO_sh_vs_nc <- list(c('knockdown_sh_vs_nc','knockdownsh.treatmentDMSO'))
res_treatment_DMSO_sh_vs_nc <- results(dds, contrast=treatment_DMSO_sh_vs_nc)

# use AnnotationDbi convert ENSEMBL id to gene symbol and ENTREZID
res_treatment_DMSO_sh_vs_nc$symbol <- mapIds(org.Hs.eg.db,
                                            keys=rownames(res_treatment_DMSO_sh_vs_nc),
                                            column="SYMBOL",
                                            keytype="ENSEMBL",
                                            multiVals="first")
res_treatment_DMSO_sh_vs_nc$entrez <- mapIds(org.Hs.eg.db,
                                            keys=rownames(res_treatment_DMSO_sh_vs_nc),
                                            column="ENTREZID",
                                            keytype="ENSEMBL",
                                            multiVals="first")
head(res_treatment_DMSO_sh_vs_nc,5)
res_treatment_DMSO_sh_vs_nc <- merge(as.data.frame(res_treatment_DMSO_sh_vs_nc), as.data.frame(normalized_counts),by="row.names",sort=FALSE)

res_treatment_DMSO_sh_vs_nc_compare <- "DMSO_sh_vs_nc"
res_treatment_DMSO_sh_vs_nc$compare <- res_treatment_DMSO_sh_vs_nc_compare
res_treatment_DMSO_sh_vs_nc$change <- as.factor(
  ifelse(
    res_treatment_DMSO_sh_vs_nc$padj<=0.05 & abs(res_treatment_DMSO_sh_vs_nc$log2FoldChange)>=1,
    ifelse(res_treatment_DMSO_sh_vs_nc$log2FoldChange>=1, "Up", "Down"),
    "NoDiff"
  )
)

write.table(res_treatment_DMSO_sh_vs_nc, file="DMSO_sh_vs_nc.tsv",quote=F, sep="\t", row.names=T, col.names=NA)


# check size factor of dds
# sizeFactors(dds) 

# counts method normalized : raw_count / sizeFactor
# head(counts(dds, normalized=TRUE),5)

# merge deg result(no gene counts) with deseq2 normalized counts
res_treatment_B14_sh_vs_nc <- merge(as.data.frame(res_treatment_B14_sh_vs_nc), as.data.frame(normalized_counts),by="row.names",sort=FALSE)

# add up and down column to result_deseq_norm_merged by logFC and padj
# padj 0.05
# logFC 1
res_treatment_B14_sh_vs_nc_compare <- "B14_sh_vs_nc"
res_treatment_B14_sh_vs_nc$compare <- res_treatment_B14_sh_vs_nc_compare
res_treatment_B14_sh_vs_nc$change <- as.factor(
  ifelse(
    res_treatment_B14_sh_vs_nc$padj<=0.05 & abs(res_treatment_B14_sh_vs_nc$log2FoldChange)>=1,
    ifelse(res_treatment_B14_sh_vs_nc$log2FoldChange>=1, "Up", "Down"),
    "NoDiff"
  )
)




res_treatment_TGF_sh_vs_nc <- results(dds, contrast=c('knockdown','sh','nc'))
res_treatment_TGF_sh_vs_nc$symbol <- mapIds(org.Hs.eg.db,
                                            keys=rownames(res_treatment_TGF_sh_vs_nc),
                                            column="SYMBOL",
                                            keytype="ENSEMBL",
                                            multiVals="first")
res_treatment_TGF_sh_vs_nc$entrez <- mapIds(org.Hs.eg.db,
                                            keys=rownames(res_treatment_TGF_sh_vs_nc),
                                            column="ENTREZID",
                                            keytype="ENSEMBL",
                                            multiVals="first")
head(res_treatment_TGF_sh_vs_nc)

res_treatment_TGF_sh_vs_nc <- merge(as.data.frame(res_treatment_TGF_sh_vs_nc), as.data.frame(normalized_counts),by="row.names",sort=FALSE)


# add up and down column to result_deseq_norm_merged by logFC and padj
# padj 0.05
# logFC 1
res_treatment_TGF_sh_vs_nc_compare <- "TGFb_sh_vs_nc"
res_treatment_TGF_sh_vs_nc$compare <- res_treatment_TGF_sh_vs_nc_compare
res_treatment_TGF_sh_vs_nc$change <- as.factor(
  ifelse(
    res_treatment_TGF_sh_vs_nc$padj<=0.05 & abs(res_treatment_TGF_sh_vs_nc$log2FoldChange)>=1,
    ifelse(res_treatment_TGF_sh_vs_nc$log2FoldChange>=1, "Up", "Down"),
    "NoDiff"
  )
)


head(res_treatment_TGF_sh_vs_nc)

# save deg with up down nodiff
write.table(res_treatment_TGF_sh_vs_nc, file="TGFb_sh_vs_nc.tsv",quote=F, sep="\t", row.names=T, col.names=NA)

# get sigDiff_deg_subset from deg result
# sigDiff was filtered by padj < 0.05 and foldchange >= 1
diff_genes <- subset(res_treatment_B14_sh_vs_nc, padj < 0.05 & abs(log2FoldChange) >= 1.0)

dim(diff_genes)

summary(diff_genes)
######################## plotting ########################
# before plot, need log transfer
# https://www.rdocumentation.org/packages/DESeq2/versions/1.12.3/topics/rlog
# blind： 
#   logical, whether to blind the transformation to the experimental design. 
#   blind=TRUE should be used for comparing samples in an manner unbiased by prior information on samples, 
#   for example to perform sample QA (quality assurance). blind=FALSE should be used for transforming data for downstream analysis,
#   where the full use of the design information should be made. blind=FALSE will skip re-estimation of the dispersion trend, 
#   if this has already been calculated. If many of genes have large differences in counts due to the experimental design, 
#   it is important to set blind=FALSE for downstream analysis.
rld <- rlog(dds, blind=FALSE)
#Get 25 top varying genes
topVarGenes <- head( order( rowVars( assay(rld) ), decreasing=TRUE ), 80)
#make a subset of the log transformed counts for just the top 25 varying genes
top25Counts<-assay(rld)[topVarGenes,]


# use ggpubr plot volcano plot
# these steps work on adding gene name label to dataframe then show on the plotting
res_treatment_B14_sh_vs_nc$logp <- -log10(res_treatment_B14_sh_vs_nc$padj+.Machine$double.xmin)
res_treatment_B14_sh_vs_nc$label <- ""
res_treatment_B14_sh_vs_nc <- res_treatment_B14_sh_vs_nc[order(res_treatment_B14_sh_vs_nc$padj, -abs(res_treatment_B14_sh_vs_nc$log2FoldChange)), ] # order by padj asc and abs(logFC) desc
up.gene <- head(res_treatment_B14_sh_vs_nc$symbol[which(res_treatment_B14_sh_vs_nc$change == "Up")],10)
down.gene <- head(res_treatment_B14_sh_vs_nc$symbol[which(res_treatment_B14_sh_vs_nc$change == "Down")],10)
up.gene
down.gene
deg.top10.genes <- c(as.character(up.gene),as.character(down.gene))
deg.top10.genes
res_treatment_B14_sh_vs_nc$label[match(deg.top10.genes,res_treatment_B14_sh_vs_nc$symbol)] <- deg.top10.genes

# plot volcano
# add label to dataframe then show on the plotting
ggscatter(
  res_treatment_B14_sh_vs_nc, 
  x="log2FoldChange",
  y="logp",
  color = "change",
  palette = c("#2f5688","#BBBBBB","#CC0000"),
  label = res_treatment_B14_sh_vs_nc$label,
  font.label = 8,
  xlab = "log2FoldChange",
  ylab = "-log10(padj)",
  title = "DEG_B14vsTGF",
  repel = T,
  label.rectangle = T,
  max.overlaps = 20,
  size = 0.8
) + 
  ylim(0, 138) + # 修改y轴比例
  xlim(-3, 5) + # 修改x轴比例
  theme_base()+ 
  # geom_point(alpha=0.3, size=1) +
  geom_hline(yintercept=-log10(0.01), lty=2, col="gray", lwd=0.5)+
  geom_vline(xintercept=c(-1, 1), lty=2, col="gray", lwd=0.5)


# good looking simple sample distance cluster
sampleDists <- dist(t(assay(rld)))

library(factoextra)
hdis <- hcut(sampleDists, k = 2, stand = TRUE)
# Visualize
fviz_dend(hdis, 
          # 加边框
          rect = TRUE,
          # 边框颜色
          rect_border="cluster",
          # 边框线条类型
          rect_lty=2,
          # 边框线条粗细
          lwd=1.2,
          # 边框填充
          rect_fill = T,
          # 字体大小
          cex = 1,
          # 字体颜色
          color_labels_by_k=T,
          # 平行放置
          horiz=T)

## heat map plotting
# heat map of sample distance
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(rld$condition, rld$type, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)

# DEG heatmap
select <- order(rowMeans(normalized_counts), decreasing=T)[1:1000]
nt <- normTransform(dds)
log2.norm.counts <- assay(nt)[select,]
df <- as.data.frame(colData(dds)[, c("name", "condition")])
pheatmap(log2.norm.counts, cluster_rows=T, show_rownames=F, cluster_cols=T, annotation_col=df, fontsize=10,display_numbers = F)

# top 50 DEGs heatmap
topVarGenes <- head(order(rowVars(assay(rld)), decreasing = TRUE), 100)
mat <- assay(rld)[ topVarGenes, ]

symbol <- mapIds(org.Hs.eg.db,
                     keys=rownames(mat),
                     column="SYMBOL",
                     keytype="ENSEMBL",
                     multiVals="first")

rownames(mat) <- symbol

mat <- mat - rowMeans(mat)
anno <- as.data.frame(colData(rld))
pheatmap(
  mat,
  show_colnames = TRUE,
  show_rownames = TRUE,
  annotation_col = anno,
  display_numbers = F,
  main='top_100_deg_genes_heatmap_B14vsDMSOvsTGF')

######################## enrichment ########################

# annotation DEG
library(DOSE)
library(topGO)
keytypes(org.Hs.eg.db) # check database ID type

# GO
df <- diff_genes[order(diff_genes$pvalue),]
ego <- enrichGO(
  # gene = row.names(diff_genes),
  gene = df$Row.names,
  OrgDb = org.Hs.eg.db,
  keyType = "ENSEMBL",
  ont = "MF",
  readable = T
)

# dot plot
dotplot(ego,font.size=10)
# bar plot
barplot(ego,showCategory = 15)
# GO graph
plotGOgraph(ego)

goplot(ego)

emapplot(ego, showCategory = 30)

#GO term与差异基因关系网络图
cnetplot(ego, showCategory = 5)

######################## KEGG ########################
geneId = res_B14vsDMSO$entrez
kegg <- enrichKEGG(
  gene = geneId,
  organism = "hsa",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  qvalueCutoff = 0.05
)

# dot plot
dotplot(kegg)
# bar plot
barplot(kegg,showCategory = 15)

# browse 
browseKEGG(kegg,"hsa285203")

#pathway关系网络图（通过差异基因关联）
emapplot(kegg,  showCategory = 30)

#pathway与差异基因关系网络图
cnetplot(kegg, showCategory = 5)

#pathway映射
browseKEGG(kegg, "hsa04934") #在pathway通路图上标记富集到的基因，会链接到KEGG官网

######################## EGF/EGFR genes DEG  ########################
library(dplyr)

# 获取标准化的dds count矩阵
normalized_counts <- counts(dds,normalized=T)

dim(normalized_counts)
dim(diff_gene_deseq2)

write.csv(normalized_counts, "normalized_counts.csv") # 生成标准化的矩阵

# 1.用diff行匹配，获取padj和logfc筛选过的normalized_counts
# 2.筛选过的normalized_counts选取含有对应egf通路基因的row，获得小的矩阵

# 筛选好的标准化的矩阵
deg_target_gene <- read.csv("./normalized_deg_df.csv", quote="\t")

head(deg_target_gene)

symbol <- mapIds(org.Hs.eg.db,
                 keys=deg_target_gene$Unnamed..0,
                 column="SYMBOL",
                 keytype="ENSEMBL",
                 multiVals="first")

rownames(deg_target_gene) <- symbol

# 第一列去掉
deg_target_gene <- select(deg_target_gene,-c(1))

deg_target_gene_for_heatmap <- deg_target_gene

library(reshape)

melted_deg_target_gene <- melt(as.matrix(deg_target_gene))

colnames(melted_deg_target_gene) <- c("gene","samplename","normalized_counts")

colData$samplename <- colnames(deg_target_gene)
melted_deg_target_gene <- merge(melted_deg_target_gene,colData)

# plotting EGF/EGFR deg genes
ggplot(melted_deg_target_gene)+
  geom_point(aes(x=gene,y=normalized_counts,color=condition))+
  # scale_y_log10()+
  xlab("Genes")+
  ylab("Normalized Counts")+
  ggtitle("EGF/EGFR Signaling Pathway Genes DEG")
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45,hjust = 1))+
  theme(plot.title = element_text(hjust = 0.5))
  
  
anno <- as.data.frame(colData(rld))
heat.colors <- brewer.pal(6,"YlOrRd")
mat_egf <- as.matrix(deg_target_gene_for_heatmap)

egf_heatmap <- pheatmap(
  mat_egf,
  color=heat.colors, 
  cluster_rows = T,show_colnames = T, show_rownames = T,
  # annotation = anno,
  border_color = NA,
  fontsize = 10,scale = "row",fontsize_row = 10,height = 20)
egf_heatmap



######################## DO ########################
doAna <- enrichDO(
  gene = geneId,
  ont = "DO",
  pvalueCutoff = 0.05,
  qvalueCutoff = 1,
  readable = TRUE,
  minGSSize = 1,
  maxGSSize = 500
)

dotplot(doAna)
barplot(doAna,showCategory = 15)

######################## DEG report ########################
degCheckFactors(normalized_counts)

res_B14vsTGF <- results(dds,contrast = B14vsTGF_contrast)
res_B14vsDMSO <- results(dds,contrast = B14vsDMSO_contrast)
res_TGFvsDMSO <- results(dds,contrast = TGFvsDMSO_contrast)

degQC(normalized_counts, group_info, pvalue = res_B14vsTGF[["pvalue"]])
degQC(normalized_counts, group_info, pvalue = res_B14vsDMSO[["pvalue"]])
degQC(normalized_counts, group_info, pvalue = res_TGFvsDMSO[["pvalue"]])

resCov <- degCovariates(log2(normalized_counts+0.5),
                        colData(dds))

cor <- degCorCov(colData(dds))


createReport(colData(dds)[["condition"]], counts(dds, normalized = TRUE),
             row.names(res_B14vsTGF)[1:20], res_B14vsTGF[["pvalue"]], path = "/Users/jplab/Desktop/RNA_downstream/raw_count_flow/DEG")

resreport <- degResults(dds = dds, name = "test", org = NULL,
                        do_go = TRUE, group = "condition", xs = "condition",
                        path_results = NULL)