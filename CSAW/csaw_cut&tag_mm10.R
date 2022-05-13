# csaw_workflow.R

# csaw workflow for ATAC-seq differential accessibility analysis, in R
# Jake Reske
# Michigan State University, 2020
# reskejak@msu.edu
# https://github.com/reskejak

# Machine-readable version of Figure 6 workflow for ATAC-seq DA analysis with csaw
# will describe methods for using 1) pre-defined peaks from MACS2 as well as 2) csaw de novo enriched window calling by local enrichment, 
# and normalization methods including 1) TMM on binned counts and 2) loess-based for trended biases

# brief aside: use gc() to help clear memory after intensive commands if crashes/errors occur

# example experimental design: n=2 mouse ATAC-seq biological replicates for two conditions: treat and control 

library(GenomicRanges)
library(csaw)

setwd('/Users/jplab/Desktop/2022-2/data/2-16_CSAW')

########################################
########################################
########################################

# starting from MACS2 filtered narrowPeak

# read replicate narrowPeak files
FOXA1_B14_rep1.peaks <- read.table("/Users/jplab/Desktop/snakepipes_Data/CUT_TAG_snakepipes/FOXA1_H3K27ac/MACS2/FOXA1_B14_1.filtered.short.BAM_peaks.narrowPeak", sep="\t")[,1:3]
FOXA1_B14_rep2.peaks <- read.table("/Users/jplab/Desktop/snakepipes_Data/CUT_TAG_snakepipes/FOXA1_H3K27ac/MACS2/FOXA1_B14_2.filtered.short.BAM_peaks.narrowPeak", sep="\t")[,1:3]
FOXA1_DMSO_rep1.peaks <- read.table("/Users/jplab/Desktop/snakepipes_Data/CUT_TAG_snakepipes/FOXA1_H3K27ac/MACS2/FOXA1_DMSO_1.filtered.short.BAM_peaks.narrowPeak", sep="\t")[,1:3]
FOXA1_DMSO_rep2.peaks <- read.table("/Users/jplab/Desktop/snakepipes_Data/CUT_TAG_snakepipes/FOXA1_H3K27ac/MACS2/FOXA1_DMSO_2.filtered.short.BAM_peaks.narrowPeak", sep="\t")[,1:3]
FOXA1_TGF_rep1.peaks <- read.table("/Users/jplab/Desktop/snakepipes_Data/CUT_TAG_snakepipes/FOXA1_H3K27ac/MACS2/FOXA1_TGF_1.filtered.short.BAM_peaks.narrowPeak", sep="\t")[,1:3]
FOXA1_TGF_rep2.peaks <- read.table("/Users/jplab/Desktop/snakepipes_Data/CUT_TAG_snakepipes/FOXA1_H3K27ac/MACS2/FOXA1_TGF_2.filtered.short.BAM_peaks.narrowPeak", sep="\t")[,1:3]
colnames(FOXA1_B14_rep1.peaks) <- c("chrom", "start", "end")
colnames(FOXA1_B14_rep2.peaks) <- c("chrom", "start", "end")
colnames(FOXA1_DMSO_rep1.peaks) <- c("chrom", "start", "end")
colnames(FOXA1_DMSO_rep2.peaks) <- c("chrom", "start", "end")
colnames(FOXA1_TGF_rep1.peaks) <- c("chrom", "start", "end")
colnames(FOXA1_TGF_rep2.peaks) <- c("chrom", "start", "end")

# read naive overlap broadPeak files
# treat.overlap.peaks <- read.table("treat_overlap_peaks.filt.broadPeak", sep="\t")[,1:3]
# control.overlap.peaks <- read.table("control_overlap_peaks.filt.broadPeak", sep="\t")[,1:3]
# colnames(treat.overlap.peaks) <- c("chrom", "start", "end")
# colnames(control.overlap.peaks) <- c("chrom", "start", "end")

# convert to GRanges objects
FOXA1_B14_rep1.peaks <- GRanges(FOXA1_B14_rep1.peaks)
FOXA1_B14_rep2.peaks <- GRanges(FOXA1_B14_rep2.peaks)
FOXA1_DMSO_rep1.peaks <- GRanges(FOXA1_DMSO_rep1.peaks)
FOXA1_DMSO_rep2.peaks <- GRanges(FOXA1_DMSO_rep2.peaks)
FOXA1_TGF_rep1.peaks <- GRanges(FOXA1_TGF_rep1.peaks)
FOXA1_TGF_rep2.peaks <- GRanges(FOXA1_TGF_rep2.peaks)


# define consensus peakset

# one method: union of all replicate peak sets for both conditions
# treat.peaks <- union(treat1.peaks, treat2.peaks)
# control.peaks <- union(control1.peaks, control2.peaks)
# all.peaks <- union(treat.peaks, control.peaks)

# another method: intersect between biological replicates; union between both experimental conditions
FOXA1_B14.peaks <- intersect(FOXA1_B14_rep1.peaks, FOXA1_B14_rep2.peaks)
FOXA1_DMSO.peaks <- intersect(FOXA1_DMSO_rep1.peaks, FOXA1_DMSO_rep2.peaks)
FOXA1_TGF.peaks <- intersect(FOXA1_TGF_rep1.peaks, FOXA1_TGF_rep2.peaks)
B14.DMSO.peaks <- union(FOXA1_B14.peaks, FOXA1_DMSO.peaks)
all.peaks <- union(B14.DMSO.peaks, FOXA1_TGF.peaks)

# yet another method: union between naive overlapping peak sets
# all.peaks <- union(treat.overlap.peaks, control.overlap.peaks)

##############################
# specify paired-end BAMs
pe.bams <- c("/Users/jplab/Desktop/snakepipes_Data/CUT_TAG_snakepipes/FOXA1_H3K27ac/short_bams/FOXA1_B14_1.short.cleaned.bam",
             "/Users/jplab/Desktop/snakepipes_Data/CUT_TAG_snakepipes/FOXA1_H3K27ac/short_bams/FOXA1_B14_2.short.cleaned.bam",
             "/Users/jplab/Desktop/snakepipes_Data/CUT_TAG_snakepipes/FOXA1_H3K27ac/short_bams/FOXA1_DMSO_1.short.cleaned.bam",
             "/Users/jplab/Desktop/snakepipes_Data/CUT_TAG_snakepipes/FOXA1_H3K27ac/short_bams/FOXA1_DMSO_2.short.cleaned.bam",
             "/Users/jplab/Desktop/snakepipes_Data/CUT_TAG_snakepipes/FOXA1_H3K27ac/short_bams/FOXA1_TGF_1.short.cleaned.bam",
             "/Users/jplab/Desktop/snakepipes_Data/CUT_TAG_snakepipes/FOXA1_H3K27ac/short_bams/FOXA1_TGF_2.short.cleaned.bam"
             )

##############################
# read genome blacklist
# attention: if blacklist bed file have 0 position start lines, when perform regionCounts(),will throw error:
# BiocParallel errors
# element index: 1, 2, 3, 4, 5, 6
# first error: position vector should be 1-based

# i convert blacklist.bed which have 0 start in line to 1 start with preffix 1based, but not other lines. but not sure will conflic with peak bed file;
blacklist <- read.table("/Users/jplab/Desktop/blacklists/1based_hg38-blacklist.v2.bed", sep="\t")
colnames(blacklist) <- c("chrom", "start", "end")
blacklist <- GRanges(blacklist)

# define read parameters
firstChrom1=1
humanChrom22=22
mouseChrom19=19

standard.chr <- paste0("chr", c(firstChrom1:humanChrom22, "X", "Y")) # only use standard chromosomes
param <- readParam(max.frag=1000, pe="both", discard=blacklist, restrict=standard.chr)

##############################
# count reads in windows specified by MACS2                                      
peak.counts <- regionCounts(pe.bams, all.peaks, param=param)

##############################
# MACS2 peaks only: filter low abundance peaks
library("edgeR")
peak.abundances <- aveLogCPM(asDGEList(peak.counts)) 
peak.counts.filt <- peak.counts[peak.abundances > -3, ] # only use peaks logCPM > -3
# few or no peaks should be removed; modify as desired

##############################

# get paired-end fragment size distribution
FOXA1_B14_1.pe.sizes <- getPESizes("/Users/jplab/Desktop/snakepipes_Data/CUT_TAG_snakepipes/FOXA1_H3K27ac/short_bams/FOXA1_B14_1.short.cleaned.bam")
FOXA1_B14_2.pe.sizes <- getPESizes("/Users/jplab/Desktop/snakepipes_Data/CUT_TAG_snakepipes/FOXA1_H3K27ac/short_bams/FOXA1_B14_2.short.cleaned.bam")
FOXA1_B14_1.pe.sizes <- getPESizes("/Users/jplab/Desktop/snakepipes_Data/CUT_TAG_snakepipes/FOXA1_H3K27ac/short_bams/FOXA1_B14_1.short.cleaned.bam")
FOXA1_B14_2.pe.sizes <- getPESizes("/Users/jplab/Desktop/snakepipes_Data/CUT_TAG_snakepipes/FOXA1_H3K27ac/short_bams/FOXA1_B14_2.short.cleaned.bam")
FOXA1_TGF_1.pe.sizes <- getPESizes("/Users/jplab/Desktop/snakepipes_Data/CUT_TAG_snakepipes/FOXA1_H3K27ac/short_bams/FOXA1_TGF_1.short.cleaned.bam")
FOXA1_TGF_2.pe.sizes <- getPESizes("/Users/jplab/Desktop/snakepipes_Data/CUT_TAG_snakepipes/FOXA1_H3K27ac/short_bams/FOXA1_TGF_2.short.cleaned.bam")
gc()
# plot
hist(FOXA1_B14_1.pe.sizes$sizes) # repeat for all replicates and conditions

# for analysis with csaw de novo enriched query windows, select a window size that is greater than the majority of fragments
# max(FOXA1_B14_1.pe.sizes$sizes) gives 150bp because i filter fragments by 150 bp

##############################
# count BAM reads in, e.g. 300 bp windows
counts <- windowCounts(pe.bams, width=200, param=param) # set width as desired from the fragment length distribution analyses

# filter uninteresting features (windows) by local enrichment
# local background estimator: 2kb neighborhood
neighbor <- suppressWarnings(resize(rowRanges(counts), width=2000, fix="center")) # change width parameter as desired
wider <- regionCounts(pe.bams, regions=neighbor, param=param) # count reads in neighborhoods
# filter.stat <- filterWindows(counts, wider, type="local") # the filterWindows() function is deprecated and has been replaced by filterWindowsLocal(). This is an archived step.
filter.stat <- filterWindowsLocal(counts, wider)
counts.local.filt <- counts[filter.stat$filter > log2(3),] # threshold of 3-fold increase in enrichment over 2kb neighborhood abundance; change as desired

###############################
# count BAM background bins (for TMM normalization)
binned <- windowCounts(pe.bams, bin=TRUE, width=10000, param=param)

##########################################
# NORMALIZATION

# method 1: MACS2 peaks only, TMM normalization based on binned counts
peak.counts.tmm <- peak.counts.filt
peak.counts.tmm <- normFactors(binned, se.out=peak.counts.tmm)

# method 2: MACS2 peaks only, csaw loess-normalization
# peak.counts.loess <- peak.counts.filt
# peak.counts.loess <- normOffsets(peak.counts.loess, se.out=TRUE) # type="loess" is now default
# from vignette: "For type="loess", a numeric matrix of the same dimensions as counts, containing the log-based offsets for use in GLM fitting."

# method 3: csaw de novo peaks by local enrichment, TMM normalization based on binned counts
# counts.local.tmm <- counts.local.filt
# counts.local.tmm <- normFactors(binned, se.out=counts.local.tmm)

# method 4: csaw de novo peaks by local enrichment, csaw loess-normalization
# counts.local.loess <- counts.local.filt
# counts.local.loess <- normOffsets(counts.local.loess, se.out=TRUE) # type="loess" is now default
# from vignette: "For type="loess", a numeric matrix of the same dimensions as counts, containing the log-based offsets for use in GLM fitting."

#########################################
# DIFFERENTIAL ACCESSIBILITY ANALYSIS

# set working windows for the desired analysis
working.windows <- peak.counts.tmm # MACS2 peaks only, standard TMM normalization based on binned counts
# working.windows <- peak.counts.loess # MACS2 peaks only, for trended biases
# working.windows <- counts.local.tmm # csaw de novo peaks by local enrichment, standard TMM normalization based on binned counts
# working.windows <- counts.local.loess # csaw de novo peaks by local enrichment, for trended biases
# SEE THE CSAW MANUAL FOR MORE INFO ON NORMALIZATION METHODS
###########

# setup design matrix
# see edgeR manual for more information
y <- asDGEList(working.windows)
colnames(y$counts) <- c("FOXA1_B14_rep1", "FOXA1_B14_rep2", "FOXA1_DMSO_rep1", "FOXA1_DMSO_rep2","FOXA1_TGF_rep1", "FOXA1_TGF_rep2")
rownames(y$samples) <- c("FOXA1_B14_rep1", "FOXA1_B14_rep2", "FOXA1_DMSO_rep1", "FOXA1_DMSO_rep2","FOXA1_TGF_rep1", "FOXA1_TGF_rep2")
y$samples$group <- c("FOXA1_B14", "FOXA1_B14", "FOXA1_DMSO", "FOXA1_DMSO","FOXA1_TGF", "FOXA1_TGF")
design <- model.matrix(~0+group, data=y$samples)
colnames(design) <- c("FOXA1_B14","FOXA1_DMSO", "FOXA1_TGF") # CONFIRM THAT THESE COLUMNS CORRECTLY ALIGN!!
# design
# IMPORTANT: the user should manually confirm that the design matrix is correctly labeled according to sample metadata!

# stabilize dispersion estimates with empirical bayes
y <- estimateDisp(y, design)
fit <- glmQLFit(y, design, robust=TRUE)

# testing for differentially-accessible windows
results <- glmQLFTest(fit, contrast=makeContrasts(FOXA1_B14-FOXA1_DMSO, levels=design))
# head(results$table)
rowData(working.windows) <- cbind(rowData(working.windows), results$table) # combine GRanges rowdata with differential statistics
# working.windows@rowRanges

# merge nearby windows
# up to "tol" distance apart: 500 bp in this case; max merged window width: 5000 bp
merged.peaks <- mergeWindows(rowRanges(working.windows), tol=500L, max.width=5000L)
# summary(width(merged.peaks$region))
# should merge some peaks; change as desired

# use most significant window as statistical representation for p-value and FDR for merged windows
tab.best <- getBestTest(merged.peaks$id, results$table)

# concatenating all relevant statistical data for final merged windows (no redundant columns)
final.merged.peaks <- GRanges(cbind(as.data.frame(merged.peaks$region), results$table[tab.best$rep.test, -4], tab.best[,-c(7:8)]))

# sort by FDR
final.merged.peaks <- final.merged.peaks[order(final.merged.peaks@elementMetadata$FDR), ]
final.merged.peaks

# filter by FDR threshold
FDR.thresh <- 0.05 # set as desired
final.merged.peaks.sig <- final.merged.peaks[final.merged.peaks@elementMetadata$FDR < FDR.thresh, ]
final.merged.peaks.sig # significant differentially-accessible windows

write.table(final.merged.peaks, "FOXA1.B14_vs_FOXA1.DMSO_csaw_DA-windows_all.txt", sep="\t", quote=F, col.names=T, row.names=F)
write.table(final.merged.peaks.sig, "FOXA1.B14_vs_FOXA1.DMSO_csaw_DA-windows_significant.txt", sep="\t", quote=F, col.names=T, row.names=F)

PVALUE.thresh <- 0.05
final.merged.peaks.pvalue.sig <- final.merged.peaks[final.merged.peaks@elementMetadata$PValue < PVALUE.thresh, ]
final.merged.peaks.pvalue.sig # significant differentially-accessible windows

write.table(final.merged.peaks.pvalue.sig, "FOXA1.B14_vs_FOXA1.DMSO_csaw_DA-windows_pvalue_significant.txt", sep="\t", quote=F, col.names=T, row.names=F)

###########################################

# Generate MA plot
library(ggplot2)

final.merged.peaks$sig <- "n.s."
final.merged.peaks$sig[final.merged.peaks$FDR < FDR.thresh] <- "significant"

ggplot(data=data.frame(final.merged.peaks),
       aes(x = logCPM, y = logFC, col = factor(sig, levels=c("n.s.", "significant")))) + 
  geom_point() + scale_color_manual(values = c("black", "red")) + 
  geom_smooth(inherit.aes=F, aes(x = logCPM, y = logFC), method = "loess") + # smoothed loess fit; can add span=0.5 to reduce computation load/time
  geom_hline(yintercept = 0) + labs(col = NULL)
