bedtools intersect -wa -u -a FOXA1_B14vsTGF_P0.05_up.bed -b FOXA1_B14vsDMSO_nothres_up.bed > FOXA1_BDT_up.bed


bedtools intersect -wa -u -a FOXA1_BDT_up.bed -b FOXA1_B14vsDMSO_nothres_up.bed > FOXA1_BDT_up.bed


bedtools intersect -wa -u -a FOXA1_BDT_up.bed -b B14.filtered.short.BAM_peaks.bed > FOXA1_BDT_up_B14ATAC.bed


bedtools intersect -wa -u -a FOXA1_BDT_up_B14ATAC.bed -b H3K27ac_B14.filtered.short.BAM.bed > FOXA1_BDT_up_B14ATAC_H3K27acB14.bed


# exclude some regions
bedtools intersect \
-a /Users/jplab/Desktop/snakepipes_Data/ATAC_seq_shFOXA1_snakepipes_results/CSAW_MACS2/FOXA1.B14_vs_ncFOXA1_B14/bed/FOXA1.B14_vs_ncFOXA1_B14_FDRsig_up.bed \
-b /Users/jplab/Desktop/snakepipes_Data/ATAC_seq_shFOXA1_snakepipes_results/MACS_allLen_noDUP_noChrM/common_open_peak/F_B_NC_B_primary.common.open.bed \
-v \
> B14_SHvsNC_SHuniqUp.bed