bedtools intersect -wa -u -a FOXA1_B14vsTGF_P0.05_up.bed -b FOXA1_B14vsDMSO_nothres_up.bed > FOXA1_BDT_up.bed



bedtools intersect -wa -u -a FOXA1_BDT_up.bed -b FOXA1_B14vsDMSO_nothres_up.bed > FOXA1_BDT_up.bed


bedtools intersect -wa -u -a FOXA1_BDT_up.bed -b B14.filtered.short.BAM_peaks.bed > FOXA1_BDT_up_B14ATAC.bed


bedtools intersect -wa -u -a FOXA1_BDT_up_B14ATAC.bed -b H3K27ac_B14.filtered.short.BAM.bed > FOXA1_BDT_up_B14ATAC_H3K27acB14.bed