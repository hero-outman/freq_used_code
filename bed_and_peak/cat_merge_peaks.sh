# merge FOXA1 B14 DMSO, but not TGF(use TGF rep2, rep1 is bad enriched) peak for deeptools computematrix
peak_folder=/Data/Projects/CUT_TAG_snakepipe/DNA_mapping/MACS2/
awk '{print $0"\t","FOXA1_B14_1-"NR}' $peak_folder'FOXA1_B14_1.filtered.short.BAM_peaks.narrowPeak' > FOXA1_B14_1.filtered.short.BAM_peaks.bed.temp && \
awk '{print $0"\t","FOXA1_B14_2-"NR}' $peak_folder'FOXA1_B14_2.filtered.short.BAM_peaks.narrowPeak' > FOXA1_B14_2.filtered.short.BAM_peaks.bed.temp && \
awk '{print $0"\t","FOXA1_DMSO_1-"NR}' $peak_folder'FOXA1_DMSO_1.filtered.short.BAM_peaks.narrowPeak' > FOXA1_DMSO_1.filtered.short.BAM_peaks.bed.temp && \
awk '{print $0"\t","FOXA1_DMSO_2-"NR}' $peak_folder'FOXA1_DMSO_2.filtered.short.BAM_peaks.narrowPeak' > FOXA1_DMSO_2.filtered.short.BAM_peaks.bed.temp && \
awk '{print $0"\t","FOXA1_TGF-"NR}' $peak_folder'FOXA1_TGF_2.filtered.short.BAM_peaks.narrowPeak' > FOXA1_TGF_2.filtered.short.BAM_peaks.bed.temp && \
cat FOXA1_B14_1.filtered.short.BAM_peaks.bed.temp FOXA1_B14_2.filtered.short.BAM_peaks.bed.temp FOXA1_DMSO_1.filtered.short.BAM_peaks.bed.temp FOXA1_DMSO_2.filtered.short.BAM_peaks.bed.temp FOXA1_TGF_2.filtered.short.BAM_peaks.bed.temp | sort -k1,1 -k2,2n | mergeBed -i stdin -o collapse -c 4 > FOXA1_B14.DMSO.TGF_peaks.bed && \
rm FOXA1_B14_1.filtered.short.BAM_peaks.bed.temp FOXA1_B14_2.filtered.short.BAM_peaks.bed.temp FOXA1_DMSO_1.filtered.short.BAM_peaks.bed.temp FOXA1_DMSO_2.filtered.short.BAM_peaks.bed.temp FOXA1_TGF_2.filtered.short.BAM_peaks.bed.temp


# merge H3K27ac B14 DMSO, but not TGF(use TGF rep2, rep1 is bad enriched) peak for deeptools computematrix
peak_folder=/Data/Projects/CUT_TAG_snakepipe/DNA_mapping/MACS2/
awk '{print $0"\t","H3K27ac_B14_1-"NR}' $peak_folder'H3K27ac_B14_1.filtered.short.BAM_peaks.narrowPeak' > H3K27ac_B14_1.filtered.short.BAM_peaks.bed.temp && \
awk '{print $0"\t","H3K27ac_B14_2-"NR}' $peak_folder'H3K27ac_B14_2.filtered.short.BAM_peaks.narrowPeak' > H3K27ac_B14_2.filtered.short.BAM_peaks.bed.temp && \
awk '{print $0"\t","H3K27ac_DMSO_1-"NR}' $peak_folder'H3K27ac_DMSO_1.filtered.short.BAM_peaks.narrowPeak' > H3K27ac_DMSO_1.filtered.short.BAM_peaks.bed.temp && \
awk '{print $0"\t","H3K27ac_DMSO_2-"NR}' $peak_folder'H3K27ac_DMSO_2.filtered.short.BAM_peaks.narrowPeak' > H3K27ac_DMSO_2.filtered.short.BAM_peaks.bed.temp && \
awk '{print $0"\t","H3K27ac_TGF_1-"NR}' $peak_folder'H3K27ac_TGF_1.filtered.short.BAM_peaks.narrowPeak' > H3K27ac_TGF_1.filtered.short.BAM_peaks.bed.temp && \
awk '{print $0"\t","H3K27ac_TGF_2-"NR}' $peak_folder'H3K27ac_TGF_2.filtered.short.BAM_peaks.narrowPeak' > H3K27ac_TGF_2.filtered.short.BAM_peaks.bed.temp && \
cat H3K27ac_B14_1.filtered.short.BAM_peaks.bed.temp H3K27ac_B14_2.filtered.short.BAM_peaks.bed.temp H3K27ac_DMSO_1.filtered.short.BAM_peaks.bed.temp H3K27ac_DMSO_2.filtered.short.BAM_peaks.bed.temp H3K27ac_TGF_1.filtered.short.BAM_peaks.bed.temp H3K27ac_TGF_2.filtered.short.BAM_peaks.bed.temp | sort -k1,1 -k2,2n | mergeBed -i stdin -o collapse -c 4 > H3K27ac_B14.DMSO.TGF_peaks.bed && \
rm H3K27ac_B14_1.filtered.short.BAM_peaks.bed.temp H3K27ac_B14_2.filtered.short.BAM_peaks.bed.temp H3K27ac_DMSO_1.filtered.short.BAM_peaks.bed.temp H3K27ac_DMSO_2.filtered.short.BAM_peaks.bed.temp H3K27ac_TGF_1.filtered.short.BAM_peaks.bed.temp H3K27ac_TGF_2.filtered.short.BAM_peaks.bed.temp
