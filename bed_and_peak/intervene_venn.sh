# check overlap between old B14 cid1 and 
intervene venn -i /Users/jplab/Desktop/DAILY_CODE_DATA/2021-05_codes/5-10/redefine/clus/cid_1.bed /Users/jplab/Desktop/DAILY_CODE_DATA/2022-6/data/6-2_pybedtools_B14TF/ncFOXA1.BvsDnoDnoT.Bup.bed --B14_uniq_bed_overlap.png

# check cid overlap with 



intervene venn -i /Users/jplab/Desktop/DAILY_CODE_DATA/2021-05_codes/5-10/redefine/clus/cid_1.bed /Users/jplab/Desktop/snakepipes_Data/ATAC_seq_shFOXA1_snakepipes_results/CSAW_MACS2/ncFOXA1_BDT/bed/oc_co_po_bed/ncFOXA1_BDT_B14UniqOpen.bed --filenames

intervene venn -i \
/Users/jplab/Desktop/snakepipes_Data/ATAC_seq_shFOXA1_snakepipes_results/MACS_allLen_noDUP_noChrM/naiveOverlap_q001/NC_T_overlap_peaks.sorted.narrowPeak \
/Users/jplab/Desktop/snakepipes_Data/ATAC_seq_shFOXA1_snakepipes_results/MACS_allLen_noDUP_noChrM/naiveOverlap_q001/NC_D_overlap_peaks.sorted.narrowPeak \
/Users/jplab/Desktop/snakepipes_Data/ATAC_seq_shFOXA1_snakepipes_results/MACS_allLen_noDUP_noChrM/naiveOverlap_q001/NC_B_overlap_peaks.sorted.narrowPeak \
--filenames

intervene venn -i \
/Users/jplab/Desktop/snakepipes_Data/ATAC_seq_shFOXA1_snakepipes_results/MACS_allLen_noDUP_noChrM/naiveOverlap_q001/NC_T_overlap_peaks.sorted.narrowPeak \
/Users/jplab/Desktop/snakepipes_Data/ATAC_seq_shFOXA1_snakepipes_results/MACS_allLen_noDUP_noChrM/naiveOverlap_q001/NC_D_overlap_peaks.sorted.narrowPeak \
/Users/jplab/Desktop/DAILY_CODE_DATA/2022-6/data/6-9_peak_venn/beds/ncFOXA1.BvsDnoT.Bup.bed \
--filenames

intervene venn -i \
/Users/jplab/Desktop/snakepipes_Data/ATAC_seq_shFOXA1_snakepipes_results/MACS_allLen_noDUP_noChrM/naiveOverlap_q001/NC_T_overlap_peaks.sorted.narrowPeak \
/Users/jplab/Desktop/snakepipes_Data/ATAC_seq_shFOXA1_snakepipes_results/MACS_allLen_noDUP_noChrM/naiveOverlap_q001/NC_D_overlap_peaks.sorted.narrowPeak \
/Users/jplab/Desktop/DAILY_CODE_DATA/2022-6/data/6-9_peak_venn/beds/ncFOXA1.navieOverlap.B.NoD.NoT.bed \
--filenames

intervene venn -i \
/Users/jplab/Desktop/snakepipes_Data/ATAC_seq_shFOXA1_snakepipes_results/MACS_allLen_noDUP_noChrM/naiveOverlap_q001/NC_T_overlap_peaks.sorted.narrowPeak \
/Users/jplab/Desktop/snakepipes_Data/ATAC_seq_shFOXA1_snakepipes_results/MACS_allLen_noDUP_noChrM/naiveOverlap_q001/NC_D_overlap_peaks.sorted.narrowPeak \
/Users/jplab/Desktop/DAILY_CODE_DATA/2022-6/data/6-9_peak_venn/beds/ncFOXA1.BvsDnoDnoT.Bup.bed \
--filenames