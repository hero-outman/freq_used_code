################################################################################
# NOTE:
# ATACorrect --peaks
# ScoreBigwig --regions
# BINDetect --peaks
# For all three tools, should use merged peaks for normalization is equal
################################################################################

# ATACorrect for all 6 conditions
ncB14_naiveoverlap_bam=/Data/Projects/ATAC_seq_shFOXA1_20220302_snakepipes/noDUP_nochrM_bam/NC-B.merged.noDup.nochrM.sortpos.bam
ncDMSO_naiveoverlap_bam=/Data/Projects/ATAC_seq_shFOXA1_20220302_snakepipes/noDUP_nochrM_bam/NC-D.merged.noDup.nochrM.sortpos.bam
ncTGF_naiveoverlap_bam=/Data/Projects/ATAC_seq_shFOXA1_20220302_snakepipes/noDUP_nochrM_bam/NC-T.merged.noDup.nochrM.sortpos.bam

shB14_naiveoverlap_bam=/Data/Projects/ATAC_seq_shFOXA1_20220302_snakepipes/noDUP_nochrM_bam/F-B.merged.noDup.nochrM.sortpos.bam
shDMSO_naiveoverlap_bam=/Data/Projects/ATAC_seq_shFOXA1_20220302_snakepipes/noDUP_nochrM_bam/F-D.merged.noDup.nochrM.sortpos.bam
shTGF_naiveoverlap_bam=/Data/Projects/ATAC_seq_shFOXA1_20220302_snakepipes/noDUP_nochrM_bam/F-T.merged.noDup.nochrM.sortpos.bam

# remember to use merged peak for normalization
merged_peaks_ncshFOXA1BDT_6_conditions=/Data/Projects/ATAC_seq_shFOXA1_20220302_snakepipes/MACS_allLen_noDUP_noChrM/naiveOverlap_q001/ncFOXA1.shFOXA1.B14.DMSO.TGF.navieoverlap001.bed


hg38_path=/Data/Ref_genome_anno/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna
blacklist_path=/Data/Ref_genome_anno/ENCODE_BLACK_LIST_REGION/hg38-blacklist.v2.bed

motif_jaspar=/Data/Ref_genome_anno/motif_database/Vertebrates_singlebatch_motifs.jaspar
# motif_homocomo=/Data/Ref_genome_anno/motif_database/meme_motif_databases/HUMAN/HOCOMOCOv11_HUMAN.meme

# ncFOXA1
TOBIAS ATACorrect \
    --bam $ncB14_naiveoverlap_bam \
    --genome $hg38_path \
    --peaks $merged_peaks_ncshFOXA1BDT_6_conditions \
    --blacklist $blacklist_path \
    --outdir ./atacorrect \
    --cores 10 \
&> B14_bias_corr_july_11.log &

TOBIAS ATACorrect \
    --bam $ncDMSO_naiveoverlap_bam \
    --genome $hg38_path \
    --peaks $merged_peaks_ncshFOXA1BDT_6_conditions \
    --blacklist $blacklist_path \
    --outdir ./atacorrect \
    --cores 10 \
&> DMSO_bias_corr_july_11.log &

TOBIAS ATACorrect \
    --bam $ncTGF_naiveoverlap_bam \
    --genome $hg38_path \
    --peaks $merged_peaks_ncshFOXA1BDT_6_conditions \
    --blacklist $blacklist_path \
    --outdir ./atacorrect \
    --cores 10 \
&> TGF_bias_corr_july_11.log &

# shFOXA1
TOBIAS ATACorrect \
    --bam $shB14_naiveoverlap_bam \
    --genome $hg38_path \
    --peaks $merged_peaks_ncshFOXA1BDT_6_conditions \
    --blacklist $blacklist_path \
    --outdir ./atacorrect \
    --cores 10 \
&> shFOXA1.B14_bias_corr_july_11.log &

TOBIAS ATACorrect \
    --bam $shDMSO_naiveoverlap_bam \
    --genome $hg38_path \
    --peaks $merged_peaks_ncshFOXA1BDT_6_conditions \
    --blacklist $blacklist_path \
    --outdir ./atacorrect \
    --cores 10 \
&> shFOXA1.DMSO_bias_corr_july_11.log &

TOBIAS ATACorrect \
    --bam $shTGF_naiveoverlap_bam \
    --genome $hg38_path \
    --peaks $merged_peaks_ncshFOXA1BDT_6_conditions \
    --blacklist $blacklist_path \
    --outdir ./atacorrect \
    --cores 10 \
&>  shFOXA1.TGF_bias_corr_july_11.log &