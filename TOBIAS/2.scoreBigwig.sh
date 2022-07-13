################################################################################
# NOTE:
# ATACorrect --peaks
# ScoreBigwig --regions
# BINDetect --peaks
# For all three tools, should use merged peaks for normalization is equal
################################################################################

# remember to use merged peak for normalization
merged_peaks_ncshFOXA1BDT_6_conditions=/Data/Projects/ATAC_seq_shFOXA1_20220302_snakepipes/MACS_allLen_noDUP_noChrM/naiveOverlap_q001/ncFOXA1.shFOXA1.B14.DMSO.TGF.navieoverlap001.bed

hg38_path=/Data/Ref_genome_anno/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna
blacklist_path=/Data/Ref_genome_anno/ENCODE_BLACK_LIST_REGION/hg38-blacklist.v2.bed


# scorebigwig
TOBIAS ScoreBigwig \
    --signal ./atacorrect/NC-B.merged.noDup.nochrM.sortpos_corrected.bw \
    --regions $merged_peaks_ncshFOXA1BDT_6_conditions \
    --output ncB14_naiveoverlap_footprint.bw \
    --cores 10 \
&> ncB14_scorebigwig_July_11.log &

TOBIAS ScoreBigwig \
    --signal ./atacorrect/NC-D.merged.noDup.nochrM.sortpos_corrected.bw \
    --regions $merged_peaks_ncshFOXA1BDT_6_conditions \
    --output ncDMSO_naiveoverlap_footprint.bw \
    --cores 10 \
&> ncDMSO_scorebigwig_July_11.log &

TOBIAS ScoreBigwig \
    --signal ./atacorrect/NC-T.merged.noDup.nochrM.sortpos_corrected.bw \
    --regions $merged_peaks_ncshFOXA1BDT_6_conditions \
    --output ncTGF_naiveoverlap_footprint.bw \
    --cores 10 \
&> ncTGF_scorebigwig_July_11.log &


TOBIAS ScoreBigwig \
    --signal ./atacorrect/F-B.merged.noDup.nochrM.sortpos_bias.bw \
    --regions $merged_peaks_ncshFOXA1BDT_6_conditions \
    --output shB14_naiveoverlap_footprint.bw \
    --cores 10 \
&> ncB14_scorebigwig_July_11.log &

TOBIAS ScoreBigwig \
    --signal ./atacorrect/F-D.merged.noDup.nochrM.sortpos_bias.bw \
    --regions $merged_peaks_ncshFOXA1BDT_6_conditions \
    --output shDMSO_naiveoverlap_footprint.bw \
    --cores 10 \
&> ncDMSO_scorebigwig_July_11.log &

TOBIAS ScoreBigwig \
    --signal ./atacorrect/F-T.merged.noDup.nochrM.sortpos_bias.bw \
    --regions $merged_peaks_ncshFOXA1BDT_6_conditions \
    --output shTGF_naiveoverlap_footprint.bw \
    --cores 10 \
&> ncTGF_scorebigwig_July_11.log &