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


motif_jaspar=/Data/Ref_genome_anno/motif_database/Vertebrates_singlebatch_motifs.jaspar
motif_homocomo=/Data/Ref_genome_anno/motif_database/meme_motif_databases/HUMAN/HOCOMOCOv11_HUMAN.meme

####################################
# ATACorrect
####################################

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


####################################
# ScoreBigwig
####################################

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


####################################
# BINDetect
# need to remove KI or GI from peak file!!! or will throw error and stop running
####################################
TOBIAS BINDetect \
    --motifs $motif_jaspar \
    --signals \
        ./scorebigwig/ncB14_naiveoverlap_footprint.bw \
        ./scorebigwig/ncDMSO_naiveoverlap_footprint.bw \
        ./scorebigwig/ncTGF_naiveoverlap_footprint.bw \
        ./scorebigwig/shB14_naiveoverlap_footprint.bw \
        ./scorebigwig/shDMSO_naiveoverlap_footprint.bw \
        ./scorebigwig/shTGF_naiveoverlap_footprint.bw \
    --genome $hg38_path \
    --peaks $merged_peaks_ncshFOXA1BDT_6_conditions \
    --outdir ./bindetect_shANDnc \
    --cond_names ncB14 ncDMSO ncTGF shB14 shDMSO shTGF \
    --cores 40 \
&> bindetect_meme_July-11.log &

#########################################################################################################
# PlotAggregate
# attention: do not use region.bed to constrain PlotAggregate, this will cause weaker signal of footprints
# if use, need to think carefully
# or choose *_all.bed to get stronger footprint
#########################################################################################################

# plot aggregate
# b14 uniq open region
TOBIAS PlotAggregate \
    --TFBS \
        ./bindetect/ELF3*/beds/*B14_bound.bed \
        ./bindetect/ELF5*/beds/*DMSO_bound.bed \
        ./bindetect/ELF5*/beds/*TGF_bound.bed \
    --signals \
        ./atacorrect/NC-B.merged.noDup.nochrM.sortpos_corrected.bw \
        ./atacorrect/NC-D.merged.noDup.nochrM.sortpos_corrected.bw \
        ./atacorrect/NC-T.merged.noDup.nochrM.sortpos_corrected.bw \
    --output \
        ncB14.DMSO.TGF.navieoverlap.merged.pdf \
    --regions \
        /Data/Projects/ATAC_seq_2020_12/2umol_start_from_raw/peak/macs_peak/clus/cid_1.bed \
        /Data/Projects/ATAC_seq_2020_12/2umol_start_from_raw/peak/macs_peak/clus/cid_3.bed \
    --share_y both \
    --plot_boundaries



TOBIAS PlotAggregate \
    --TFBS \
        ./bindetect/jaspar/FOS_*/beds/FOS*_all.bed \
        ./bindetect/jaspar/FOSJUND_*/beds/FOSJUND*_all.bed \
        ./bindetect/jaspar/JUNvar.2_*/beds/JUNvar.2*_all.bed \
        ./bindetect/jaspar/FOSL1JUN_*/beds/FOSL1JUN*_all.bed \
        ./bindetect/jaspar/Smad2Smad3_*/beds/Smad2Smad3*_all.bed \
        ./bindetect/jaspar/FOSL1_*/beds/FOSL1*_all.bed \
        ./bindetect/jaspar/BATF3_*/beds/BATF3*_all.bed \
        ./bindetect/jaspar/FOSL2_*/beds/FOSL2*_all.bed \
        ./bindetect/jaspar/FOSL2JUN_*/beds/FOSL2JUN*_all.bed \
        ./bindetect/jaspar/BATFJUN_*/beds/BATFJUN*_all.bed \
        ./bindetect/jaspar/BATF_*/beds/BATF*_all.bed \
        ./bindetect/jaspar/FOSL1JUNB_*/beds/FOSL1JUNB*_all.bed \
        ./bindetect/jaspar/FOSL2JUNB_*/beds/FOSL2JUNB*_all.bed \
        ./bindetect/jaspar/FOSBJUNB_*/beds/FOSBJUNB*_all.bed \
        ./bindetect/jaspar/FOSJUN_*/beds/FOSJUN*_all.bed \
        ./bindetect/jaspar/FOSL2JUND_*/beds/FOSL2JUND*_all.bed \
        ./bindetect/jaspar/JUNB_*/beds/JUNB*_all.bed \
        ./bindetect/jaspar/FOSJUNB_*/beds/FOSJUNB*_all.bed \
        ./bindetect/jaspar/JUND_*/beds/JUND*_all.bed \
        ./bindetect/jaspar/FOSL1JUND_*/beds/FOSL1JUND*_all.bed \
        ./bindetect/jaspar/JUNJUNB_*/beds/JUNJUNB*_all.bed \
    --TFBS-labels \
        FOS \
        FOSJUND \
        JUNvar.2 \
        FOSL1JUN \
        Smad2Smad3 \
        FOSL1 \
        BATF3 \
        FOSL2 \
        FOSL2JUN \
        BATFJUN \
        BATF \
        FOSL1JUNB \
        FOSL2JUNB \
        FOSBJUNB \
        FOSJUN \
        FOSL2JUND \
        JUNB \
        FOSJUNB \
        JUND \
        FOSL1JUND \
        JUNJUNB \
    --signals \
        ./atacorrect/NC-B.merged.noDup.nochrM.sortpos_corrected.bw \
        ./atacorrect/NC-D.merged.noDup.nochrM.sortpos_corrected.bw \
        ./atacorrect/NC-T.merged.noDup.nochrM.sortpos_corrected.bw \
    --signal-labels \
        ncFOXA1.B1.14 \
        ncFOXA1.DMSO \
        ncFOXA1.TGF-beta \
    --output \
        ncB14.DMSO.TGF.navieoverlap.merged_B14closed_top20.png \
    --share_y both \
    --plot_boundaries