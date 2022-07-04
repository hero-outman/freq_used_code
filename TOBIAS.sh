# ATACorrect for all 3 conditions
ncB14_naiveoverlap_bam=/Data/Projects/ATAC_seq_shFOXA1_20220302_snakepipes/noDUP_nochrM_bam/NC-B.merged.noDup.nochrM.sortpos.bam
ncDMSO_naiveoverlap_bam=/Data/Projects/ATAC_seq_shFOXA1_20220302_snakepipes/noDUP_nochrM_bam/NC-D.merged.noDup.nochrM.sortpos.bam
ncTGF_naiveoverlap_bam=/Data/Projects/ATAC_seq_shFOXA1_20220302_snakepipes/noDUP_nochrM_bam/NC-T.merged.noDup.nochrM.sortpos.bam

# merge B D T 3 conditions peak file
# merge H3K27ac B14 DMSO, but not TGF(use TGF rep2, rep1 is bad enriched) peak for deeptools computematrix
cat /Data/Projects/ATAC_seq_shFOXA1_20220302_snakepipes/MACS_allLen_noDUP_noChrM/naiveOverlap_q001/NC_B_overlap_peaks.narrowPeak \
 /Data/Projects/ATAC_seq_shFOXA1_20220302_snakepipes/MACS_allLen_noDUP_noChrM/naiveOverlap_q001/NC_D_overlap_peaks.narrowPeak \
 /Data/Projects/ATAC_seq_shFOXA1_20220302_snakepipes/MACS_allLen_noDUP_noChrM/naiveOverlap_q001/NC_T_overlap_peaks.narrowPeak \
 | sort -k1,1 -k2,2n | mergeBed -i stdin -o collapse -c 4 > ncB14.DMSO.TGF.navieoverlap.merged.bed &

merged_bed=/Data/Projects/ATAC_seq_shFOXA1_20220302_snakepipes/TOBIAS/ncB14.DMSO.TGF.navieoverlap.merged.bed

hg38_path=/Data/Ref_genome_anno/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna
blacklist_path=/Data/Ref_genome_anno/ENCODE_BLACK_LIST_REGION/hg38-blacklist.v2.bed

motif_jaspar=/Data/Ref_genome_anno/motif_database/Vertebrates_singlebatch_motifs.jaspar
motif_homocomo=/Data/Ref_genome_anno/motif_database/meme_motif_databases/HUMAN/HOCOMOCOv11_HUMAN.meme

# for B14
TOBIAS ATACorrect \
    --bam $ncB14_naiveoverlap_bam \
    --genome $hg38_path \
    --peaks $merged_bed \
    --blacklist $blacklist_path \
    --outdir ./atacorrect \
    --cores 20 \
&> B14_bias_corr_6_10.log &

TOBIAS ATACorrect \
    --bam $ncDMSO_naiveoverlap_bam \
    --genome $hg38_path \
    --peaks $merged_bed \
    --blacklist $blacklist_path \
    --outdir ./atacorrect \
    --cores 20 \
&> DMSO_bias_corr_6_10.log &

TOBIAS ATACorrect \
    --bam $ncTGF_naiveoverlap_bam \
    --genome $hg38_path \
    --peaks $merged_bed \
    --blacklist $blacklist_path \
    --outdir ./atacorrect \
    --cores 20 \
&> TGF_bias_corr_6_10.log &



# scorebigwig
# for B14 
TOBIAS ScoreBigwig \
    --signal ./atacorrect/NC-B.merged.noDup.nochrM.sortpos_corrected.bw \
    --regions $merged_bed \
    --output ncB14_naiveoverlap_footprint.bw \
    --cores 20 \
&> ncB14_scorebigwig_6_10.log &

TOBIAS ScoreBigwig \
    --signal ./atacorrect/NC-D.merged.noDup.nochrM.sortpos_corrected.bw \
    --regions $merged_bed \
    --output ncDMSO_naiveoverlap_footprint.bw \
    --cores 20 \
&> ncDMSO_scorebigwig_6_10.log &

TOBIAS ScoreBigwig \
    --signal ./atacorrect/NC-T.merged.noDup.nochrM.sortpos_corrected.bw \
    --regions $merged_bed \
    --output ncTGF_naiveoverlap_footprint.bw \
    --cores 20 \
&> ncTGF_scorebigwig_6_10.log &


# need to remove KI or GI from peak file!!! or will throw error and stop running
# Bindetect
# B14 vs TGFb MEME database
TOBIAS BINDetect \
    --motifs $motif_jaspar \
    --signals \
        ./scorebigwig/ncB14_naiveoverlap_footprint.bw \
        ./scorebigwig/ncDMSO_naiveoverlap_footprint.bw \
        ./scorebigwig/ncTGF_naiveoverlap_footprint.bw \
    --genome $hg38_path \
    --peaks $merged_bed \
    --outdir ./bindetect \
    --cond_names ncB14 ncDMSO ncTGF \
    --cores 20 \
&> bindetect_meme_6-6.log &

#########################################################################################################
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