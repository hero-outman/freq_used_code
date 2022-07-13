#########################################################################################################
# PlotAggregate 
#
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