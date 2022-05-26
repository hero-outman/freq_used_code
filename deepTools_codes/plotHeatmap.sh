plotHeatmap \
    -m FOXA1_BDT_refpoint.mat.gz \
    -out SMAD2_BDT_onepeak_test.png \
    --sortUsing sum --startLabel "Peak Start" \
    --endLabel "Peak End" \
    --xAxisLabel "" \
    --regionsLabel "peak_region" \
    --colorMap coolwarm \
    --samplesLabel "B14_SMAD2" "DMSO_SMAD2" "TGF_SMAD2"

# use kmeans: need to give same number of regionsLabel with kmeans numbers
plotHeatmap \
    -m /Users/jplab/Desktop/snakepipes_Data/CUT_TAG_snakepipes/FOXA1_H3K27ac/heatmap_profile/FOXA1_BDT.center.3k.repMerged.mat.gz \
    -out FOXA1_BDT.center.3k.kmeans4.png \
    --sortUsing sum --startLabel "Peak Start" \
    --endLabel "Peak End" \
    --xAxisLabel "" \
    --regionsLabel "peak_region_1" "peak_region_2" "peak_region_3" "peak_region_4" \
    --colorMap coolwarm \
    --kmeans 4
    # --samplesLabel "FOXA1_B14.1" "FOXA1_B14.2" "FOXA1_DMSO.1" "FOXA1_DMSO.2" "FOXA1_TGF.1" "FOXA1_TGF.2"

# with more params
plotHeatmap \
    -m /Users/jplab/Desktop/DAILY_CODE_DATA/2022-5/data/5-25_profileheatmap_motifHeatmap/shFOXA1.BDT.oc_co_po.center.3k.mat.gz \
    -out shFOXA1.BDT.oc_co_po.png \
    --sortUsing sum --startLabel "Peak Start" \
    --endLabel "Peak End" \
    --xAxisLabel "" \
    --regionsLabel "common open" "TGF.uniq.open" "B14.uniq.open" \
    --samplesLabel "shFOXA1.B14.rp1" "shFOXA1.B14.rp2" "shFOXA1.DMSO.rp1" "shFOXA1.DMSO.rp2" "shFOXA1.TGF.rp1" "shFOXA1.TGF.rp2" \
    --colorMap coolwarm \
    --heatmapHeight 50


plotHeatmap \
    -m /Users/jplab/Desktop/DAILY_CODE_DATA/2022-5/data/5-25_profileheatmap_motifHeatmap/ncFOXA1.BDT.oc_co_po.center.3k.mat.gz \
    -out ncFOXA1.BDT.oc_co_po.png \
    --sortUsing sum --startLabel "Peak Start" \
    --endLabel "Peak End" \
    --xAxisLabel "" \
    --regionsLabel "common open" "TGF.uniq.open" "B14.uniq.open" \
    --samplesLabel "ncFOXA1.B14.rp1" "ncFOXA1.B14.rp2" "ncFOXA1.DMSO.rp1" "ncFOXA1.DMSO.rp2" "ncFOXA1.TGF.rp1" "ncFOXA1.TGF.rp2" \
    --colorMap coolwarm \
    --heatmapHeight 50    
