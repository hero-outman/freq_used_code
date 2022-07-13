# remember to use merged peak for normalization
merged_peaks_ncshFOXA1BDT_6_conditions=/Data/Projects/ATAC_seq_shFOXA1_20220302_snakepipes/MACS_allLen_noDUP_noChrM/naiveOverlap_q001/ncFOXA1.shFOXA1.B14.DMSO.TGF.navieoverlap001.bed

hg38_path=/Data/Ref_genome_anno/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna
blacklist_path=/Data/Ref_genome_anno/ENCODE_BLACK_LIST_REGION/hg38-blacklist.v2.bed

motif_jaspar=/Data/Ref_genome_anno/motif_database/Vertebrates_singlebatch_motifs.jaspar
# motif_homocomo=/Data/Ref_genome_anno/motif_database/meme_motif_databases/HUMAN/HOCOMOCOv11_HUMAN.meme

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
