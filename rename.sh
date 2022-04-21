rename -v fq.gz fastq.gz /Data/RAW_DATA/CUT_TAG_2021-06-02/Rawdata/*/FOXA1*fq.gz

rename -v fq.gz fastq.gz /Data/RAW_DATA/CUT_TAG_2021-06-02/Rawdata/*/H3K27*fq.gz

rename -v FOXA1_TGF_2 Foxa1TgfRep2 /Data/Projects/CUT_TAG_snakepipe/DNA_mapping/filtered_bam/FOXA1*

rename -v Foxa1B14Rep2 FOXA1_B14_2 /Data/Projects/CUT_TAG_snakepipe/DNA_mapping/filtered_bam/Foxa1*

rename -v Foxa1TgfRep1 FOXA1_TGF_1 /Data/Projects/CUT_TAG_snakepipe/DNA_mapping/filtered_bam/Foxa1*



rename -v _1.fq.gz _R1.fq.gz ./*_1.fq.gz

rename -v _2.fq.gz _R2.fq.gz ./*_2.fq.gz