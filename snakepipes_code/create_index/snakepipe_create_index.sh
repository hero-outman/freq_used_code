# !/bin/bash

# snakepipe hg19 create index


# -v voerbose
# --local run workflow locally; default: jobs are submitted to Slurm queue (default: 'False') 
# Genome: name to save as prefix for index files, hg19 or GRCh38 all ok
createIndices \
    -o /Data/Ref_genome_anno/Homo_sapiens/hg19/snakepipes_index \
    --local \
    --genome /Data/Ref_genome_anno/Homo_sapiens/hg19/hg19_fa_and_gtf_file_20211207/hg19.p13.plusMT.no_alt_analysis_set.fa.gz \
    --gtf /Data/Ref_genome_anno/Homo_sapiens/hg19/hg19_fa_and_gtf_file_20211207/gtf/hg19.ensGene.gtf.gz \
    --blacklist /Data/Ref_genome_anno/ENCODE_BLACK_LIST_REGION/hg38-blacklist.v2.bed \
    hg19 &> snakepipes_createindex.log &