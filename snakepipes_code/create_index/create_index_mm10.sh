createIndices \
    -j 30 \
    --DAG \
    -o /Data/Ref_genome_anno/Mus_musculus/mm10/snakepipes_index \
    --local \
    --genome /Data/Ref_genome_anno/Mus_musculus/mm10/Mus_musculus_Ensemble_96.fa.gz \
    --gtf /Data/Ref_genome_anno/Mus_musculus/mm10/Mus_musculus_Ensemble_96.gtf.gz \
    --blacklist /Data/Ref_genome_anno/Mus_musculus/mm10/mm10.blacklist.bed \
    --rmskURL /Data/Ref_genome_anno/Mus_musculus/mm10/rmsk_mm10.txt \
    mm10.96 &> snakepipes_createindex.log &