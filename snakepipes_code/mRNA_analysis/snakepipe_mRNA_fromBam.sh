# check if ./temp exists, if not, create it 
if [ ! -d ./temp ]; then
    mkdir ./temp
fi

which_hg=hg38
# which_hg=hs37d5
nohup mRNA-seq \
    -v \
    --snakemakeOptions='--conda-prefix /home/sunchu/.conda/envs/snakePipes/lib/python3.9/site-packages/snakePipes/shared/rules/envs/' \
    --local \
    -j 20 \
    --featureCountsOptions '-T 8 -t exon -g gene_id -p -B -C' \
    --sampleSheet '/Data/Projects/RNA_seq_FOXA1KO_scramble_snakepipe/sample_sheet.tsv' \
    -i /Data/Projects/RNA_seq_FOXA1KO_scramble_20211213/bam \
    -o /Data/Projects/RNA_seq_FOXA1KO_scramble_snakepipe \
    --fromBAM \
    --bamExt .sortpos.bam \
    $which_hg &> snakepipe_mRNA_fromBam.out &