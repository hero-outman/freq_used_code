# genome need created and replaced in the config file: use snakepipes info to show
# issue: defalut config will create a temp folder in /data/temp/ , but there is no such path, so the flow will throw error and stop work
# i update default config change temp folder to ./temp/ ; so, need create temp folder in the workding dir!
# this snakepipe is ugly unstable, so better download error log and use vscode to check;
# todo: use snakemake to make a pipeline myself;

# check if ./temp exists, if not, create it
if [ ! -d ./temp ]; then
    mkdir ./temp
fi

sample_sheet=$1
if [ -z "$sample_sheet" ]; then
    echo "NO sample sheet provided, differentail expression analysis will not be performed!"
    
    mRNA-seq \
    -v \
    --local \
    -j 20 \
    --ext '.fq.gz' --reads '_R1' '_R2' \
    --fastqc \
    --trim --trimmer cutadapt --trimmerOptions '-a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT -m 20 -q 20' \
    --aligner HISAT2 --alignerOptions '-p 20 -t --mm --no-unal --no-mixed --no-discordant' \
    --featureCountsOptions '-T 8 -t exon -g gene_id -p -B -C' \
    -i /Data/Projects/mRNA_seq_FOXA1_KD_NC_2022221_snakepipe/fastq \
    -o /Data/Projects/mRNA_seq_FOXA1_KD_NC_2022221_snakepipe/ \
    --dnaContam \
    hg38 &> snakepipe_mRNA.out &
fi

if [ ! -z "$sample_sheet" ]; then
    echo "Differentail expression analysis will be performed with sample sheet: $sample_sheet !!!"
    
    mRNA-seq \
    -v \
    --local \
    -j 20 \
    --ext '.fq.gz' --reads '_R1' '_R2' \
    --fastqc \
    --trim --trimmer cutadapt --trimmerOptions '-a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT -m 20 -q 20' \
    --aligner HISAT2 --alignerOptions '-p 20 -t --mm --no-unal --no-mixed --no-discordant' \
    --featureCountsOptions '-T 8 -t exon -g gene_id -p -B -C' \
    -i /Data/Projects/mRNA_seq_FOXA1_KD_NC_2022221_snakepipe/fastq \
    -o /Data/Projects/mRNA_seq_FOXA1_KD_NC_2022221_snakepipe/ \
    --sampleSheet $sample_sheet \
    --dnaContam \
    hg38 &> snakepipe_mRNA.out &
fi