# input have to be gzipped: fastq.gz; or will throw error!
# use pigz for compression, gzip too slow! 
# genome need created and replaced in the config file: use snakepipes info to show
# issue: defalut config will create a temp folder in /data/temp/ , but there is no such path, so the flow will throw error and stop work
# i update default config change temp folder to ./temp/ ; so, need create temp folder in the workding dir!
# this snakepipe is ugly unstable, so better download error log and use vscode to check;
# todo: use snakemake to make a pipeline myself;

# need to create this folder to solve snakepipe running issue
work_dir=/Data/Projects/ATAC_seq_shFOXA1_20220302_snakepipes
mkdir -p $work_dir'/temp/'
echo 'folder '$work_dir'/temp/ was created!'

DNA-mapping \
    -v --DAG -j 20 --local \
    -i $work_dir'/fastq/' -o $work_dir \
    --ext ".fq.gz" \
    --trim --trimmer cutadapt --trimmerOptions '-a CTGTCTCTTATACACATCTCCGAGCCCACGAGAC -A CTGTCTCTTATACACATCTGACGCTGCCGACGA' \
    --fastqc \
    --bwBinSize '1' \
    --plotFormat 'pdf' \
    --aligner 'Bowtie2' \
    --alignerOpts '-p 20 -t --mm --very-sensitive --no-unal --no-mixed --no-discordant -X2000' \
    --qualimap \
    --dedup \
    --mapq '3' \
    --insertSizeMax '2000' \
    --properPairs \
    hg38 &> snakepipe_DNAmapping.out &