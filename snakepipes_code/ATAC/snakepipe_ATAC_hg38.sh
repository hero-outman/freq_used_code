DNA_mapping_wk_dir=/Data/Projects/ATAC_seq_shFOXA1_20220302_snakepipes
ref_genome=hg38
minFragLength=20
maxFragLength=2000

sample_sheet=$1
if [ -z "$sample_sheet" ]; then
    echo "please provide sample sheet"
    exit 1
fi

start_from_filteredBam=$2
if [ -z "$start_from_filteredBam" ]; then
    echo "start_from_filteredBam file, yes or no"
    exit 1
fi

work_dir=./
temp_folder=$work_dir'temp/'
mkdir -p $temp_folder
echo 'folder '$temp_folder' was created!'

ATAC-seq \
    -d $DNA_mapping_wk_dir \
    -v -j 20 --local --DAG \
    --peakCaller MACS2 \
    --maxFragmentSize $maxFragLength --minFragmentSize $minFragLength --qval 0.001 \
    --sampleSheet $sample_sheet \
    --FDR 0.05 --LFC 1 \
    $ref_genome  &> 'snake.atac_lengthmax'$maxFragLength'min'$minFragLength.log
