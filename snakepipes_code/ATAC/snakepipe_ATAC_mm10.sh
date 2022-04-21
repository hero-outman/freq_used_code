DNA_mapping_wk_dir=/Data/Projects/cw/20220214_combine_RNA_ATAC/ATAC_AB_snakepipes
ref_genome=mm10
cores_num=20

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
    -v -j $cores_num --local --DAG \
    --peakCaller MACS2 \
    --maxFragmentSize 150 --minFragmentSize 0 --qval 0.001 \
    --sampleSheet $sample_sheet \
    --FDR 0.05 --LFC 1 \
    $ref_genome  &> snake.atac_maxlength150.log
