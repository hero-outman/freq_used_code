merged_peak_folder=/Data/Projects/CUT_TAG_snakepipe/DNA_mapping/MACS2/merged_peaks/
finalbigwig_folder=/Data/Projects/CUT_TAG_snakepipe/DNA_mapping/deepTools_ATAC/bamCompare/

peak_region=H3K27ac_B14.DMSO.TGF_peaks.bed
sample_name=`basename $peak_region _peaks.bed`

cores=20

computeMatrix reference-point \
  -R $merged_peak_folder${peak_region} \
  -S $finalbigwig_folder'H3K27ac_B14_1.filtered.bw' \
    $finalbigwig_folder'H3K27ac_B14_2.filtered.bw' \
    $finalbigwig_folder'H3K27ac_DMSO_1.filtered.bw' \
    $finalbigwig_folder'H3K27ac_DMSO_2.filtered.bw' \
    $finalbigwig_folder'H3K27ac_TGF_1.filtered.bw' \
    $finalbigwig_folder'H3K27ac_TGF_2.filtered.bw' \
  --referencePoint center \
  --skipZeros \
  -bl /Data/Ref_genome_anno/ENCODE_BLACK_LIST_REGION/hg38-blacklist.v2.bed \
  --missingDataAsZero \
  -o ${sample_name}.center.3k.mat.gz \
  -a 3000 -b 3000 \
  --outFileNameMatrix ${sample_name}.center.3k..tab \
  --outFileSortedRegions ${sample_name}.center.3k.bed \
  -p ${cores}
