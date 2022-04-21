merged_peak_folder=/Data/Projects/CUT_TAG_snakepipe/DNA_mapping/MACS2/merged_peaks/
finalbigwig_folder=/Data/Projects/CUT_TAG_snakepipe/DNA_mapping/deepTools_ATAC/bamCompare/

computeMatrix scale-regions \
  -R $merged_peak_folder'FOXA1_B14.DMSO.TGF_peaks.bed' \
  -S $finalbigwig_folder'FOXA1_B14_1.filtered.bw' \
    $finalbigwig_folder'FOXA1_B14_2.filtered.bw' \
    $finalbigwig_folder'FOXA1_DMSO_1.filtered.bw' \
    $finalbigwig_folder'FOXA1_DMSO_2.filtered.bw' \
    $finalbigwig_folder'FOXA1_TGF_2.filtered.bw' \
  --skipZeros \
  -bl /Data/Ref_genome_anno/ENCODE_BLACK_LIST_REGION/hg38-blacklist.v2.bed \
  --missingDataAsZero \
  -o FOXA1_BDT_mat.gz \
  -a 5000 -b 5000 \
  --outFileNameMatrix FOXA1_BDT_mat.tab \
  --outFileSortedRegions FOXA1_BDT_mat.bed \
  -p 20

computeMatrix scale-regions \
  -R $merged_peak_folder'H3K27ac_B14.DMSO.TGF_peaks.bed' \
  -S $finalbigwig_folder'H3K27ac_B14_1.filtered.bw' \
    $finalbigwig_folder'H3K27ac_B14_2.filtered.bw' \
    $finalbigwig_folder'H3K27ac_DMSO_1.filtered.bw' \
    $finalbigwig_folder'H3K27ac_DMSO_2.filtered.bw' \
    $finalbigwig_folder'H3K27ac_TGF_1.filtered.bw' \
    $finalbigwig_folder'H3K27ac_TGF_2.filtered.bw' \
  --skipZeros \
  -bl /Data/Ref_genome_anno/ENCODE_BLACK_LIST_REGION/hg38-blacklist.v2.bed \
  --missingDataAsZero \
  -o H3K27ac_BDT.mat.gz \
  -a 5000 -b 5000 \
  --outFileNameMatrix H3K27ac_BDT.tab \
  --outFileSortedRegions H3K27ac_BDT.bed \
  -p 20