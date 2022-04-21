merged_peak_folder=/Data/Projects/CUT_TAG_snakepipe/DNA_mapping/MACS2/merged_peaks/
finalbigwig_folder=/Data/Projects/CUT_TAG_snakepipe/DNA_mapping/deepTools_ATAC/bamCompare/

computeMatrix reference-point \
  -R $merged_peak_folder'H3K27ac_B14.DMSO.TGF_peaks.bed' \
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
  -o H3K27ac_BDT.center.3k.mat.gz \
  -a 3000 -b 3000 \
  --outFileNameMatrix H3K27ac_BDT.center.3k..tab \
  --outFileSortedRegions H3K27ac_BDT.center.3k.bed \
  -p 20

computeMatrix reference-point \
  -R $merged_peak_folder'H3K27ac_B14.DMSO.TGF_peaks.bed' \
  -S $finalbigwig_folder'H3K27ac_B14_1.filtered.bw' \
    $finalbigwig_folder'H3K27ac_B14_2.filtered.bw' \
    $finalbigwig_folder'H3K27ac_DMSO_1.filtered.bw' \
    $finalbigwig_folder'H3K27ac_DMSO_2.filtered.bw' \
    $finalbigwig_folder'H3K27ac_TGF_1.filtered.bw' \
    $finalbigwig_folder'H3K27ac_TGF_2.filtered.bw' \
  --referencePoint TSS \
  --skipZeros \
  -bl /Data/Ref_genome_anno/ENCODE_BLACK_LIST_REGION/hg38-blacklist.v2.bed \
  --missingDataAsZero \
  -o H3K27ac_BDT.TSS.3k.mat.gz \
  -a 3000 -b 3000 \
  --outFileNameMatrix H3K27ac_BDT.TSS.3k..tab \
  --outFileSortedRegions H3K27ac_BDT.TSS.3k.bed \
  -p 20
