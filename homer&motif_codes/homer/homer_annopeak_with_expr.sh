# new clus bed got from hogwarts
# use homer to bind expr data
annotatePeaks.pl  /Users/jplab/Desktop/hogwarts_atac_codes/redefined/clus/cid_1.bed \
    hg38 \
    -gene /Users/jplab/Desktop/Projects/RNA_seq_workflow/downstream/deg_2021/DEG_B14_TGFb_ALL_filtered_expr.tsv \
    -m /Users/jplab/Desktop/hogwarts_atac_codes/downtream/anno_homer/interested_motif/FOXM1.motif > B14_open_foxm1_anno_expr.tsv 

# no expr data
# with motif file to show binding sequence is or not in peak region
annotatePeaks.pl /Users/jplab/Desktop/2022-2/data/2-16_CSAW/FOXA1_BDT_up_B14ATAC_H3K27acB14.bed \
    hg38 \
    -m /Users/jplab/Desktop/hogwarts_atac_codes/downtream/anno_homer/interested_motif/FOXM1.motif /Users/jplab/Downloads/FOXA1.jaspar /Users/jplab/Downloads/CTCF.jaspar \
    -nmotifs \
    -mbed FOXM1_FOXA1_CTCF.bed \
    -mask \
    -annStats anno_info.tsv > FOXA1_BDT_up_B14ATAC_H3K27acB14.tsv

# only annotate peak region, no motif, no expr
annotatePeaks.pl -mask /Users/jplab/Desktop/hogwarts_atac_codes/redefined/clus/cid_1.bed hg38 > B14_open_foxm1_num.tsv