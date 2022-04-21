##############################################################################################################################

# 2021-08-27

##############################################################################################################################

blacklist_path=/Data/Ref_genome_anno/ENCODE_BLACK_LIST_REGION/hg38-blacklist.v2.bed
hg38_path=/Data/Ref_genome_anno/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna
hg38_chrom_sizes=/Data/Ref_genome_anno/hg38.chrom.sizes


# 0.build link
ln -s /Data/RAW_DATA/CUT_TAG_2021-06-02/Rawdata/*/*gz /Data/Projects/CUT_TAG_2021-06-02/fastq

# 1. qc on rawdata
fastqc -t 20  */*gz -o /Users/jplab/Desktop/cut_tag_data/ANNO_XS05KF2019110639_PM-XS05KF2019110639-10_BHFMK3CCX2_2021-05-30_01-50-52/qc/raw_qc

# 2. cut adaptor and low quality data
####################
trimmomatic.sh
####################
echo 'begin to cut adaptor of ' $1 && \
java -jar /Data/Tools/Trimmomatic-0.39/trimmomatic-0.39.jar \
PE -phred33 \
-threads 10 \
$2 $3 $1'_cut_R1.fq.gz' $1'_cut.unpaired_R1.fq.gz' $1'_cut_R2.fq.gz' $1'_cut.unpaired_R2.fq.gz' \
ILLUMINACLIP:/Data/Tools/Trimmomatic-0.39/adapters/NexteraPE-PE.fa:2:30:10:2:true LEADING:3 SLIDINGWINDOW:4:15 MINLEN:30 \
&& echo $1 ' finished' &

####################
cut_adaptor_loop.sh
####################
for f in  ./*.fq.gz
do
    root=`basename $f`

    if [[ $root == *R1.fq* ]]
    then # PE
	bf=`echo $root | sed -r 's/_R1.fq//g' | sed -r s/.gz//g` # control_a
        tt=`echo $bf_cut_R*.fq.gz` # control_a.bam
        p1=`echo $f` # control_a_cut_R1.fq.gz
	p2=`echo $f | sed 's/\_R1/_R2/g'` # change control_a_cut_R1.fq.gz to control_a_cut_R2.fq.gz
        if ! [ -f $tt ] # -f check a file if exist; if output.bam don.t exist, do align; if script won't work,possible target bam already exist(like 0 size)
        then
            echo PE.p33 ... $tt $bf
            s_dt=$(date)
            echo $s_dt $bf' start......'
            # $bf for sample name and out put name
            # $p1 for pair1 $p2 for pair2
            # in do_bt2.pe.hg38.sh script use $1 $2 $2 to catch arguments
            bash ./trimmomatic.sh $bf $p1 $p2
            f_dt=$(date)
            echo $f_dt $bf' finished......'
            sleep 1
        fi
    else # SE
        if [[ $root != *R2.fq* ]]
        then
	    bf=`echo $root | sed -r 's/.fq//g' | sed -r s/.gz//g`
            tt=`echo $bf_cut_R*.fq.gz`
            if ! [ -f $tt ]
            then
                echo 'SE (no run) ...' $tt $bf
    	    fi
        fi
    fi
done   

# qc on trimmed data
fastqc -t 10 *fq.gz -o /Users/jplab/Desktop/cut_tag_data/ANNO_XS05KF2019110639_PM-XS05KF2019110639-10_BHFMK3CCX2_2021-05-30_01-50-52/qc/cut_qc

# 3. align(use mac flow, not seacer flow)
# samtools view -q 35
####################
bowtie2_align.sh
####################
opts='-p 20 -t --mm --very-sensitive --no-unal --no-mixed --no-discordant -I 10 -X 2000'
BOWTIE2_INDEXES='/Data/Ref_genome_anno/seqs_for_alignment_pipelines.ucsc_ids/hg38.bowtie_index'
echo 'start to align '$1 && bowtie2 $opts -x $BOWTIE2_INDEXES/hg38.bowtie_index -1 $2 -2 $3 | grep -E -v 'chrM|chrUn|random|RANDOM|KI|GL|chrEBV' | samtools view -@ 20 -q 35 -b -F 1804 -f 2 | samtools sort -@ 20 -n >$1.bam && echo $1' alignment done'


####################
bowtie_loop.sh
####################
# supports SE p33, p64, PE p33
for f in  ./*cut_R*.fq.gz
do
    root=`basename $f`

    if [[ $root == *_cut_R1.fq* ]]
    then # PE
	bf=`echo $root | sed -r 's/_cut_R1.fq//g' | sed -r s/.gz//g` # control_a
        tt=`echo $bf.sam` # control_a.bam
        p1=`echo $f` # control_a_cut_R1.fq.gz
	p2=`echo $f | sed 's/\_cut_R1/_cut_R2/g'` # change control_a_cut_R1.fq.gz to control_a_cut_R2.fq.gz
        if ! [ -f $tt ] # -f check a file if exist; if output.bam don.t exist, do align; if script won't work,possible target bam already exist(like 0 size)
        then
            echo PE.p33 ... $tt $bf
            s_dt=$(date)
            echo $s_dt $bf' start......'
            # $bf for sample name and out put name
            # $p1 for pair1 $p2 for pair2
            # in do_bt2.pe.hg38.sh script use $1 $2 $2 to catch arguments
            bash ./bowtie2_align.sh $bf $p1 $p2
            f_dt=$(date)
            echo $f_dt $bf' finished......'
            sleep 1
        fi
    else # SE
        if [[ $root != *_cut_R2.fq* ]]
        then
	    bf=`echo $root | sed -r 's/.fq//g' | sed -r s/.gz//g`
            tt=`echo $bf.sam`
            if ! [ -f $tt ]
            then
                echo 'SE (no run) ...' $tt $bf
                #qsub -N bt2se.$bf -v inp=$f,out=$bf do_bt2.se.sh

    	    fi
        fi
    fi
done

# 4.1 sort each bam and index
####################
sort_bam.sh
####################
for f in *.bam
do	
	bf=`basename $f`
	cf=`echo $bf | sed 's/.bam//g'`	
	echo $f
	echo 'command :' samtools sort -@ 20 $bf -o $cf.sortpos.bam 
    time samtools sort -@ 20 $bf -o $cf.sortpos.bam &> $cf.sortpos.log &
    sleep 2
done 

####################
index_bam.sh
####################
for f in *sortpos.bam
do	
	bf=`basename $f`
	cf=`echo $bf | sed 's/.sortpos.bam//g'`	
	echo $f
    time samtools index $bf &> $cf-bam_index.log &
    sleep 2
done

# 4.2 delete unsorted bam
rm *a.bam *b.bam

# 5. deeptools
########################################################
# deeptools
# 1. multiBamSummary for npz
# 2. plotCorrelation and pca
# 3. plotHeatmap
# 4. plotCoverage
# 5. bamPEFragmentSize for insert size
########################################################

########################################################
# deeptools_multiBamSummary.sh
# 1. multiBamSummary for npz
########################################################
# use deeptools multiBamSummary check bam files
# output are: compressed numpy matrix; raw_counts; scaling factor;
# bams must be indexed first
multiBamSummary \
    bins \
    --bamfiles *sortpos.bam \
    --blackListFileName $blacklist_path \
    --numberOfProcessors max/2 \
    --minMappingQuality 10 \
    -o ATAC_B14_titrition_readcounts.npz \
    --outRawCounts ATAC_B14_titrition_readCounts.tsv \
    --scalingFactors ATAC_B14_titrition_scaling_factor.tsv

########################################################
# deeptools_plotCorrelationPca.sh
# 2. plotCorrelation and pca
########################################################
# correlation and pca
plotCorrelation \
    -in ATAC_B14_titrition_readcounts.npz \
    --corMethod pearson \
    --skipZeros \
    --removeOutliers \
    --plotTitle "B14 titrition Spearman Correlation of Read Counts" \
    --whatToPlot scatterplot -o scatterplot_PearsonCorr_bigwigScores.png \
    --outFileCorMatrix PearsonCorr_bigwigScores.tsv

plotCorrelation \
    -in ATAC_B14_titrition_readcounts.npz \
    --corMethod spearman --skipZeros \
    --removeOutliers \
    --plotTitle "B14 titrition Spearman Correlation of Read Counts" \
    --whatToPlot heatmap \
    --colorMap RdYlBu \
    --plotNumbers \
    -o heatmap_SpearmanCorr_readCounts.png \
    --outFileCorMatrix SpearmanCorr_readCounts.tsv

plotPCA \
    -in ATAC_B14_titrition_readcounts.npz \
    -o PCA_readCounts.png \
    -T "PCA of B14 titrition read counts"

# heatmap
plotHeatmap \
    -m ATAC_B14_titrition_readcounts.npz \
    -out ATAC_B14_titrition.png \
    --colorMap RdBu \
    --whatToShow 'heatmap and colorbar' \
    --zMin -3 --zMax 3 \
    --kmeans 4

########################################################
# deeptools
# 4. plotCoverage
########################################################
# plot coverage for all samples
plotCoverage \
    -b *sortpos.bam \
    --plotFile ATAC_B14_titrition_coverage_all \
    --smartLabels \
    -n 1000000 \
    --plotTitle "ATAC_B14_titrition_coverage" \
    --outRawCounts ATAC_B14_titrition_coverage_all.tsv \
    --ignoreDuplicates \
    --minMappingQuality 10 \
    --blackListFileName $blacklist_path\
    -p 40

########################################################
# deeptools
# 5. bamPEFragmentSize for insert size
########################################################
# bampe_frag_size for all samples
bamPEFragmentSize \
    -hist ATAC_B14_titrition_fragmentSize_all.png \
    -T "Fragment size of PE all ATAC B14.titrition data" \
    -b ./*sortpos.bam -p 20 -bl $blacklist_path \
    --table ATAC_B14_titrition_frag_length.tsv

# bampe_frag_size and for each sample
blacklist_path=/Data/Ref_genome_anno/ENCODE_BLACK_LIST_REGION/hg38-blacklist.v2.bed
for f in *sortpos.bam
do
    echo $f
    bn=`basename $f`
    sn=`echo $bn | sed 's/.sortpos.bam//g'`
  
    bamPEFragmentSize \
        -hist $sn'_ATAC_titrition_fragmentSize.png' \
        -T "Fragment size of PE "$sn" ATAC data" \
        -b $bn \
        -p 20 \
        -bl $blacklist_path \
        --table $sn'_ATAC_titrition_frag_length.tsv'
done

# 6. merge reps bam(mainly for bigwig)
# FOXA1_B14_1.sortpos.bam FOXA1_B14_2.sortpos.bam
for f in  ./*sortpos.bam
do
    bn=`basename $f`

    if [[ $bn == *_a.sortpos.bam ]]
    then # PE
	bf=`echo $bn | sed -r 's/.sortpos.bam//g'` # FOXA1_B14_1
    gn=`echo $bn | sed -r 's/_a.sortpos.bam//g' | sed -r 's/_b.sortpos.bam//g'` # FOXA1_B14
    tt=`echo $gn.merged.bam` # control_a.bam
    p1=`echo $f` # FOXA1_B14_1.sortpos.bam
	p2=`echo $f | sed 's/\_a.sortpos/_b.sortpos/g'` # change control_a_cut_R1.fq.gz to control_a_cut_R2.fq.gz
        if ! [ -f $tt ] # -f check a file if exist; if output.bam don.t exist, do align; if script won't work,possible target bam already exist(like 0 size)
        then
            # $bf for sample name and out put name
            # $p1 for pair1 $p2 for pair2
            # in do_bt2.pe.hg38.sh script use $1 $2 $2 to catch arguments
            samtools merge -@ 20 $tt $p1 $p2
            # echo $tt $p1 $p2
            echo $bf' finished......'
            sleep 1
        fi
    else # SE
        echo 'data is not pair-end mod, please check'
    fi
done

# sort and index merged bam
for f in ./*merged.bam
do
    echo $f
    bn=`basename $f`
    sn=`echo $bn | sed 's/.merged.bam//g'`
    # after merge rep bam
    # bin size
    samtools sort -@ 20 $bn -o $sn'.merged.sort.pos.bam' && samtools index -@ 20 $sn'.merged.sort.pos.bam' && echo $sn ' sort and index done'
done

# delete merged bam, keep merged.sort.pos.bam
rm *merged.bam

# 7. convert merged bam to bigwig
for f in ./*merged.sort.pos.bam
do
    echo $f
    bn=`basename $f`
    sn=`echo $bn | sed 's/.merged.sort.pos.bam//g'`
    # after merge rep bam
    # bin size
    bamCoverage \
        --bam $bn \
        --outFileName $sn'.bw' \
        --binSize 1 \
        --normalizeUsing RPGC \
        --effectiveGenomeSize 2913022398 \
        --blackListFileName /Data/Ref_genome_anno/ENCODE_BLACK_LIST_REGION/hg38-blacklist.v2.bed \
        -p 40 &> $sn'bamcoverage.log'
done  

# 8. sort bam(already sort by pos) by name for bam to bed
# name sorted file can not be indexed
for f in *sortpos.bam
do
    bf=`basename $f`
    cf=`echo $bf | sed 's/.sortpos.bam//g'`
    echo $f
    echo 'command :' samtools sort -@ 20 -n $bf -o $cf.sortname.bam
    time samtools sort -@ 20 -n $bf -o $cf.sortname.bam &
    sleep 2
done


# 8. bam to bed 
# !!! important: bam file used by bedtools have to sorted by name!!! not pos !!!
# or will throw millions waring messages and can not got qcbc results
for f in ./*.sortpos.bam
do
    bn=`basename $f` # bn B14_2_a.bam
    rn=`echo $bn | sed s/.sortpos.bam/.bed.gz/g` # rn B14_2_a.bed.gz
    if ! [ -f $rn ] # already existed bed file will not be  handled
    then
	echo $rn
    # f B14_2_a.bam
    # rn B14_2_a.bed.gz
    bash bam2bed.sh $f $rn &
	sleep 2
    fi
done

####################
# bam2bed.sh
####################
# -q 35
# THe PE with QC variant
# bam_pe to bed
# bedpe : Write BAM alignments in BEDPE format. Only one alignment from paired-end reads will be reported.
# mate1 : When writing BEDPE (-bedpe) format, always report mate one as the first BEDPE “block”.
samtools view -@ 20 -q 35 -b $1 | bedtools bamtobed -mate1 -bedpe > $2.temp 

cat $2.temp | awk 'BEGIN{OFS="\t"}{print $1,$2,$4,$6,$9,$10}' |sort | uniq -c | awk 'BEGIN{mt=0;m0=0;m1=0;m2=0} ($1==1){m1=m1+1} ($1==2){m2=m2+1} {m0=m0+1} {mt=mt+$1} END{printf "TotalReadPairs=%d\nDistinctReadPairs=%d\nOneReadPair=%d\nTwoReadPairs=%d\nNRF=%f\nPBC1=%f\nPBC2=%f\n",mt,m0,m1,m2,m0/mt,m1/m0,m1/m2}' > $2.qcbc

echo '# QCBC adivce' >>$2.qcbc
echo '# PBC1 is the primary measure. Provisionally,' >>$2.qcbc 
echo '# 0-0.5 is severe bottlenecking' >>$2.qcbc
echo '# 0.5-0.8 is moderate bottlenecking' >>$2.qcbc
echo '# 0.8-0.9 is mild bottlenecking' >>$2.qcbc
echo '# 0.9-1.0 is no bottlenecking' >>$2.qcbc

# Actual kept BED
# Remove the qula scores so you can make is unqiue
# '!x[$0]++': This one-liner removes duplicate lines from text input without pre-sorting.
cat $2.temp | awk '{FS=OFS="\t"} { if ($9 =="+") {$2=$2; print $1,$2,$6,".",0,$9}  else if ($9 == "-") {$3=$3; print $1,$5,$3,".",0,$9}}' | awk '!x[$0]++' - | gzip >$2
rm $2.temp

# 9. remove black list region from bed
for f in ./*bed.gz
do
    echo $f
    bn=`basename $f`
    sn=`echo $bn | sed 's/.bed.gz//g'`
  
    bedtools intersect -v -a $f -b $blacklist_path > $sn.blacklistremoved.bed
done

# remove origin bed file
rm *.bed.gz

# 10. bed to flat
import sys, glbase3; 
glbase3.bed_to_flat(
    ['B14_0_4_a.blacklistremoved.bed','B14_0_4_b.blacklistremoved.bed'], 
    'B14_0_4.flat', 
    name='B14_0_4', 
    isPE=True, gzip=False);

glbase3.bed_to_flat(
    ['B14_0_8_a.blacklistremoved.bed','B14_0_8_b.blacklistremoved.bed'], 
    'B14_0_8.flat', 
    name='B14_0_8', 
    isPE=True, gzip=False);    

glbase3.bed_to_flat(
    ['B14_1_2_a.blacklistremoved.bed','B14_1_2_b.blacklistremoved.bed'], 
    'B14_1_2.flat', 
    name='B14_1_2', 
    isPE=True, gzip=False);   

glbase3.bed_to_flat(
    ['B14_1_6_a.blacklistremoved.bed','B14_1_6_b.blacklistremoved.bed'], 
    'B14_1_6.flat', 
    name='B14_1_6', 
    isPE=True, gzip=False);

# 11. macs call peak
for f in ./*.blacklistremoved.bed
do
    echo $f
    bn=`basename $f`
    sn=`echo $bn | sed 's/.blacklistremoved.bed//g'`
  
    macs2 callpeak -q 0.05 -t $f -n $sn -g hs -f BEDPE --outdir ./peak &> ./peak/$sn.err
done

# 12. redefine peak
import os
from glbase3 import *

peaks = [
    genelist(filename="./B14_0_4_a_peaks.narrowPeak", format=format.bed, gzip=False),
    genelist(filename="./B14_0_4_b_peaks.narrowPeak", format=format.bed, gzip=False),
    genelist(filename="./B14_0_8_a_peaks.narrowPeak", format=format.bed, gzip=False),
    genelist(filename="./B14_0_8_b_peaks.narrowPeak", format=format.bed, gzip=False),
    # genelist(filename="./B14_1_2_a_peaks.narrowPeak", format=format.bed, gzip=False),
    genelist(filename="./B14_1_2_b_peaks.narrowPeak", format=format.bed, gzip=False),
    # genelist(filename="./B14_1_6_a_peaks.narrowPeak", format=format.bed, gzip=False),
    genelist(filename="./B14_1_6_b_peaks.narrowPeak", format=format.bed, gzip=False),
    ]

gl = glglob()
superset = gl.chip_seq_cluster(list_of_peaks=peaks)
superset.sort('loc')

# Output a redefined peaklist for each inputted flat file:
flats = [
    flat_track(filename='../B14_0_4.flat'),
    flat_track(filename='../B14_0_8.flat'),
    flat_track(filename='../B14_1_2.flat'),
    flat_track(filename='../B14_1_6.flat'),
    ]

# todo : Z_threshold ? how to decide?
rets = gl.redefine_peaks(superset, flats, filename='B14_titrition_z1.8thres_models', Z_threshold=1.8)

for f in rets:
    rets[f].saveBED('%s.bed' % f.replace(' ', '_'), uniqueID=True)
    rets[f].save('%s.glb' % f.replace(' ', '_'))


# try to group peak
import glob
import os
from glbase3 import * 

config.draw_mode = "pdf"

# merged flats
trks = [
    flat_track(filename='../FOXA1_B14.merge.flat'),
    flat_track(filename='../FOXA1_DMSO.merge.flat'),
    flat_track(filename='../FOXA1_TGF.merge.flat')
    ]
peaks = [    
    genelist(filename="FOXA1_B14.merge.bed", format=format.minimal_bed, gzip=False),
    genelist(filename="FOXA1_DMSO.merge.bed", format=format.minimal_bed, gzip=False),
    genelist(filename="FOXA1_TGF.merge.bed", format=format.minimal_bed, gzip=False),
    ]
gl = glglob()
# # two steps:
# # 1.chip_seq_cluster_heatmap return a genlist 
# 2. use gl to process this list to expression object
ret = gl.chip_seq_cluster_heatmap(peaks, trks, "heatmap.pdf",
    cache_data="data.bin",
    normalise=True,
    imshow=True,
    size=[10,27], # figSize
    pileup_distance=1000,
    bins=200, read_extend=0) 
gl.chip_seq_cluster_pileup(filename="clus/redefine_clusters.png")

for cid in ret:
    print("cid:", cid, "len:", len(ret[cid]["genelist"]))
    print("cid:", cid, "genlist:", ret[cid]["genelist"])
    ret[cid]["genelist"].saveBED(filename="clus/cid_%s.bed" % cid, uniqueID=True)


# motif    
findMotifsGenome.pl B14_0_4.bed hg38 ../motif/B14_0_4 -mask -S 15 -dumpFasta -p 20 &> B14_0_4.err &
findMotifsGenome.pl B14_0_8.bed hg38 ../motif/B14_0_8 -mask -S 15 -dumpFasta -p 20 &> B14_0_8.err &
findMotifsGenome.pl B14_1_2.bed hg38 ../motif/B14_1_2 -mask -S 15 -dumpFasta -p 20 &> B14_1_2.err &
findMotifsGenome.pl B14_1_6.bed hg38 ../motif/B14_1_6 -mask -S 15 -dumpFasta -p 20 &> B14_1_6.err &
