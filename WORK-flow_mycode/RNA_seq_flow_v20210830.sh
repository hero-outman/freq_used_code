##############################################################################################################################

# RNA_seq_flow_v20210830

##############################################################################################################################

blacklist_path=/Data/Ref_genome_anno/ENCODE_BLACK_LIST_REGION/hg38-blacklist.v2.bed
hg38_path=/Data/Ref_genome_anno/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna
hg38_chrom_sizes=/Data/Ref_genome_anno/hg38.chrom.sizes


# 0.build link
ln -s /Data/RAW_DATA/RNA-seq_spheroid_2021-8/Rawdata/*/*fq.gz /Data/Projects/RNA-seq_spheroid_2021-8/fastq


# 1. qc on rawdata
fastqc -t 20  ./*fq.gz -o /Users/jplab/Desktop/cut_tag_data/ANNO_XS05KF2019110639_PM-XS05KF2019110639-10_BHFMK3CCX2_2021-05-30_01-50-52/qc/raw_qc

# 2. cut adaptor and low quality data
####################
trimmomatic.sh
####################
echo 'begin to cut adaptor of ' $1 && \
java -jar /Data/Tools/Trimmomatic-0.39/trimmomatic-0.39.jar \
PE -phred33 \
-threads 10 \
$2 $3 $1'_cut_R1.fq.gz' $1'_cut.unpaired_R1.fq.gz' $1'_cut_R2.fq.gz' $1'_cut.unpaired_R2.fq.gz' \
ILLUMINACLIP:/Data/Tools/Trimmomatic-0.39/adapters/TruSeq3-PE.fa:2:30:10:2:true LEADING:3 SLIDINGWINDOW:4:15 MINLEN:30 \
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
hisat2_align_v20210830.sh
####################
# -p : threads to use
# -t : show time costing
# -x : ref genome
opts='-p 8 -t --mm --no-unal --no-mixed --no-discordant'
hisat2_index='/Data/Ref_genome_anno/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.hisat2_index/'
# keep chrM?
remove_sequence='chrUn|random|RANDOM|KI|GL|chrEBV'

echo 'start to align: ' $1 ' at ' $(date) && hisat2 $opts -x $hisat2_index -1 $2 -2 $3 | grep -E -v $remove_sequence | samtools view -@ 20 -q 25 -b -F 1804 -f 2 | samtools sort -@ 20 -n >$1.bam && echo $1' alignment done, end at: '  $(date)


####################
hisat2_align_loop.sh
####################
# supports SE p33, p64, PE p33
for f in  ./*cut_R*.fq.gz
do
    root=`basename $f`

    if [[ $root == *_cut_R1.fq* ]]
    then # PE
	bf=`echo $root | sed -r 's/_cut_R1.fq//g' | sed -r s/.gz//g` # control_a
        tt=`echo $bf.bam` # control_a.bam
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
            bash ./hisat2_align_v20210830.sh $bf $p1 $p2
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
    --numberOfProcessors max/2 \
    --minMappingQuality 10 \
    -o ../deeptools/2umol_spheriod_readcounts.npz \
    --outRawCounts ../deeptools/2umol_spheriod_readcounts.tsv \
    --scalingFactors ../deeptools/2umol_spheriod_readcounts.tsv

########################################################
# deeptools_plotCorrelationPca.sh
# 2. plotCorrelation and pca
########################################################
# correlation and pca
plotCorrelation \
    -in 2umol_spheriod_readcounts.npz \
    --corMethod pearson \
    --skipZeros \
    --removeOutliers \
    --plotTitle "2umol_spheriod RNA seq Spearman Correlation of Read Counts" \
    --whatToPlot scatterplot -o 2umol_spheriod_scatterplot_PearsonCorr_bigwigScores.png \
    --outFileCorMatrix 2umol_spheriod_PearsonCorr_bigwigScores.tsv

plotCorrelation \
    -in 2umol_spheriod_readcounts.npz \
    --corMethod spearman --skipZeros \
    --removeOutliers \
    --plotTitle "2umol_spheriod RNA seq Spearman Correlation of Read Counts" \
    --whatToPlot heatmap \
    --colorMap RdYlBu \
    --plotNumbers \
    -o heatmap_SpearmanCorr_readCounts.png \
    --outFileCorMatrix SpearmanCorr_readCounts.tsv

plotPCA \
    -in 2umol_spheriod_readcounts.npz \
    -o PCA_readCounts.png \
    -T "PCA of 2umol_spheriod RNAread counts"

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
    --plotFile 2umol_RNA_spheroid_coverage \
    --smartLabels \
    -n 1000000 \
    --plotTitle "2umol_RNA_spheroid_coverage" \
    --outRawCounts 2umol_RNA_spheroid_coverage.tsv \
    --ignoreDuplicates \
    --minMappingQuality 10 \
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
for f in *sortpos.bam
do
    echo $f
    bn=`basename $f`
    sn=`echo $bn | sed 's/.sortpos.bam//g'`
  
    bamPEFragmentSize \
        -hist $sn'_RNA_spheroid_fragmentSize.png' \
        -T "Fragment size of PE "$sn" RNA_spheroid data" \
        -b $bn \
        -p 20 \
        --table $sn'_RNA_spheroid_frag_length.tsv'
done

########################################################
# featurecounts.sh
########################################################
# count reads
opts='-T 8 -t exon -g gene_id -p -B -C'
anno_path=/Data/Ref_genome_anno/formal_bak/Homo_sapiens.GRCh38.103.gtf.gz
featureCounts $opts -a $anno_path -o 2umol.RNA.spheroid.readcounts.tsv *sortpos.bam &> counts.log
