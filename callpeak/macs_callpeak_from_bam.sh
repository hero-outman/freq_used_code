fileList=("B14_2_a.filtered.bam" "B14_2_b.filtered.bam" "control_a.filtered.bam" "control_b.filtered.bam" "TGF_a.filtered.bam" "TGF_b.filtered.bam")

for f in $fileList
do
    bn=`basename $f`
    rn=`echo $bn | sed s/.bam//g`

	echo $rn
    # f B14_2_a.bam
    # rn B14_2_a
    macs2 callpeak -q 0.01 -t '/Data/Projects/ATAC_snakepipe/DNA_mapping/filtered_bam/'$f -n $rn -g hs -f BAMPE --outdir /Data/Projects/ATAC_snakepipe/DNA_mapping/MACS2/allFraglength/ &> $rn.err
	sleep 1
    
done





for f in /Data/Projects/ATAC_snakepipe/DNA_mapping/filtered_bam/great_than_150bp/*.bam
do
    bn=`basename $f`
    rn=`echo $bn | sed s/.bam//g`

	echo $rn
    # f B14_2_a.bam
    # rn B14_2_a
    macs2 callpeak -q 0.01 -t $f -n $rn -g hs -f BAMPE --outdir /Data/Projects/ATAC_snakepipe/DNA_mapping/MACS2/gt150bp_fraglength/ &> $rn.err
	sleep 2
    
done



suffix='.noDup.nochrM.bam'
for f in *$suffix
do	
	bf=`basename $f`
	sample_name=`echo $bf | sed 's/'$suffix'//g'`	
	echo $f
    
    macs2 callpeak \
        -q 0.01 \
        -t $f \
        -f BAMPE \
        -g hs \
        -n $sample_name \
        --outdir $sample_name &> $sample_name'-macs2.log' &

    sleep 2
done

# pool
macs2 callpeak -q 0.01 \
    -t /Users/jplab/Desktop/snakepipes_Data/CUT_TAG_snakepipes/FOXA1_H3K27ac/noDUP_nochrM_bam/FOXA1_B14_1.noDup.nochrM.sortpos.bam /Users/jplab/Desktop/snakepipes_Data/CUT_TAG_snakepipes/FOXA1_H3K27ac/noDUP_nochrM_bam/FOXA1_B14_2.noDup.nochrM.sortpos.bam \
    -f BAMPE \
    -g hs \
    -n FOXA1_B14_noDup.nochrM.pool_peaks.narrowPeak \
    --outdir FOXA1_B14_noDup.nochrM.pool_peaks &> FOXA1_B14_pool_peaks-macs2.log & 