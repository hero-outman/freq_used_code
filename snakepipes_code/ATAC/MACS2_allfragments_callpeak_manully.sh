for f in /Data/Projects/CUT_TAG_snakepipe/DNA_mapping/filtered_bam/*.bam
do
    bn=`basename $f`
    rn=`echo $bn | sed s/.bam//g`

	echo $rn
    # f B14_2_a.bam
    # rn B14_2_a.bed.gz
    macs2 callpeak -q 0.01 -t $f -n $rn -g hs -f BAMPE --outdir /Data/Projects/CUT_TAG_snakepipe/DNA_mapping/MACS2_manual &> $rn.err
	sleep 2
    
done