suffix='.filtered.bam'
# exclude_chr='chrM\|KI*\|GL*'
exclude_chr='chrM'

for f in *$suffix
do	
	bf=`basename $f`
	sample_name=`echo $bf | sed 's/'$suffix'//g'`	
	echo $f
    samtools idxstats $f | cut -f 1 | grep -v $exclude_chr | xargs samtools view -q 30 -b $f > $sample_name'.noDup.nochrM.bam'
    sleep 2
done