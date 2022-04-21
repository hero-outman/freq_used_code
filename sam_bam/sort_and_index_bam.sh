# can not run together
for f in *.noDup.nochrM.bam
do	
	bf=`basename $f`
	cf=`echo $bf | sed 's/.noDup.nochrM.bam//g'`	
	echo $f
	echo 'command :' samtools sort -@ 20 $bf -o $cf.noDup.nochrM.sortpos.bam 
    time samtools sort -@ 20 $bf -o $cf.noDup.nochrM.sortpos.bam &
    sleep 2
done
echo 'All bam files are sorted'


for f in *.noDup.nochrM.sortpos.bam
do	
	bf=`basename $f`
	cf=`echo $bf | sed 's/.noDup.nochrM.sortpos.bam//g'`	
	echo $f
    time samtools index $bf &
    sleep 2
done
echo 'All bam files are indexed'

echo 'delete unsorted bam'
rm *a.bam *b.bam