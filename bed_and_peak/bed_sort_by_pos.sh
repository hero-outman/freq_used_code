for f in /Users/jplab/Desktop/snakepipes_Data/ATAC_seq_shFOXA1_snakepipes_results/MACS_allLen_noDUP_noChrM/naiveOverlap_q001/*.narrowPeak
do	
	bf=`basename $f`
	cf=`echo $bf | sed 's/.narrowPeak//g'`	
	echo $f
	echo 'command :' bedtools sort -i $f > ${cf}.sorted.narrowPeak
    time bedtools sort -i $f > ${cf}.sorted.narrowPeak &
    sleep 2
done
echo 'All bed files are sorted'