#!/bin/sh

# get 2 bed files's common open region(NOT the final common open region)
# like navieoverlap, if 2 regions's overlap > 0.5, keep the region.
# !!! NOT the final common open region, becasue these region include both open but one condition is more enrich than other condition, these region will be removed in next step !

##### part1: get the primary common open region
if [ "$#" -eq 4 ]
    then
            # if [ -d $var1 ]
            # then
            # echo directory ${var1} exist
            # else
            # echo Directory ${var1} Does not exists
            # fi
            # if [ -d $var2 ]
            # then
            # echo directory ${var2} exist
            # else
            # echo Directory ${var2} Does not exists
            # fi
        echo 'arguments are ok'
    else
    echo ''
    echo "Please give 4 arguments: 1.bed file 2.bed file 3.conditionA_enriched_bed file 4.conditionB_enriched_bed file"
    echo ''
    exit 1
    fi


BED1=$1
BED2=$2
enriched_region_in_condition1=$3
enriched_region_in_condition2=$4
file_name_1=`basename $BED1`
file_name_2=`basename $BED2`
sample_name1=`echo $file_name_1 | sed 's/_overlap_peaks.narrowPeak//g'`
sample_name2=`echo $file_name_2 | sed 's/_overlap_peaks.narrowPeak//g'`

# narrowPeak
intersectBed -wo \
-a ${BED1} -b ${BED2} | \
awk 'BEGIN{FS="\t";OFS="\t"}{s1=$3-$2; s2=$13-$12; if (($21/s1 >= 0.5) || ($21/s2 >= 0.5)) {print $0}}' | \
cut -f 1-10 | sort | uniq > ${sample_name1}_${sample_name2}_overlap.narrowPeak.temp &&
# sleep 10
echo ${sample_name1}_${sample_name2}_overlap.narrowPeak.temp 'has regions: ' | wc -l ${sample_name1}_${sample_name2}_overlap.narrowPeak.temp &&

intersectBed -wo \
-a ${BED2} -b ${BED1} | \
awk 'BEGIN{FS="\t";OFS="\t"}{s1=$3-$2; s2=$13-$12; if (($21/s1 >= 0.5) || ($21/s2 >= 0.5)) {print $0}}' | \
cut -f 1-10 | sort | uniq > ${sample_name2}_${sample_name1}_overlap.narrowPeak.temp &&
# sleep 10
echo ${sample_name2}_${sample_name1}_overlap.narrowPeak.temp 'has regions: ' | wc -l ${sample_name2}_${sample_name1}_overlap.narrowPeak.temp && 

# merge book-end regions
cat ${sample_name1}_${sample_name2}_overlap.narrowPeak.temp \
${sample_name2}_${sample_name1}_overlap.narrowPeak.temp \
| sort -k1,1 -k2,2n | mergeBed -i stdin -o collapse -c 4 > ${sample_name1}_${sample_name2}_primary.common.open.bed &&
# sleep 10
echo ${sample_name1} ${sample_name2} 'has common open regions: ' | wc -l ${sample_name1}_${sample_name2}_primary.common.open.bed &&

rm ${sample_name1}_${sample_name2}_overlap.narrowPeak.temp ${sample_name2}_${sample_name1}_overlap.narrowPeak.temp &&

echo 'Primary Common Open region Done!!! Change the file name as needed !!!'
echo 'This is NOT the final common open region !!! '
echo 'Becasue these region include both open but one condition is more enrich than other condition, these region will be removed in next step'

##### part2: get the final common open region
# remove both open but in conditionA more enriched region
bedtools intersect \
-a ${sample_name1}_${sample_name2}_primary.common.open.bed \
-b ${enriched_region_in_condition1} -v \
> ${sample_name1}_${sample_name2}_primary.common.open.notEnrichedInCondition1.bed.temp &&

# remove both open but in conditionB more enriched region, this is common open region
bedtools intersect \
-a ${sample_name1}_${sample_name2}_primary.common.open.notEnrichedInCondition1.bed.temp \
-b ${enriched_region_in_condition2} -v \
> ${sample_name1}_${sample_name2}_primary.common.open.notEnrichedInCondition1.notEnrichedInCondition2.bed.temp &&

# rename file, remove temp file
mv ${sample_name1}_${sample_name2}_primary.common.open.notEnrichedInCondition1.notEnrichedInCondition2.bed.temp ${sample_name1}_${sample_name2}_final.common.open.bed &&

rm ${sample_name1}_${sample_name2}_primary.common.open.notEnrichedInCondition1.bed.temp

echo ${sample_name1} ${sample_name2} 'has final common open regions: ' | wc -l ${sample_name1}_${sample_name2}_final.common.open.bed &&
echo 'Done!!! Final common open region file done!'

# need to check sh and nc more enrichment  region by intersect with final common open region file 
# what if sh or nc more enriched region intersect primary common open region file? just open in sh or nc ? need to test.