# note: command is alignmentSieve, not alignmentSieve.py !!!
alignmentSieve \
    -v -p 20 \
    -b B14_2_a.filtered.bam -o ./great_than_150bp/B14_2_a.filtered_gt150bp.bam \
    --minMappingQuality 20 \
    --minFragmentLength 150 \
    --maxFragmentLength 0 \
    --filterMetrics B14_2_a_filtered_log.txt &