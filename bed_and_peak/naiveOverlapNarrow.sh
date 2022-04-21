#!/bin/bash

##########################

# naiveOverlapNarrow
# November 2019
# Jake Reske
# Michigan State University
# reskejak@msu.edu
# https://github.com/reskejak

# Computing ENCODE-defined naive overlapping narrowPeak set from two ATAC or ChIP individual replicate narrowPeak and replicate-pooled narrowPeak sets
# Definition from Kundaje et al.: "Find pooled peaks that overlap Rep1 and Rep2 where overlap is defined as the fractional overlap wrt any one of the overlapping peak pairs >= 0.5".
# usage: naiveOverlapNarrow treat1_peaks.narrowPeak treat2_peaks.narrowPeak treat_pool_peaks.narrowPeak > treat_overlap_peaks.narrowPeak
# ensure to firstly make executable by: chmod u+x naiveOverlapNarrow

# dependencies: bedtools (for intersectBed / bedtools intersect)

##########################

REP1=$1
REP2=$2
POOL=$3

# narrowPeak
intersectBed -wo \
-a ${POOL} -b ${REP1} | \
awk 'BEGIN{FS="\t";OFS="\t"}{s1=$3-$2; s2=$13-$12; if (($21/s1 >= 0.5) || ($21/s2 >= 0.5)) {print $0}}' | \
cut -f 1-10 | sort | uniq | \
intersectBed -wo \
-a stdin -b ${REP2} | \
awk 'BEGIN{FS="\t";OFS="\t"}{s1=$3-$2; s2=$13-$12; if (($21/s1 >= 0.5) || ($21/s2 >= 0.5)) {print $0}}' | \
cut -f 1-10 | sort | uniq

# intersectBed -wo will combine two peak in one line and with overlaped bp number.
# $21 is overlap bp number.
# s1=$3-$2; s2=$13-$12 are the length of two peaks intersected by bedtools

# intersectBed -wo -a /Users/jplab/Desktop/DAILY_CODE_DATA/2022-4/code/4-8_test_macs_pooled_peak/B14_2_pool.filtered.demo.narrowPeak -b /Users/jplab/Desktop/DAILY_CODE_DATA/2022-4/code/4-8_test_macs_pooled_peak/FOXA1_B14_1_demo.narrowPeak | \
# awk 'BEGIN{FS="\t";OFS="\t"}{s1=$3-$2; s2=$13-$12; if (($21/s1 >= 0.5) || ($21/s2 >= 0.5)) {print $0}}' | cut -f 1-10

# output：
# chr1	1079569	1080350	B14_2_pool.filtered_peak_17	286	.	8.41655	32.9733	28.6788	580	chr1	1079580	1080302	FOXA1_B14_1_peak_8	163	.	8.19621	20.7926	16.3178	135	722