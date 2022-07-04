# if use  -size default(default 200), use summit.bed
# if use  -size given, use narrowpeak
# -norevopp and both strands means find motif only on pos strand and both strands

# only pos strand and default size
findMotifsGenome.pl /Users/jplab/Desktop/snakepipes_Data/ATAC_seq_shFOXA1_snakepipes_results/MACS_allLen_noDUP_noChrM/F-B1/F-B1_summits.bed hg38 F-B1.posStrand.sizedefault -norevopp -mask -p 20 -S 15 &> F-B1.posStrand.sizedefault.log

# only pos strand and given size
findMotifsGenome.pl /Users/jplab/Desktop/snakepipes_Data/ATAC_seq_shFOXA1_snakepipes_results/MACS_allLen_noDUP_noChrM/F-B1/F-B1_peaks.narrowPeak hg38 F-B1.posStrand.sizeGiven -size given -norevopp -mask -p 20 -S 15 &> F-B1.posStrand.sizeGiven.log

# only pos strand and default size and use narrowpeak
findMotifsGenome.pl /Data/Projects/ATAC_seq_shFOXA1_20220302_snakepipes/MACS_allLen_noDUP_noChrM/F-B1/F-B1_peaks.narrowPeak hg38 F-B1.posStrand.narrowPeak.sizedefault -norevopp -mask -p 20 -S 15 &> F-B1.posStrand.narrowPeak.sizedefault.log

# both strands and default size
findMotifsGenome.pl /Users/jplab/Desktop/snakepipes_Data/ATAC_seq_shFOXA1_snakepipes_results/MACS_allLen_noDUP_noChrM/F-B1/F-B1_summits.bed hg38 F-B1.bothStrands.sizedefault  -mask -p 20 -S 15 &> F-B1.bothStrands.sizedefault.log

# both strands and given size
findMotifsGenome.pl /Users/jplab/Desktop/snakepipes_Data/ATAC_seq_shFOXA1_snakepipes_results/MACS_allLen_noDUP_noChrM/F-B1/F-B1_peaks.narrowPeak hg38 F-B1.bothStrands.sizeGiven -size given  -mask -p 20 -S 15 &> F-B1.bothStrands.sizeGiven.log
