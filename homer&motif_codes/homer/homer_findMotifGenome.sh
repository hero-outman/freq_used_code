# if use  -size default(default 200), use summit.bed
# if use  -size given, use narrowpeak
# -norevopp and both strands means find motif only on pos strand and both strands

# only pos strand and default size
findMotifsGenome.pl /Users/jplab/Desktop/snakepipes_Data/ATAC_seq_shFOXA1_snakepipes_results/MACS_allLen_noDUP_noChrM/F-B1/F-B1_summits.bed hg38 F-B1.posStrand.sizedefault -norevopp -mask -p 20 -S 15 &> F-B1.posStrand.sizedefault.log
findMotifsGenome.pl /Users/jplab/Desktop/snakepipes_Data/ATAC_seq_shFOXA1_snakepipes_results/MACS_allLen_noDUP_noChrM/F-B2/F-B2_summits.bed hg38 F-B2.posStrand.sizedefault -norevopp -mask -p 20 -S 15 &> F-B2.posStrand.sizedefault.log
findMotifsGenome.pl /Users/jplab/Desktop/snakepipes_Data/ATAC_seq_shFOXA1_snakepipes_results/MACS_allLen_noDUP_noChrM/F-D1/F-D1_summits.bed hg38 F-D1.posStrand.sizedefault -norevopp -mask -p 20 -S 15 &> F-D1.posStrand.sizedefault.log
findMotifsGenome.pl /Users/jplab/Desktop/snakepipes_Data/ATAC_seq_shFOXA1_snakepipes_results/MACS_allLen_noDUP_noChrM/F-D2/F-D2_summits.bed hg38 F-D2.posStrand.sizedefault -norevopp -mask -p 20 -S 15 &> F-D2.posStrand.sizedefault.log
findMotifsGenome.pl /Users/jplab/Desktop/snakepipes_Data/ATAC_seq_shFOXA1_snakepipes_results/MACS_allLen_noDUP_noChrM/F-T1/F-T1_summits.bed hg38 F-T1.posStrand.sizedefault -norevopp -mask -p 20 -S 15 &> F-T1.posStrand.sizedefault.log
findMotifsGenome.pl /Users/jplab/Desktop/snakepipes_Data/ATAC_seq_shFOXA1_snakepipes_results/MACS_allLen_noDUP_noChrM/F-T2/F-T2_summits.bed hg38 F-T2.posStrand.sizedefault -norevopp -mask -p 20 -S 15 &> F-T2.posStrand.sizedefault.log
findMotifsGenome.pl /Users/jplab/Desktop/snakepipes_Data/ATAC_seq_shFOXA1_snakepipes_results/MACS_allLen_noDUP_noChrM/NC-B1/NC-B1_summits.bed hg38 NC-B1.posStrand.sizedefault -norevopp -mask -p 20 -S 15 &> NC-B1.posStrand.sizedefault.log
findMotifsGenome.pl /Users/jplab/Desktop/snakepipes_Data/ATAC_seq_shFOXA1_snakepipes_results/MACS_allLen_noDUP_noChrM/NC-B2/NC-B2_summits.bed hg38 NC-B2.posStrand.sizedefault -norevopp -mask -p 20 -S 15 &> NC-B2.posStrand.sizedefault.log
findMotifsGenome.pl /Users/jplab/Desktop/snakepipes_Data/ATAC_seq_shFOXA1_snakepipes_results/MACS_allLen_noDUP_noChrM/NC-D1/NC-D1_summits.bed hg38 NC-D1.posStrand.sizedefault -norevopp -mask -p 20 -S 15 &> NC-D1.posStrand.sizedefault.log
findMotifsGenome.pl /Users/jplab/Desktop/snakepipes_Data/ATAC_seq_shFOXA1_snakepipes_results/MACS_allLen_noDUP_noChrM/NC-D2/NC-D2_summits.bed hg38 NC-D2.posStrand.sizedefault -norevopp -mask -p 20 -S 15 &> NC-D2.posStrand.sizedefault.log
findMotifsGenome.pl /Users/jplab/Desktop/snakepipes_Data/ATAC_seq_shFOXA1_snakepipes_results/MACS_allLen_noDUP_noChrM/NC-T1/NC-T1_summits.bed hg38 NC-T1.posStrand.sizedefault -norevopp -mask -p 20 -S 15 &> NC-T1.posStrand.sizedefault.log
findMotifsGenome.pl /Users/jplab/Desktop/snakepipes_Data/ATAC_seq_shFOXA1_snakepipes_results/MACS_allLen_noDUP_noChrM/NC-T2/NC-T2_summits.bed hg38 NC-T2.posStrand.sizedefault -norevopp -mask -p 20 -S 15 &> NC-T2.posStrand.sizedefault.log

# only pos strand and given size
findMotifsGenome.pl /Users/jplab/Desktop/snakepipes_Data/ATAC_seq_shFOXA1_snakepipes_results/MACS_allLen_noDUP_noChrM/F-B1/F-B1_peaks.narrowPeak hg38 F-B1.posStrand.sizeGiven -size given -norevopp -mask -p 20 -S 15 &> F-B1.posStrand.sizeGiven.log
findMotifsGenome.pl /Users/jplab/Desktop/snakepipes_Data/ATAC_seq_shFOXA1_snakepipes_results/MACS_allLen_noDUP_noChrM/F-B2/F-B2_peaks.narrowPeak hg38 F-B2.posStrand.sizeGiven -size given -norevopp -mask -p 20 -S 15 &> F-B2.posStrand.sizeGiven.log
findMotifsGenome.pl /Users/jplab/Desktop/snakepipes_Data/ATAC_seq_shFOXA1_snakepipes_results/MACS_allLen_noDUP_noChrM/F-D1/F-D1_peaks.narrowPeak hg38 F-D1.posStrand.sizeGiven -size given -norevopp -mask -p 20 -S 15 &> F-D1.posStrand.sizeGiven.log
findMotifsGenome.pl /Users/jplab/Desktop/snakepipes_Data/ATAC_seq_shFOXA1_snakepipes_results/MACS_allLen_noDUP_noChrM/F-D2/F-D2_peaks.narrowPeak hg38 F-D2.posStrand.sizeGiven -size given -norevopp -mask -p 20 -S 15 &> F-D2.posStrand.sizeGiven.log
findMotifsGenome.pl /Users/jplab/Desktop/snakepipes_Data/ATAC_seq_shFOXA1_snakepipes_results/MACS_allLen_noDUP_noChrM/F-T1/F-T1_peaks.narrowPeak hg38 F-T1.posStrand.sizeGiven -size given -norevopp -mask -p 20 -S 15 &> F-T1.posStrand.sizeGiven.log
findMotifsGenome.pl /Users/jplab/Desktop/snakepipes_Data/ATAC_seq_shFOXA1_snakepipes_results/MACS_allLen_noDUP_noChrM/F-T2/F-T2_peaks.narrowPeak hg38 F-T2.posStrand.sizeGiven -size given -norevopp -mask -p 20 -S 15 &> F-T2.posStrand.sizeGiven.log
findMotifsGenome.pl /Users/jplab/Desktop/snakepipes_Data/ATAC_seq_shFOXA1_snakepipes_results/MACS_allLen_noDUP_noChrM/NC-B1/NC-B1_peaks.narrowPeak hg38 NC-B1.posStrand.sizeGiven -size given -norevopp -mask -p 20 -S 15 &> NC-B1.posStrand.sizeGiven.log
findMotifsGenome.pl /Users/jplab/Desktop/snakepipes_Data/ATAC_seq_shFOXA1_snakepipes_results/MACS_allLen_noDUP_noChrM/NC-B2/NC-B2_peaks.narrowPeak hg38 NC-B2.posStrand.sizeGiven -size given -norevopp -mask -p 20 -S 15 &> NC-B2.posStrand.sizeGiven.log
findMotifsGenome.pl /Users/jplab/Desktop/snakepipes_Data/ATAC_seq_shFOXA1_snakepipes_results/MACS_allLen_noDUP_noChrM/NC-D1/NC-D1_peaks.narrowPeak hg38 NC-D1.posStrand.sizeGiven -size given -norevopp -mask -p 20 -S 15 &> NC-D1.posStrand.sizeGiven.log
findMotifsGenome.pl /Users/jplab/Desktop/snakepipes_Data/ATAC_seq_shFOXA1_snakepipes_results/MACS_allLen_noDUP_noChrM/NC-D2/NC-D2_peaks.narrowPeak hg38 NC-D2.posStrand.sizeGiven -size given -norevopp -mask -p 20 -S 15 &> NC-D2.posStrand.sizeGiven.log
findMotifsGenome.pl /Users/jplab/Desktop/snakepipes_Data/ATAC_seq_shFOXA1_snakepipes_results/MACS_allLen_noDUP_noChrM/NC-T1/NC-T1_peaks.narrowPeak hg38 NC-T1.posStrand.sizeGiven -size given -norevopp -mask -p 20 -S 15 &> NC-T1.posStrand.sizeGiven.log
findMotifsGenome.pl /Users/jplab/Desktop/snakepipes_Data/ATAC_seq_shFOXA1_snakepipes_results/MACS_allLen_noDUP_noChrM/NC-T2/NC-T2_peaks.narrowPeak hg38 NC-T2.posStrand.sizeGiven -size given -norevopp -mask -p 20 -S 15 &> NC-T2.posStrand.sizeGiven.log

# only pos strand and default size and use narrowpeak
findMotifsGenome.pl /Data/Projects/ATAC_seq_shFOXA1_20220302_snakepipes/MACS_allLen_noDUP_noChrM/F-B1/F-B1_peaks.narrowPeak hg38 F-B1.posStrand.narrowPeak.sizedefault -norevopp -mask -p 20 -S 15 &> F-B1.posStrand.narrowPeak.sizedefault.log
findMotifsGenome.pl /Data/Projects/ATAC_seq_shFOXA1_20220302_snakepipes/MACS_allLen_noDUP_noChrM/F-B2/F-B2_peaks.narrowPeak hg38 F-B2.posStrand.narrowPeak.sizedefault -norevopp -mask -p 20 -S 15 &> F-B2.posStrand.narrowPeak.sizedefault.log
findMotifsGenome.pl /Data/Projects/ATAC_seq_shFOXA1_20220302_snakepipes/MACS_allLen_noDUP_noChrM/F-D1/F-D1_peaks.narrowPeak hg38 F-D1.posStrand.narrowPeak.sizedefault -norevopp -mask -p 20 -S 15 &> F-D1.posStrand.narrowPeak.sizedefault.log
findMotifsGenome.pl /Data/Projects/ATAC_seq_shFOXA1_20220302_snakepipes/MACS_allLen_noDUP_noChrM/F-D2/F-D2_peaks.narrowPeak hg38 F-D2.posStrand.narrowPeak.sizedefault -norevopp -mask -p 20 -S 15 &> F-D2.posStrand.narrowPeak.sizedefault.log
findMotifsGenome.pl /Data/Projects/ATAC_seq_shFOXA1_20220302_snakepipes/MACS_allLen_noDUP_noChrM/F-T1/F-T1_peaks.narrowPeak hg38 F-T1.posStrand.narrowPeak.sizedefault -norevopp -mask -p 20 -S 15 &> F-T1.posStrand.narrowPeak.sizedefault.log
findMotifsGenome.pl /Data/Projects/ATAC_seq_shFOXA1_20220302_snakepipes/MACS_allLen_noDUP_noChrM/F-T2/F-T2_peaks.narrowPeak hg38 F-T2.posStrand.narrowPeak.sizedefault -norevopp -mask -p 20 -S 15 &> F-T2.posStrand.narrowPeak.sizedefault.log
findMotifsGenome.pl /Data/Projects/ATAC_seq_shFOXA1_20220302_snakepipes/MACS_allLen_noDUP_noChrM/NC-B1/NC-B1_peaks.narrowPeak hg38 NC-B1.posStrand.narrowPeak.sizedefault -norevopp -mask -p 20 -S 15 &> NC-B1.posStrand.narrowPeak.sizedefault.log
findMotifsGenome.pl /Data/Projects/ATAC_seq_shFOXA1_20220302_snakepipes/MACS_allLen_noDUP_noChrM/NC-B2/NC-B2_peaks.narrowPeak hg38 NC-B2.posStrand.narrowPeak.sizedefault -norevopp -mask -p 20 -S 15 &> NC-B2.posStrand.narrowPeak.sizedefault.log
findMotifsGenome.pl /Data/Projects/ATAC_seq_shFOXA1_20220302_snakepipes/MACS_allLen_noDUP_noChrM/NC-D1/NC-D1_peaks.narrowPeak hg38 NC-D1.posStrand.narrowPeak.sizedefault -norevopp -mask -p 20 -S 15 &> NC-D1.posStrand.narrowPeak.sizedefault.log
findMotifsGenome.pl /Data/Projects/ATAC_seq_shFOXA1_20220302_snakepipes/MACS_allLen_noDUP_noChrM/NC-D2/NC-D2_peaks.narrowPeak hg38 NC-D2.posStrand.narrowPeak.sizedefault -norevopp -mask -p 20 -S 15 &> NC-D2.posStrand.narrowPeak.sizedefault.log
findMotifsGenome.pl /Data/Projects/ATAC_seq_shFOXA1_20220302_snakepipes/MACS_allLen_noDUP_noChrM/NC-T1/NC-T1_peaks.narrowPeak hg38 NC-T1.posStrand.narrowPeak.sizedefault -norevopp -mask -p 20 -S 15 &> NC-T1.posStrand.narrowPeak.sizedefault.log
findMotifsGenome.pl /Data/Projects/ATAC_seq_shFOXA1_20220302_snakepipes/MACS_allLen_noDUP_noChrM/NC-T2/NC-T2_peaks.narrowPeak hg38 NC-T2.posStrand.narrowPeak.sizedefault -norevopp -mask -p 20 -S 15 &> NC-T2.posStrand.narrowPeak.sizedefault.log

# both strands and default size
findMotifsGenome.pl /Users/jplab/Desktop/snakepipes_Data/ATAC_seq_shFOXA1_snakepipes_results/MACS_allLen_noDUP_noChrM/F-B1/F-B1_summits.bed hg38 F-B1.bothStrands.sizedefault  -mask -p 20 -S 15 &> F-B1.bothStrands.sizedefault.log
findMotifsGenome.pl /Users/jplab/Desktop/snakepipes_Data/ATAC_seq_shFOXA1_snakepipes_results/MACS_allLen_noDUP_noChrM/F-B2/F-B2_summits.bed hg38 F-B2.bothStrands.sizedefault  -mask -p 20 -S 15 &> F-B2.bothStrands.sizedefault.log
findMotifsGenome.pl /Users/jplab/Desktop/snakepipes_Data/ATAC_seq_shFOXA1_snakepipes_results/MACS_allLen_noDUP_noChrM/F-D1/F-D1_summits.bed hg38 F-D1.bothStrands.sizedefault  -mask -p 20 -S 15 &> F-D1.bothStrands.sizedefault.log
findMotifsGenome.pl /Users/jplab/Desktop/snakepipes_Data/ATAC_seq_shFOXA1_snakepipes_results/MACS_allLen_noDUP_noChrM/F-D2/F-D2_summits.bed hg38 F-D2.bothStrands.sizedefault  -mask -p 20 -S 15 &> F-D2.bothStrands.sizedefault.log
findMotifsGenome.pl /Users/jplab/Desktop/snakepipes_Data/ATAC_seq_shFOXA1_snakepipes_results/MACS_allLen_noDUP_noChrM/F-T1/F-T1_summits.bed hg38 F-T1.bothStrands.sizedefault  -mask -p 20 -S 15 &> F-T1.bothStrands.sizedefault.log
findMotifsGenome.pl /Users/jplab/Desktop/snakepipes_Data/ATAC_seq_shFOXA1_snakepipes_results/MACS_allLen_noDUP_noChrM/F-T2/F-T2_summits.bed hg38 F-T2.bothStrands.sizedefault  -mask -p 20 -S 15 &> F-T2.bothStrands.sizedefault.log
findMotifsGenome.pl /Users/jplab/Desktop/snakepipes_Data/ATAC_seq_shFOXA1_snakepipes_results/MACS_allLen_noDUP_noChrM/NC-B1/NC-B1_summits.bed hg38 NC-B1.bothStrands.sizedefault  -mask -p 20 -S 15 &> NC-B1.bothStrands.sizedefault.log
findMotifsGenome.pl /Users/jplab/Desktop/snakepipes_Data/ATAC_seq_shFOXA1_snakepipes_results/MACS_allLen_noDUP_noChrM/NC-B2/NC-B2_summits.bed hg38 NC-B2.bothStrands.sizedefault  -mask -p 20 -S 15 &> NC-B2.bothStrands.sizedefault.log
findMotifsGenome.pl /Users/jplab/Desktop/snakepipes_Data/ATAC_seq_shFOXA1_snakepipes_results/MACS_allLen_noDUP_noChrM/NC-D1/NC-D1_summits.bed hg38 NC-D1.bothStrands.sizedefault  -mask -p 20 -S 15 &> NC-D1.bothStrands.sizedefault.log
findMotifsGenome.pl /Users/jplab/Desktop/snakepipes_Data/ATAC_seq_shFOXA1_snakepipes_results/MACS_allLen_noDUP_noChrM/NC-D2/NC-D2_summits.bed hg38 NC-D2.bothStrands.sizedefault  -mask -p 20 -S 15 &> NC-D2.bothStrands.sizedefault.log
findMotifsGenome.pl /Users/jplab/Desktop/snakepipes_Data/ATAC_seq_shFOXA1_snakepipes_results/MACS_allLen_noDUP_noChrM/NC-T1/NC-T1_summits.bed hg38 NC-T1.bothStrands.sizedefault  -mask -p 20 -S 15 &> NC-T1.bothStrands.sizedefault.log
findMotifsGenome.pl /Users/jplab/Desktop/snakepipes_Data/ATAC_seq_shFOXA1_snakepipes_results/MACS_allLen_noDUP_noChrM/NC-T2/NC-T2_summits.bed hg38 NC-T2.bothStrands.sizedefault  -mask -p 20 -S 15 &> NC-T2.bothStrands.sizedefault.log

# both strands and given size
findMotifsGenome.pl /Users/jplab/Desktop/snakepipes_Data/ATAC_seq_shFOXA1_snakepipes_results/MACS_allLen_noDUP_noChrM/F-B1/F-B1_peaks.narrowPeak hg38 F-B1.bothStrands.sizeGiven -size given  -mask -p 20 -S 15 &> F-B1.bothStrands.sizeGiven.log
findMotifsGenome.pl /Users/jplab/Desktop/snakepipes_Data/ATAC_seq_shFOXA1_snakepipes_results/MACS_allLen_noDUP_noChrM/F-B2/F-B2_peaks.narrowPeak hg38 F-B2.bothStrands.sizeGiven -size given  -mask -p 20 -S 15 &> F-B2.bothStrands.sizeGiven.log
findMotifsGenome.pl /Users/jplab/Desktop/snakepipes_Data/ATAC_seq_shFOXA1_snakepipes_results/MACS_allLen_noDUP_noChrM/F-D1/F-D1_peaks.narrowPeak hg38 F-D1.bothStrands.sizeGiven -size given  -mask -p 20 -S 15 &> F-D1.bothStrands.sizeGiven.log
findMotifsGenome.pl /Users/jplab/Desktop/snakepipes_Data/ATAC_seq_shFOXA1_snakepipes_results/MACS_allLen_noDUP_noChrM/F-D2/F-D2_peaks.narrowPeak hg38 F-D2.bothStrands.sizeGiven -size given  -mask -p 20 -S 15 &> F-D2.bothStrands.sizeGiven.log
findMotifsGenome.pl /Users/jplab/Desktop/snakepipes_Data/ATAC_seq_shFOXA1_snakepipes_results/MACS_allLen_noDUP_noChrM/F-T1/F-T1_peaks.narrowPeak hg38 F-T1.bothStrands.sizeGiven -size given  -mask -p 20 -S 15 &> F-T1.bothStrands.sizeGiven.log
findMotifsGenome.pl /Users/jplab/Desktop/snakepipes_Data/ATAC_seq_shFOXA1_snakepipes_results/MACS_allLen_noDUP_noChrM/F-T2/F-T2_peaks.narrowPeak hg38 F-T2.bothStrands.sizeGiven -size given  -mask -p 20 -S 15 &> F-T2.bothStrands.sizeGiven.log
findMotifsGenome.pl /Users/jplab/Desktop/snakepipes_Data/ATAC_seq_shFOXA1_snakepipes_results/MACS_allLen_noDUP_noChrM/NC-B1/NC-B1_peaks.narrowPeak hg38 NC-B1.bothStrands.sizeGiven -size given  -mask -p 20 -S 15 &> NC-B1.bothStrands.sizeGiven.log
findMotifsGenome.pl /Users/jplab/Desktop/snakepipes_Data/ATAC_seq_shFOXA1_snakepipes_results/MACS_allLen_noDUP_noChrM/NC-B2/NC-B2_peaks.narrowPeak hg38 NC-B2.bothStrands.sizeGiven -size given  -mask -p 20 -S 15 &> NC-B2.bothStrands.sizeGiven.log
findMotifsGenome.pl /Users/jplab/Desktop/snakepipes_Data/ATAC_seq_shFOXA1_snakepipes_results/MACS_allLen_noDUP_noChrM/NC-D1/NC-D1_peaks.narrowPeak hg38 NC-D1.bothStrands.sizeGiven -size given  -mask -p 20 -S 15 &> NC-D1.bothStrands.sizeGiven.log
findMotifsGenome.pl /Users/jplab/Desktop/snakepipes_Data/ATAC_seq_shFOXA1_snakepipes_results/MACS_allLen_noDUP_noChrM/NC-D2/NC-D2_peaks.narrowPeak hg38 NC-D2.bothStrands.sizeGiven -size given  -mask -p 20 -S 15 &> NC-D2.bothStrands.sizeGiven.log
findMotifsGenome.pl /Users/jplab/Desktop/snakepipes_Data/ATAC_seq_shFOXA1_snakepipes_results/MACS_allLen_noDUP_noChrM/NC-T1/NC-T1_peaks.narrowPeak hg38 NC-T1.bothStrands.sizeGiven -size given  -mask -p 20 -S 15 &> NC-T1.bothStrands.sizeGiven.log
findMotifsGenome.pl /Users/jplab/Desktop/snakepipes_Data/ATAC_seq_shFOXA1_snakepipes_results/MACS_allLen_noDUP_noChrM/NC-T2/NC-T2_peaks.narrowPeak hg38 NC-T2.bothStrands.sizeGiven -size given  -mask -p 20 -S 15 &> NC-T2.bothStrands.sizeGiven.log

