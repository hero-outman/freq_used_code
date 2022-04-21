# merge sample reps bam
samtools merge -@ 20 A8301.merged.bam A8301_a.bam A8301_b.bam
samtools merge -@ 20 B14_2.merged.bam B14_2_a.bam B14_2_b.bam
samtools merge -@ 20 control.merged.bam control_a.bam control_b.bam
samtools merge -@ 20 Eli_2.merged.bam Eli_2_a.bam Eli_2_b.bam
samtools merge -@ 20 TGF.merged.bam TGF_a.bam TGF_b.bam