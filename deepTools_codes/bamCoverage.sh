bamCoverage \
        --bam FOXA1_B14.merged.sortpos.bam \
        --outFileName ./FOXA1_B14.merged.dedup.seqDeptNorm.bw \
        --binSize 1 \
        --normalizeUsing RPGC \
        --effectiveGenomeSize 2913022398 \
        --blackListFileName /Data/Ref_genome_anno/ENCODE_BLACK_LIST_REGION/hg38-blacklist.v2.bed \
        --extendReads \
        --smoothLength 3 \
        --centerReads \
        --ignoreDuplicates \
        -p 20 &


bamCoverage \
        --bam FOXA1_DMSO.merged.sortpos.bam \
        --outFileName ./FOXA1_DMSO.merged.dedup.seqDeptNorm.bw \
        --binSize 1 \
        --normalizeUsing RPGC \
        --effectiveGenomeSize 2913022398 \
        --blackListFileName /Data/Ref_genome_anno/ENCODE_BLACK_LIST_REGION/hg38-blacklist.v2.bed \
        --extendReads \
        --smoothLength 3 \
        --centerReads \
        --ignoreDuplicates \
        -p 20 &    

bamCoverage \
        --bam FOXA1_TGF_2.filtered.bam \
        --outFileName ./FOXA1_TGF.dedup.seqDeptNorm.bw \
        --binSize 1 \
        --normalizeUsing RPGC \
        --effectiveGenomeSize 2913022398 \
        --blackListFileName /Data/Ref_genome_anno/ENCODE_BLACK_LIST_REGION/hg38-blacklist.v2.bed \
        --extendReads \
        --smoothLength 3 \
        --centerReads \
        --ignoreDuplicates \
        -p 40 &

bamCoverage \
        --bam H3K27ac_B14.merged.sortpos.bam \
        --outFileName ./H3K27ac_B14.merged.dedup.seqDeptNorm.bw \
        --binSize 1 \
        --normalizeUsing RPGC \
        --effectiveGenomeSize 2913022398 \
        --blackListFileName /Data/Ref_genome_anno/ENCODE_BLACK_LIST_REGION/hg38-blacklist.v2.bed \
        --extendReads \
        --smoothLength 3 \
        --centerReads \
        --ignoreDuplicates \
        -p 20 &                    

bamCoverage \
        --bam H3K27ac_DMSO.merged.sortpos.bam \
        --outFileName ./H3K27ac_DMSO.merged.dedup.seqDeptNorm.bw \
        --binSize 1 \
        --normalizeUsing RPGC \
        --effectiveGenomeSize 2913022398 \
        --blackListFileName /Data/Ref_genome_anno/ENCODE_BLACK_LIST_REGION/hg38-blacklist.v2.bed \
        --extendReads \
        --smoothLength 3 \
        --centerReads \
        --ignoreDuplicates \
        -p 20 &   

bamCoverage \
        --bam H3K27ac_TGF.merged.sortpos.bam  \
        --outFileName ./H3K27ac_TGF.merged.dedup.seqDeptNorm.bw \
        --binSize 1 \
        --normalizeUsing RPGC \
        --effectiveGenomeSize 2913022398 \
        --blackListFileName /Data/Ref_genome_anno/ENCODE_BLACK_LIST_REGION/hg38-blacklist.v2.bed \
        --extendReads \
        --smoothLength 3 \
        --centerReads \
        --ignoreDuplicates \
        -p 20 &