library(hgu133a.db)
library(AnnotationDbi)
library(dplyr)


df_B14_sh_up <- read.table("/Users/jplab/Desktop/2021-12/affymetrix_to_ensg/TGF_kd_up.tsv", header=FALSE, sep ="\t")
df_B14_sh_down <- read.table("/Users/jplab/Desktop/2021-12/affymetrix_to_ensg/TGF_kd_down.tsv", header=FALSE, sep ="\t")

# convert
df_B14_sh_up$probe_id <- mapIds(hgu133a.db,
                                keys=dplyr::pull(df_B14_sh_up, 1),
                                column="PROBEID",
                                keytype="SYMBOL",
                                multiVals="first")

df_B14_sh_down$probe_id <- mapIds(hgu133a.db,
                                keys=dplyr::pull(df_B14_sh_down, 1),
                                column="PROBEID",
                                keytype="SYMBOL",
                                multiVals="first")
# remove NA
df_B14_sh_up <- subset(df_B14_sh_up, !is.na(probe_id))
df_B14_sh_down <- subset(df_B14_sh_down, !is.na(probe_id))

# remove dup
df_B14_sh_up = df_B14_sh_up[!duplicated(df_B14_sh_up$probe_id),]
df_B14_sh_down = df_B14_sh_down[!duplicated(df_B14_sh_down$probe_id),]

write.table(df_B14_sh_up, file="TGF_sh_up_withProbeID.tsv",quote=F, sep="\t", row.names=T, col.names=T)
write.table(df_B14_sh_down, file="TGF_sh_down_withProbeID.tsv",quote=F, sep="\t", row.names=T, col.names=T)

