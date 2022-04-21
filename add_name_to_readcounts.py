"""
add gene symbol by hg38
"""
from os import chdir
from glbase3 import *

chdir('/Users/jplab/Desktop/RNA_seq_FOXA1KO_scramble_20211213/downstream/data/raw_count_flow/deg/')

config.draw_mode = ['pdf']

raw_expn = expression(filename="./FOXA1_scramble_deseq2Normed.tsv", format={"force_tsv": True, "skiplines": 0, "ensg": 0}, expn="column[1:]")
print(raw_expn.getConditionNames())

raw_expn.sort_conditions()

ensg = glload('/Users/jplab/Desktop/RNA_seq_flow_v20210830/RNA_downstream/ref_genome_file/hg38_ensembl_v95_ensg_tes.glb')

# mean_replicates will merge reps to there mean value
arr = raw_expn.mean_replicates(
        ['shFOXA1_B14.rep1', 'shFOXA1_B14.rep2'], ['shFOXA1_DMSO.rep1', 'shFOXA1_DMSO.rep2'], ['shFOXA1_TGFb.rep1', 'shFOXA1_TGFb.rep2'], 
        ['nc_B14.rep1', 'nc_B14.rep2'], ['nc_DMSO.rep1', 'nc_DMSO.rep2'], ['nc_TGFb.rep1', 'nc_TGFb.rep2'],
        threshold=0.8
        )      

cond_order = ['nc_B14.rep1', 'nc_DMSO.rep1', 'nc_TGFb.rep1',  'shFOXA1_B14.rep1', 'shFOXA1_DMSO.rep1', 'shFOXA1_TGFb.rep1']



#arr = arr.sliceConditions(cond_order)

# print("Num_conditions = ", len(arr.getConditionNames()))

# print('arr conditoins : ', arr.getConditionNames())

# # print('arr method : ', dir(arr))

# # after merge, conditon name seems not proper, set the new one
# # old names: 'A8301_2.0.rp1', 'B14_2.0.rp1', 'control.rp1', 'Eli_2.0.rp1', 'TGFb.rp1'
# # new names to set :'A8301_2.0', 'B14_2.0', 'control_2.0', 'Eli_2.0', 'TGFb_2.0'
arr.setConditionNames(['nc_B14', 'nc_DMSO', 'nc_TGFb',  'shFOXA1_B14', 'shFOXA1_DMSO', 'shFOXA1_TGFb'])
print(arr.getConditionNames())

# save reps merged
mapped = ensg.map(key="ensg", genelist=arr)
mapped.save("FOXA1KO_scramble_repMerged_withGeneName.glb") # These are the final tables.
mapped.saveTSV("FOXA1KO_scramble_repMerged_withGeneName.tsv")

# save gc_normed
mapped = ensg.map(key="ensg", genelist=raw_expn)
# mapped.setConditionNames([i.replace("_", " ") for i in mapped.getConditionNames()])
mapped.save("FOXA1KO_scramble_withGeneName.glb")
mapped.saveTSV("FOXA1KO_scramble_withGeneName.tsv")