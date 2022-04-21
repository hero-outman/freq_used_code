"""
code from post_norm, for reps merging
"""

import sys, os, glob
from glbase3 import *
config.draw_mode = ['pdf']

# gc normalized: 2umol_rawtags_gc_normed.tsv
raw_expn = expression(filename="/Users/jplab/Desktop/rna_seq_new/deg/B14_TGFb_normalized_by_DEG_log.tsv", format={"force_tsv": True, "skiplines": 0, "ensg": 0}, expn="column[1:]")

# Sort the Condition Names
raw_expn.sort_conditions()

# user_path = os.path.expanduser("~")
# ensg = glload(os.path.join(user_path, "hg38", "hg38_ensembl_v95_ensg_tes.glb"))

ensg = glload('../hg38_ensembl_v95_ensg_tes.glb')

# conditions: A8301_2.0.rp1	A8301_2.0.rp2	B14_2.0.rp1	B14_2.0.rp2	control.rp1	control.rp2	Eli_2.0.rp1	Eli_2.0.rp2	TGFb.rp1	TGFb.rp2
# mean_replicates will merge reps to there mean value
arr = raw_expn.mean_replicates(
    ['B14_2.0.rp1', 'B14_2.0.rp2'], ['TGFb.rp1', 'TGFb.rp2'],
    output_pears="2umol_Pearson_correlations.tsv",
    pearson_hist="2umol_Pearson_hist.png",
    threshold=0.8)

#cond_order = ['hESC shLuc', 'hESC shHnrnpu#1', 'hESC shHnrnpu#2']
#arr.setConditionNames(shared.pretify_sample_names(arr.getConditionNames()))
#print(arr.getConditionNames())
#arr = arr.sliceConditions(cond_order)

print("Num_conditions = ", len(arr.getConditionNames()))

print('arr conditoins : ', arr.getConditionNames())

# print('arr method : ', dir(arr))

# after merge, conditon name seems not proper, set the new one
# old names: 'B14_2.0.rp1', 'control.rp1','TGFb.rp1'
# new names to set :'B14_2.0', 'control_2.0', 'TGFb_2.0'
arr.setConditionNames(['B14_2.0', 'TGFb_2.0'])

mapped = ensg.map(key="ensg", genelist=arr)
mapped.save("B14_TGFb_merged_normalized_by_DEG_log.glb") # These are the final tables.
mapped.saveTSV("B14_TGFb_merged_normalized_by_DEG_log.tsv")