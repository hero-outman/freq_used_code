from glob import glob
import sys, os,uuid
from glbase3 import *

config.draw_mode = "pdf"

# load expression .glb file
# expn = glload('./genes_ntc_expression.glb')

# # min_expression, number_of_conditions ： 
# # filter genes by a minimum_expression value in at least number_of_conditions
# expn.filter_low_expressed(100, 2)

# expn.strip_errs()

# # log for high-low expression
# expn.log(2,.1)

# # convert to zscore
# expn.row_Z()

# # ['A8301_2.0', 'B14_2.0', 'control_2.0', 'Eli_2.0', 'TGFb_2.0']
# print(expn.getConditionNames())
# print('expn',expn)

# should i do this 
# expn = expn.norm_multi_fc({'B14_2.0.rp1': ['B14_2.0.rp2']})


# zscore
expn = expression(
    filename='/Users/jplab/Desktop/2021-06_data/B14_DMSO_TGF_Rawcounts_Removelowexpr.tsv', 
    format={'force_tsv': True, 'ensg': 0, 'skiplines':0}, expn='column[1:]',
    cond_names=[
    'mean_B14_2.0',
    'mean_DMSO_2.0',
    'mean_TGFb_2.0',
    ]
)

print(expn)


# return a copy of the expression-data, but only containing 
# the condition names specified in conditions
expn = expn.sliceConditions(
    [
        'mean_B14_2.0',
        'mean_DMSO_2.0',
        'mean_TGFb_2.0'
    ]
)


# return a copy of the expression-data, but only containing 
# the condition names specified in conditions
# expn = expn.sliceConditions(
#     [
#         'B14_2.0',
#         'control_2.0',
#         'TGFb_2.0'
#     ]
# )

# todo : what's 200 meaning, how to set conditions ?
# expn.filter_low_expressed(200, 2)

# hg38
hg38_path = '/Users/jplab/Desktop/cut_tag_20210611/filter_region/hg38_ensembl_v95_enst.glb'
hg38 = glload(hg38_path)

# use diff range to find gene and peak
range = [10000,20000,50000]
for r in range:
    print('current range is ',r)

    # B14 open, Control close, TGF close clus bed
    clus_name = 'B14'
    # peaks0 = genelist('./cid_1.bed', format=format.minimal_bed) 
    # peaks1 = genelist('./B14_2.bed', format=format.minimal_bed) 
    peaks2 = genelist('/Users/jplab/Desktop/cut_tag_20210611/peak/narrowpeak/B14.FOXA1.DMSO.TGF.Uniq.bed', format=format.bed) 


    peaks2 = hg38.annotate(
        genelist=peaks2, 
        # distance=2000, 
        distance=r, 
        closest_only=False,
        image_filename = 'cid_1_%s_tss.png' % r
        )
    peaks2 = peaks2.removeDuplicates('ensg')

    # peaks1 = hg38.annotate(genelist=peaks1, distance=2000, closest_only=False)
    # peaks1 = peaks1.removeDuplicates('ensg')

    # peaks2 = hg38.annotate(genelist=peaks2, distance=2000, closest_only=False)
    # peaks2 = peaks2.removeDuplicates('ensg')

    peak_list = [
        # peaks0,
        # peaks1,
        peaks2
        ]

    expn.sort_sum_expression()


    for gl in peak_list:
        o = peaks2.frequencyAgainstArray(
            filename="%s-%s.pdf" % (gl.name, r),
            expression=expn,
            match_key="ensg",
            # bracket = None, wont work, throw an error
            # bracket=[4, 16],
            bracket=[-3, 3],
            spline_interpolate=False,  # right plot type
            # step_style = True,
            # draw_frames=True,
            window=2000,
            # size=[10,25],
        )

        # print('o',o)

        # use diff range to find overlap
        o.saveTSV("{0}_peak_overlap_{1}.tsv".format(clus_name,r))

