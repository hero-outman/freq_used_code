from glob import glob
import sys, os,uuid
from glbase3 import *

config.draw_mode = "pdf"

# load expression .glb file
expn = glload('/Users/jplab/Desktop/cut_tag_20210611/filter_region/genes_ntc_expression_2021.glb')

# min_expression, number_of_conditions ： 
# filter genes by a minimum_expression value in at least number_of_conditions
expn.filter_low_expressed(100, 2)

expn.strip_errs()

# log for high-low expression
expn.log(2,.1)

# convert to zscore
expn.row_Z()

# ['A8301_2.0', 'B14_2.0', 'control_2.0', 'Eli_2.0', 'TGFb_2.0']
print(expn.getConditionNames())


# return a copy of the expression-data, but only containing 
# the condition names specified in conditions
expn = expn.sliceConditions(
    [
        'B14_2.0',
        'control_2.0',
        'TGFb_2.0'
    ]
)

print('expn',expn)

# hg38
hg38 = glload("/Users/jplab/Desktop/cut_tag_20210611/filter_region/hg38_ensembl_v95_enst.glb")

# # use diff range to find gene and peak
range = [2000,5000,10000,20000,50000]
for r in range:
    print('current range is ',r)

    # B14 open, Control close, TGF close clus bed
    clus_name = 'B14_DMSO_TGF'
    # clus_name = 'TGF'
   
    peaks = genelist('/Users/jplab/Desktop/cut_tag_20210611/filter_region/6-23/B14_FOXA1.ATACUniq.K27ac.bed', format=format.bed) 


    peaks = hg38.annotate(
        genelist=peaks, 
        distance=r, 
        # closest_only=False,
        closest_only=True,
        image_filename = 'B14_FOXA1.ATACUniq.K27ac_%s_tss.png' % r
        )
    peaks = peaks.removeDuplicates('ensg')

    peak_list = [
        peaks
        ]

    expn.sort_sum_expression()


    for gl in peak_list:
        o = gl.frequencyAgainstArray(
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

