import sys, os, glob, gzip

from glbase3 import *

'''
1. read expression from genes which have TF target squence
2. trans expression  value to Z-score
'''

# get expression obj from tsv
# esng: geneid index, start from 0
# expn: expresion value data column index, start from 0
# skiplines : skip rows, index, start from 0
# cond_names : 
# expn1 = expression(
#     filename='TGFbopen.10kb.anno.withExpr.tsv', 
#     format={'force_tsv': True, 'ensg': 19, 'skiplines': 0}, expn='column[21:]',
#     cond_names=[
#     'mean_B14_2.0',
#     'mean_control_2.0',
#     'mean_TGFb_2.0',
#     ]
# )

expn1 = expression(
    filename='/Users/jplab/Desktop/cut_tag_20210611/filter_region/6-22/Foxa1_B14_ATACB14_K27ac_K4me_anno.xlsx', 
    format={'force_tsv': True, 'ensg': 20, 'skiplines':1}, expn='column[21:26]',
    cond_names=[
    'mean_B14_2.0',
    'mean_control_2.0',
    'mean_TGFb_2.0',
    ]
)
print('expn1',expn1)

# a method in genelist: https://github.com/oaxiom/glbase3/blob/37b8ddf5a48e32354369893ec349e60d4ccb7ca1/genelist.py
'''
        **Purpose**
            remove the duplicates in the list and returns a new list;
            keeps the first example it finds
            This will only delete duplicates within the 'key'. For example,
            these three entries in a genelist:
            1: name: Stat3, score: 20, splicing: canonical
            2: name: Stat3, score: 30, splicing: alternate
            3: name: Nanog, score: 40, splicing: alternate
            gl = gl.removeDuplicates("name")
            will give:
            1: name: Stat3, score: 20, splicing: canonical
            3: name: Nanog, score: 40, splicing: alternate
            whilst
            gl = gl.removeDuplicates("splicing")
            will result in:
            1: name: Stat3, score: 20, splicing: canonical
            2: name: Stat3, score: 30, splicing: alternate
        **Arguments**
            key
                The key in which to make search for duplicates.
        **Returns**
            The new list with the duplicates removed.
'''

expn1 = expn1.removeDuplicates('ensg')

# float to integer
expn1.coerce(int)

# first log2 transform
# second zscore
expn1.log(2,.1)

# trans to z-score
expn1.row_Z()

# save to tsv for plotting 
# expn1.saveTSV("foxa1_TGFbopen_10kb_zscore.tsv")
expn1.saveTSV("SMAD2_TGFbopen_10kb_log_zscore.tsv")