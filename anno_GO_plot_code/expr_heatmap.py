# heatmap: row name display gene symbol
import glob, os, math
import matplotlib.colors as colors
from glbase3 import *
import matplotlib.cm as cm
from  matplotlib.colors import LinearSegmentedColormap


os.chdir('/Users/jplab/Desktop/2021-7/data/7-12_km_plot')

config.draw_mode = 'pdf'

expn = expression(
    # filename='/Users/jplab/Desktop/2021-7/data/7-9_heatmap/deg_metastasis_positive_fig.tsv' 
    filename='/Users/jplab/Desktop/2021-7/data/7-9_heatmap/deg_metastaic_neg.tsv', 
    format={'force_tsv': True, 'ensg': 0, 'skiplines':0,'name':7}, expn='[column[17], column[18],column[19],column[20],column[23],column[24]]',
    cond_names=[
    'B14_2.0.rep1','B14_2.0.rep2',
    'DMSO_2.0.rep1','DMSO_2.0.rep2',
    'TGFb.rep1','TGFbrep2'
    ]
)

expn.row_Z()

print(expn)

highlights_positive=['LGALS1','IFI27','IFI6','SH3BGRL3','TFF1','TNNT1','REG4','B2M','PRSS3','SEPW1','MT-CO1','MT-CO2','MT-CO3','MT-CYB','MUC5AC','CEACAM6','PRKCDBP','FN1','CLU','TSPAN8','DPYSL3','CXCL5','TAGLN2']
highlights_negative=['FOSB','JUNB','CASC15','RIMKLB','KRT7','NFKBIA','EFNA1','CCL20','CXCL1','CXCL2','CXCL8','ID3','ASS1','SLPI','FXYD2','CFD','RPS4Y1','KRT17']

# This dictionary defines the colormap
# cdict = {'red':  ((0.0, 0.0, 0.0),   # no red at 0
#                   (0.5, 1.0, 1.0),   # all channels set to 1.0 at 0.5 to create white
#                   (1.0, 0.8, 0.8)),  # set to 0.8 so its not too bright at 1

#         'green': ((0.0, 0.8, 0.8),   # set to 0.8 so its not too bright at 0
#                   (0.5, 1.0, 1.0),   # all channels set to 1.0 at 0.5 to create white
#                   (1.0, 0.0, 0.0)),  # no green at 1

#         'blue':  ((0.0, 0.0, 0.0),   # no blue at 0
#                   (0.5, 1.0, 1.0),   # all channels set to 1.0 at 0.5 to create white
#                   (1.0, 0.0, 0.0))   # no blue at 1
#        }

c = ["darkgreen","green","palegreen","black", "lightcoral","red","darkred"]
v = [0,.15,.4,.5,0.6,.9,1.]
l = list(zip(v,c))
cmap=LinearSegmentedColormap.from_list('rg',l, N=256)       

# Create the colormap using the dictionary
# GnRd = colors.LinearSegmentedColormap('GnRd', cdict)

res = expn.heatmap(
#     filename='positive_highlights',
    # filename='negative_highlights',
    filename='positive_zscore',
    # log=True,
    row_label_key='name',
    size=[9, 18], 
    # bracket=[-20,20],
    row_cluster=True, col_cluster=False, imshow=False,
    heat_wid=0.06, 
    # cmap=cm.GnRd,  # matplotlib's color map
    cmap=cmap,
    border=True,
    row_font_size=7, heat_hei=0.008*len(expn), grid=True,
    draw_numbers=False, 
    # draw_numbers_fmt='*',
    # draw_numbers_threshold=10**0.01, # 3 = 0.01,
    draw_numbers_font_size=6,
    # highlights=highlights_negative
    # highlights=highlights_positive
    ) # 1.30 = 0.05
