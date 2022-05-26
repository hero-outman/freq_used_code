##############################################################################################################################
# part 2/2 of motif heatmap
# plot heatmap of merged motif pvalue

# Attention: path and folder can not contains: . 
# eg: '/Users/jplab/Desktop/DAILY_CODE_DATA/2022-4/data/4-13_ATAC.SHNC.motif.heatmap/', 
# this will bring issue that heatmap file name '.png or .pdf', and can not get heatmap saved.
###############################################################################################################################
from glbase3 import *
import matplotlib.cm as cm
config.draw_mode = "pdf"

merged_motif_pvalue_path = '/Users/jplab/Desktop/DAILY_CODE_DATA/2022-5/data/5-25_profileheatmap_motifHeatmap/merged_known_motifs_pvalue.tsv'
save_to_path = '/Users/jplab/Desktop/DAILY_CODE_DATA/2022-5/data/5-25_profileheatmap_motifHeatmap/'


exp = expression(filename=merged_motif_pvalue_path,
                 format={'force_tsv': True, 'skiplines': 0, 'motif': 0}, expn='column[1:]')
print(exp)
print(exp.getConditionNames())

# exp=exp.sliceConditions(['CD34NvsMEFvsESC-all-open_groups.bed', 'CD34NvsMEFvsESC-MEF-open_groups.bed', 'CD34NvsMEFvsESC-CD34neg-open_groups.bed', 'CD34NvsMEFvsESC-ESC-open_groups.bed','CD34NvsMEFvsESC-CD34negESC-open_groups.bed', 'CD34NvsMEFvsESC-CD34negMEF-open_groups.bed',  'CD34NvsMEFvsESC-ESCMEF-open_groups.bed'])

exp_not_sogood = exp.deepcopy()

good_list = []
for i in exp:
    # sum of log10(pvalue) of each row
    if sum(i['conditions']) > 50: #andrew use 100
        i['group'] = 'good'
        good_list.append(i)
exp.load_list(good_list)
print(exp)
exp.saveTSV(save_to_path+'all_merge_motif_with_goodLabel.tsv')

# all motif
# -log10(0.05) -> pvalue 0.05
draw_numbers_threshold_p005 = 1.30103
# -log10(0.01) -> pvalue 0.01
draw_numbers_threshold_p001 = 2
# -log10(10^-20) -> pvalue 1e-20
draw_numbers_threshold_p20 = 20
draw_numbers_threshold_p50 = 50
draw_numbers_threshold_p100 = 100
draw_numbers_threshold_p150 = 150
draw_numbers_threshold_p200 = 200

colbar_label = '-log10(pvalue)' + '\n' + '* indicates p < 1e-' + str(draw_numbers_threshold_p20) + '\n' + '\n' \
                'MOTIFS with P-value > 0.01 and ' + '\n' + 'target/background <1.5 and ' + '\n' + 'target < = 2% '+ '\n' +'will give them -log10(pvalue) to 0'

# exp_not_sogood.heatmap(filename=save_to_path+'all_merge_motif_pvalue_heatmap', row_label_key='motif', 
#             bracket=[-10,380],
#             cmap=cm.Reds, row_cluster=True, col_cluster=False, col_font_size=8,
#             draw_numbers=True, draw_numbers_fmt='*', draw_numbers_threshold=draw_numbers_threshold_p001, grid=True,
#             heat_wid=0.03*len(exp_not_sogood.getConditionNames()), heat_hei=0.002*len(exp_not_sogood), row_font_size=4,figsize=(9,20),
#             colbar_label=colbar_label
#             )

# good motif
exp.heatmap(filename=save_to_path+'all_merge_motif_pvalue_good_heatmap', row_label_key='motif', 
            bracket=[-10,380],
            cmap=cm.Reds, row_cluster=True, col_cluster=False, col_font_size=8,
            draw_numbers=True, draw_numbers_fmt='*', draw_numbers_threshold=draw_numbers_threshold_p001, grid=True,
            heat_wid=0.03*len(exp.getConditionNames()), heat_hei=0.008*len(exp), row_font_size=4,figsize=(9,15),
            colbar_label=colbar_label
            )
