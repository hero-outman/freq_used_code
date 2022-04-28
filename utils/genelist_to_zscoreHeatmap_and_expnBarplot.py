from cmath import log
from operator import ge
import pandas as pd
import numpy as np
from glbase3 import *
import os, subprocess,sys
from alive_progress import alive_bar

os.chdir('/Users/jplab/Desktop/DAILY_CODE_DATA/2022-4/data/4-27_FOX_family_mRNA_shnc')
config.draw_mode = "png"

markers = {
    # 'endosomal_enrty':['ACTR2','ACTR3','ARPC3','ARPC4','RAB7A','UVRAG','CCZ1B'],
    # 'spike_cleavage_and_membrance_fusion':['ATP6AP1','ATP6AP2','ATP6V1A','ATP6V1B2','ATP6V0B','ATP6V0C','ATP6V1C1','ATP6V1E1','ATP6V0D1','ATP6V1G1','TMEM199','ATP6V1H','CTSL','TOR1AIP1'],
    # 'endosome_recycling':['VPS26A','VPS29','VPS35','SNX27','PIK3C3','WDR81','ACP5','COMMD2','COMMD3','COMMD3-BMI1','COMMD4'],
    # 'Translation':['SLTM','SPEN'],
    # 'Golgi':['PPID','CHST14'],
    # 'endoplasmic_reticulum':['DPM3','ERMP1']
    'FOX family in mRNA_sh_nc_Data':{'FOXA1', 'FOXA2', 'FOXC1', 'FOXF1', 'FOXF2', 'FOXJ1', 'FOXJ2', 'FOXK2', 'FOXM1', 'FOXN1', 'FOXN3', 'FOXN4', 'FOXO1', 'FOXO3', 'FOXP1', 'FOXP2', 'FOXP3', 'FOXP4', 'FOXRED1', 'FOXRED2'}
}

num_terms = len(markers)

# B vs D vs T
deseq2Norm_counts = '/Users/jplab/Desktop/snakepipes_Data/results_mRNA_seq_FOXA1_KD_NC_2022221/featureCounts/counts_deseq2Normed_symbol.tsv'

df_deseq2Norm_counts = pd.read_csv(deseq2Norm_counts,header=0, sep='\t')

columns = df_deseq2Norm_counts.columns.to_list()
start_offset = columns.index('geneid')
end_offset = columns.index('baseMean')
gene_name_offset = columns.index('external_gene_name')
sample_names = columns[start_offset+1:end_offset]
sample_numbers = len(sample_names)

def create_folder_ifNotExists(folder_path):
    CHECK_FOLDER = os.path.isdir(folder_path)

    if not CHECK_FOLDER:
        os.makedirs(folder_path)
        print("created folder : ", folder_path)

    else:
        print(folder_path, "folder already exists.")

with alive_bar(num_terms) as bar:

    # match to markers
    for k,v in markers.items():

        # gene name upper case
        up = []
        for e in v:
            e = e.upper()
            up.append(e)

        marker_type = k
        marker_value = up

        marker_expn = df_deseq2Norm_counts[df_deseq2Norm_counts["external_gene_name"].isin(marker_value)]

        temp_fileName = '{}_marker_expn.temp.tsv'.format(marker_type)
        marker_expn.to_csv(temp_fileName,index=False,header=True,sep='\t')


        expn = expression(
                # loadable_list=marker_expn,
                filename= temp_fileName,
                format={'force_tsv': True, 'name':gene_name_offset, 'ensg':start_offset,'skiplines': 0}, 
                expn='column['+str(start_offset+1)+':'+str(end_offset)+']',
                # cond_names=[
                # # 'B1.14_0.4_rp1', 'B1.14_0.4_rp2', 'B1.14_0.8_rp1', 'B1.14_0.8_rp2', 'B1.14_1.2_rp1', 'B1.14_1.2_rp2', 'B1.14_1.6_rp1', 'B1.14_1.6_rp2'
                # ]
                )


        # expn.setConditionNames([
        # 'B1.14_NC',
        # 'DMSO_NC',
        # 'TGFbeta_NC',
        # 'B1.14_FOXA1_KD',
        # 'DMSO_FOXA1_KD',
        # 'TGFb_FOXA1_KD'
        # ])   

        # cond_order1 = [
        #     'B1.14_NC',
        #     'B1.14_FOXA1_KD',
        #     'DMSO_NC',
        #     'DMSO_FOXA1_KD',
        #     'TGFbeta_NC',
        #     'TGFb_FOXA1_KD'
        # ]

        # cond_order2 = [
        #     'B1.14_NC',
        #     'DMSO_NC',
        #     'TGFbeta_NC',
        #     'B1.14_FOXA1_KD',
        #     'DMSO_FOXA1_KD',
        #     'TGFb_FOXA1_KD'
        # ]  

        # expn1 = expn.sliceConditions(cond_order1) 
        # expn2 = expn.sliceConditions(cond_order2)

        gene_list = up
        expn_barplot = './'+marker_type+'_expnBarPlot/'
        create_folder_ifNotExists(expn_barplot)

        color_sample = ['#F08080','#82E0AA']
        colors_arrange = np.repeat(color_sample, len(sample_names)/2).tolist()
        for g in gene_list:
            expn.barh_single_item(value=g, key="name", filename=expn_barplot+"%s.png" % g, bar_cols=colors_arrange,
                # fold_change = True,
                size=(16, 12), vert_space=0.8, #tree=tree,
                title_fontsize=15, yticklabel_fontsize=15, xticklabel_fontsize=15)
      
        expn.drawBarChart(gene_symbols=gene_list, filename=expn_barplot+"%s.png" % 'genelist', key="name", bar_cols=colors_arrange,labels=None,size=(18, 15),vert_space=0.6)

        expn.coerce(int)

        # handle long list heatmap height issue
        expn_length = len(expn)
        print(expn_length)

        sample_name = marker_type
        title_sample_name = sample_name

        if len(sample_name) > 50:
            title_sample_name = sample_name[0:25] + '-' + '\n' + sample_name[25:50] + '-' + '\n' + sample_name[50:]
            
        
        fig_width = 15
        figsize=(fig_width, 25)
        heat_hei = 0.01 * expn_length
        grid=True
        if expn_length <= 50:
            heat_hei = 0.03 * expn_length
            figsize=(fig_width, 7)
        if expn_length >= 120:
            heat_hei = 0.0045 * expn_length
        if expn_length >= 150:
            heat_hei = 0.0045 * expn_length
            figsize=(fig_width, 26)
        if expn_length >= 250:
            heat_hei = 0.002 * expn_length
            figsize=(fig_width, 35)
        if expn_length >= 500:
            heat_hei = 0.0013 * expn_length
            figsize=(fig_width, 48)
            grid=False
        if expn_length >= 800: 
            heat_hei = 0.0009 * expn_length
            figsize=(fig_width, 50)
            grid=False
        if expn_length >= 2000: 
            heat_hei = 0.00025 * expn_length
            figsize=(fig_width, 200)
            grid=False 

        heatmap_filename1 = './heatmap/'+sample_name+'_heatmap.expression'
        create_folder_ifNotExists('./heatmap/')
        expn.heatmap(
        filename=heatmap_filename1,
        log=False,
        # bracket = (0.0, 2000 ), 
        heat_wid=0.25, 
        # heat_wid=0.03*len(expn.getConditionNames()),
        col_cluster=False, 
        heat_hei=heat_hei, 
        # heat_hei=0.008*len(expn),
        colbar_label="Deseq2 normed expression: \n"+title_sample_name, 
        row_cluster=True,
        optimal_ordering = True,
        figsize=figsize,
        # highlights=highlights,
        grid=grid
        )



        expn.row_Z()

        

        heatmap_filename2 = './heatmap/'+sample_name+'_heatmap'
        create_folder_ifNotExists('./heatmap/')
        expn.heatmap(
        filename=heatmap_filename2, 
        # bracket = (-2.0, 2.0), 
        heat_wid=0.25, 
        # heat_wid=0.03*len(expn.getConditionNames()),
        col_cluster=False, 
        heat_hei=heat_hei, 
        # heat_hei=0.008*len(expn),
        colbar_label="Z-score: \n"+title_sample_name, 
        row_cluster=True,
        optimal_ordering = True,
        figsize=figsize,
        # highlights=highlights,
        grid=grid
        )



        rm_cmd = 'rm {0}'.format(temp_fileName) 

        exit_status = subprocess.call(rm_cmd, shell=True)

        bar()
        
        if exit_status == 1:
            print('temp file "{0}" failed to remove. Oh no, Barbecueeeeee le....'.format(marker_type))
        if exit_status == 0:
            print('{0} zscore heatmap and expn barplot plot successfully, great!'.format(marker_type))
        
            