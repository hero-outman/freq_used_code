
import glob, sys, os, math
from glbase3 import *
import matplotlib.cm as cm

os.chdir('/Users/jplab/Desktop/ATAC_data/annotation/GO')

config.draw_mode = 'pdf'

format = {'force_tsv': True, 'pvalue': 1, 'name': 0}

def get_clus_number(s):
    # A8301_2_a-ann20000kb-nonamedupes.tsv
    # A8301_2-ann2000kb-nonamedupes.tsv
    # return '_'.join(os.path.split(filename)[1].replace('-nonamedupes.tsv', '').split('_')[1:3])
    return '_'.join(os.path.split(filename)[1].replace('-nonamedupes.tsv', '').split('_')[0:])

main_cluster_membership = {}

for ont in ('BP', 'CC', 'MF'):
    # GO_filenames = list(sorted(glob.glob('tabs/tab{0}_*2000k*.tsv'.format(ont))))
    # cid_7-ann20000kb-nonamedupes.tsv
    GO_filenames = list(sorted(glob.glob('/Users/jplab/Desktop/ATAC_data/annotation/GO/tabs/tab{0}_*2000k*.tsv'.format(ont))))
    num_clusters = len(GO_filenames)
    clus_order = []

    go_store = {}
    main_cluster_membership = {}
    for clus_number, filename in enumerate(GO_filenames):
        go = glgo(filename=filename, format=format)
        if not go:
            continue

        clus_name = get_clus_number(filename)
        print(clus_number, clus_name)

        clus_order.append(clus_name)

        #go.sort('1')
        top5 = go[0:25]

        for item in top5:
            if item['pvalue'] < 0.01:
                if item['name'] not in go_store:
                    go_store[item['name']] = [-1] * num_clusters

                go_store[item['name']][clus_number-1] = -math.log10(item['pvalue'])

    # fill in the holes:
    for clus_number, filename in enumerate(GO_filenames):
        go = glgo(filename=filename, format=format)
        if not go:
            continue

        clus_name = get_clus_number(filename)
        for k in go_store:
            this_k = go.get(key='name', value=k, mode='lazy') # by default
            if this_k:
                #print(k, clus_number)
                go_store[k][clus_number-1] = -math.log10(float(this_k[0]['pvalue']))

    newe = []

    for k in go_store:
        newe.append({'name': k, 'conditions': go_store[k]})

    goex = expression(loadable_list=newe, cond_names=clus_order)

    #goex = goex.sliceConditions(clus_order)
    goex = goex.filter_low_expressed(1.9, 1)

    res = goex.heatmap(filename='atmap_big_%s.png' % ont,
        size=[9, 18], bracket=[1,5],
        row_cluster=True, col_cluster=False, imshow=False,
        heat_wid=0.06, cmap=cm.Reds, border=True,
        row_font_size=7, heat_hei=0.008*len(goex), grid=True,
        draw_numbers=True, draw_numbers_fmt='*',
        draw_numbers_threshold=10**0.01, # 3 = 0.01,
        draw_numbers_font_size=6) # 1.30 = 0.05


    print('\n'.join(reversed(res['reordered_rows'])))
