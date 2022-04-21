import pandas as pd
import glob, os, subprocess, sys, logging
import gseapy as gp
from gseapy.plot import gseaplot
from alive_progress import alive_bar


# working dir
wk_dir = '/Users/jplab/Desktop/cw/GSEA/'
os.chdir(wk_dir)

# input file path
deg_for_gsea_file_path = '/Users/jplab/Desktop/cw/*txt'
deg_files = glob.glob(deg_for_gsea_file_path)

# set GSEA min gene set size
min_gene_set_size = 5
# GSEA method
def do_gseapy(deg_file_list, gmt_file):

    gsea_result_path = wk_dir
    statu_flag = False

    for index, f_deg in enumerate(deg_files):

        f_name = os.path.basename(f_deg).split('_gseapy_input.txt')[0]
        df_deg = pd.read_csv(f_deg, header=0, sep='\t').drop_duplicates(subset=['NAME'])
        gmt_file_name = os.path.basename(gmt_file)
        # check gmt file is a temp file or not; if is a temp file, then remove prefix 'temp_';
        # an temp example: temp_c5.all.v7.5.1.symbols.gmt

        if gmt_file_name.startswith('temp_'): 
            gmt_file_name = gmt_file_name.split('_')[-1]
            # print('gmt name use for fold: {0}'.format(gmt_file_name))

        # this_cls = cls_file_list[index]
        this_cls = df_deg.columns.values.tolist()

        # 1 is the maximum number of splits to perform
        # this_cls = [i.split('_1', 1)[0] for i in this_cls]
        # this_cls = [i.split('_2', 1)[0] for i in this_cls]
        this_cls = [i.split('1', 1)[0] for i in this_cls]
        this_cls = [i.split('2', 1)[0] for i in this_cls]
        this_cls = this_cls[2:]
        print('this_cls: {0}, using gmt file: {1}, using input: {2}'.format(this_cls, gmt_file_name, f_name))

        # print(f'deg file name:{df_deg}, cls file name: {this_cls}')

        # codes for debug
        # with open(this_cls) as c:
        #     file = c.readlines()
        #     sample_name = file[1].lstrip("# ").strip('\n').split(" ")
        # print(sample_name)
        # print(len(sample_name))

        try:
            gs_res = gp.gsea(data=df_deg,
                            gene_sets=gmt_file,
                            min_size=min_gene_set_size,
                            max_size=500,
                            cls=this_cls,
                            outdir=gsea_result_path+gmt_file_name+'/' + f_name + '/',
                            no_plot=False,
                            # set permutation_type to phenotype if samples >=15
                            permutation_type='gene_set',
                            method='log2_ratio_of_classes',
                            permutation_num=1000,
                            processes=16, seed=7,
                            format='png',
                            figsize=[6.5, 6]
                            )
            statu_flag = True                
        except Exception as e:
            statu_flag = False
            print(str(e))
            sys.exit('Oh no, Barbecue le....')

    return gs_res, statu_flag            


# multifile
# gmt_file_path = '/Users/jplab/Desktop/freq_used_file/gene_sets/gmt/'
# gmt_files = glob.glob(gmt_file_path+'c*.v7.5.1.symbols.gmt')
gmt_files = [
#     '/Users/jplab/Desktop/freq_used_file/gene_sets/gmt/c2.all.v7.5.1.symbols.gmt',
#     '/Users/jplab/Desktop/freq_used_file/gene_sets/gmt/c7.all.v7.5.1.symbols.gmt',
]

# single file
# gmt_files = glob.glob('/Users/jplab/Desktop/freq_used_file/gene_sets/gmt/*.symbols.gmt')
gmt_files = glob.glob('/Users/jplab/Desktop/freq_used_file/gene_sets/gmt/h.all.v7.5.1.symbols.gmt')


for gmt_file in gmt_files:

    file_name = os.path.basename(gmt_file).split('_geneset.gmt')[0]

    # check how many lines in gmt file, if more than 2000, split it by 2000; or GSEApy will crash;
    num_lines = sum(1 for line in open(gmt_file))
    each_split_num = 500
    temp_file_format = 'temp_'+file_name

    if num_lines <= each_split_num:
        # print('gmt file: {0} is less than {1}, doing GSEApy now!'.format(gmt_file, each_split_num))

        gsea_res, status =  do_gseapy(deg_files, gmt_file)
        # status  = 'test less than {0}, do gsea; gmt file use: {1}'.format(each_split_num,gmt_file)
        # print(status)
        # print(gsea_res)

        if status :
            print('GSEApy successed on {0} gmt!'.format(gmt_file))
        else:
            print('GSEApy failed on {0} gmt!'.format(gmt_file))
    else:
        print('GSEApy will crash, I will split the gmt file by {0} lines or less'.format(each_split_num))

        split_cmd = 'split -l {0} {1} {2}'.format(each_split_num, gmt_file,temp_file_format) 

        exit_status = subprocess.call(split_cmd, shell=True)
    
        if exit_status == 1:
            print('large gmt file "{0}" failed to split. Oh no, Barbecueeeeee le....'.format(gmt_file))
        if exit_status == 0:
            print('large gmt file "{0}" was successed split into {1} lines each file, great!'.format(gmt_file,each_split_num))

        # split command will add nasty suffix to file extension, so need to reformat it
        # os.rename('a.txt', 'b.kml')
        temp_gmt_files = glob.glob('./'+temp_file_format+'??')
        for temp_gmt in temp_gmt_files:
            nasty_suffix = temp_gmt.split(file_name)[-1]
            # print('nasty suffix: {0}'.format(nasty_suffix))
            os.rename(temp_gmt, 'temp_'+nasty_suffix+'_'+file_name)
            # print('new file name: {0}'.format('temp_'+nasty_suffix+'_'+file_name))

        # after split gmt files to temp file, perform gsea analysis
        temp_gmt_files = glob.glob('./temp_*.gmt')

        total_temp_gmt_files = len(temp_gmt_files)

        logging.basicConfig(level=logging.WARN)
        logger = logging.getLogger('now handling {0} gmt files'.format(gmt_file))

        # do gsea with progress bar
        with alive_bar(total_temp_gmt_files) as bar:
        
            for temp_gmt in temp_gmt_files:
                gsea_res, status =  do_gseapy(deg_files, temp_gmt)
                # status  = 'test more than {0}, do gsea; gmt file use: {1}'.format(each_split_num,temp_gmt)
                # print(status)

                bar() 

                if status == True :
                    print('GSEApy successed on {0} gmt!'.format(temp_gmt))
                else:
                    print('GSEApy failed on {0} gmt!'.format(temp_gmt))
   

        # # rm temp_ gmt_file
        rm_temp_cmd = 'rm {0}'.format('./temp_*') 
        rm_exit_status = subprocess.call(rm_temp_cmd, shell=True)
        if exit_status == 1:
            print('temp gmt file failed to delete')
        if exit_status == 0:
            print('temp gmt files removed')    