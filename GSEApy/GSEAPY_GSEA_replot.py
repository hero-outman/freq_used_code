import pandas as pd
import glob, os, subprocess, sys, logging
import gseapy as gp

from gseapy.plot import gseaplot

from alive_progress import alive_bar

wk_dir = '/Users/jplab/Desktop/DAILY_CODE_DATA/2022-3/data/3-18_GSEApy_replot/'
os.chdir(wk_dir)


gsea_edb_path = '/Users/jplab/gsea_home/output/mar17/FOXA1_B14_NCvsKD.Gsea.1647516108502/'

# run command inside python console
rep = gp.replot(
    indir=gsea_edb_path, 
    outdir=wk_dir,
    format='png',
    verbose=True)