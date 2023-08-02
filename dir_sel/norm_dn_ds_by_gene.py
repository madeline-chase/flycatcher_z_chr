import pandas as pd
import numpy as np
import argparse
from sys import argv

parser = argparse.ArgumentParser(description = 'Calculate pN/pS on a gene by gene basis')
parser.add_argument('-d', '--dn_ds', dest = 'dn_ds', help = 'Fasta formatted file containing coding sequences. CDS must be correct orientation.')
parser.add_argument('-o' '--out', dest = 'out', help = 'Bed12 formatted file containing coordinates of coding sequences. 0-based.')
args = parser.parse_args()


dn_ds_file = args.dn_ds
out_file = args.out

dn_ds = pd.read_csv(dn_ds_file, sep = '\t', names = ['gene', 'dn_s_count', 'dn_s_norm','dn_w_count','dn_w_norm', 'ds_s_count','ds_s_norm','ds_w_count','ds_w_norm'])


dn_ds['ds_GC_wt'] = dn_ds['ds_s_norm']/(dn_ds['ds_s_norm']+dn_ds['ds_w_norm'])
dn_ds['ds_AT_wt'] = dn_ds['ds_w_norm']/(dn_ds['ds_s_norm']+dn_ds['ds_w_norm'])

dn_ds['dn_GC_wt'] = dn_ds['dn_s_norm']/(dn_ds['dn_s_norm']+dn_ds['dn_w_norm'])
dn_ds['dn_AT_wt'] = dn_ds['dn_w_norm']/(dn_ds['dn_s_norm']+dn_ds['dn_w_norm'])

dn_ds['ds_GC'] = dn_ds['ds_s_count'] * dn_ds['ds_GC_wt']
dn_ds['ds_AT'] = dn_ds['ds_w_count'] * dn_ds['ds_AT_wt']

dn_ds['dn_GC'] = dn_ds['dn_s_count'] * dn_ds['dn_GC_wt']
dn_ds['dn_AT'] = dn_ds['dn_w_count'] * dn_ds['dn_AT_wt']

dn_ds['ds_tot'] = dn_ds['ds_GC'] + dn_ds['ds_AT']

dn_ds['dn_tot'] = dn_ds['dn_GC'] + dn_ds['dn_AT']

dn_ds['ds_norm_tot'] = dn_ds['ds_s_norm']+dn_ds['ds_w_norm']
dn_ds['dn_norm_tot'] = dn_ds['dn_s_norm']+dn_ds['dn_w_norm']

dn_ds['dn_gene'] = dn_ds['dn_tot']/dn_ds['dn_norm_tot']
dn_ds['ds_gene'] = dn_ds['ds_tot']/dn_ds['ds_norm_tot']

dn_ds['dn_ds_gene'] = dn_ds['dn_gene']/dn_ds['ds_gene']

dn_ds = dn_ds[dn_ds['dn_s_norm']>=0]
dn_ds = dn_ds[dn_ds['dn_w_norm']>=0]
dn_ds = dn_ds[dn_ds['ds_s_norm']>=0]
dn_ds = dn_ds[dn_ds['ds_w_norm']>=0]


dn_ds.to_csv(out_file, sep = '\t', index = False, na_rep='NULL', columns=['gene','ds_tot','dn_tot','ds_norm_tot','dn_norm_tot','dn_gene','ds_gene','dn_ds_gene'])