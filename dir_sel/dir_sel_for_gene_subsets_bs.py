import pandas as pd
import random
import numpy as np
import argparse
from sys import argv
import os

parser = argparse.ArgumentParser(description = 'Estimate overall pN/pS, dN/dS and optionally alpha and omega a, for a given set of genes, as well as jackknife resampled standard error.')
parser.add_argument('-d', '--dn_ds', dest = 'dn_ds', help = 'Fasta formatted file containing coding sequences. CDS must be correct orientation.')
parser.add_argument('-p', '--pn_ps', dest = 'pn_ps', help = 'Fasta formatted file containing coding sequences. CDS must be correct orientation.')
parser.add_argument('-g','--genes',dest = 'genes', help = '')
parser.add_argument('-f4', dest='freq_file_four', help = 'File with derived allele counts for neutral sites')
parser.add_argument('-f0', dest='freq_file_zero', help = 'File with derived allele counts for selected sites')
parser.add_argument('-n', dest='n_alleles', help = 'Number of alleles sampled', type = int)
parser.add_argument('-s0', dest='conf_zero', help = '')
parser.add_argument('-s1', dest='conf_one', help = '')
parser.add_argument('-a', dest='conf_alpha', help = '')
parser.add_argument('--dfe', dest = 'dfe_path', help = '')
parser.add_argument('-b', dest = 'boot', help = 'Number of bootstrap iterations', type = int)
parser.add_argument('-o' '--out', dest = 'out', help = 'Prefix to use for outfiles')
args = parser.parse_args()

dn_ds_file = args.dn_ds
pn_ps_file = args.pn_ps
gene_file = args.genes
out_pre = args.out

freq_file_four = args.freq_file_four 
freq_file_zero = args.freq_file_zero 
n_alleles = args.n_alleles 
site0_conf = args.conf_zero 
site1_conf = args.conf_one 
dfe_path = args.dfe_path
alpha_conf = args.conf_alpha if args.conf_alpha is not None else None

boot = args.boot


dn_ds_df = pd.read_csv(dn_ds_file, sep = '\t')
pn_ps_df = pd.read_csv(pn_ps_file, sep = '\t')

gene_df = pd.read_csv(gene_file, names = ['gene'], sep = '\t')

gene_list_d = list(gene_df['gene'])
gene_list_p = list(gene_df['gene'])

dn_count_obs = np.sum(dn_ds_df['dn_tot'][dn_ds_df['gene'].isin(gene_list_d)])
dn_norm_obs = np.sum(dn_ds_df['dn_norm_tot'][dn_ds_df['gene'].isin(gene_list_d)])
dn_obs = dn_count_obs/dn_norm_obs

ds_count_obs = np.sum(dn_ds_df['ds_tot'][dn_ds_df['gene'].isin(gene_list_d)])
ds_norm_obs = np.sum(dn_ds_df['ds_norm_tot'][dn_ds_df['gene'].isin(gene_list_d)])
ds_obs = ds_count_obs/ds_norm_obs

dn_ds_obs = dn_obs/ds_obs

pn_obs = np.sum(pn_ps_df['het_zero'][pn_ps_df['gene'].isin(gene_list_p)])/np.sum(pn_ps_df['len_zero'][pn_ps_df['gene'].isin(gene_list_p)])
ps_obs = np.sum(pn_ps_df['het_four'][pn_ps_df['gene'].isin(gene_list_p)])/np.sum(pn_ps_df['len_four'][pn_ps_df['gene'].isin(gene_list_p)])
pn_ps_obs = pn_obs/ps_obs



if freq_file_four is not None:
    freq_counts_four = pd.read_csv(freq_file_four, sep = '\t', names = ['scaff','start','stop','n_chr','n_der','gene'])
    freq_counts_zero = pd.read_csv(freq_file_zero, sep = '\t', names = ['scaff','start','stop','n_chr','n_der','gene'])

    count_vec_four = []
    count_vec_zero = []

    four_counts_sub = freq_counts_four[freq_counts_four['gene'].isin(gene_list_p)]
    zero_counts_sub = freq_counts_zero[freq_counts_zero['gene'].isin(gene_list_p)]
    for i in range(0,n_alleles+1):
        count_vec_four.append(len(four_counts_sub[four_counts_sub['n_der']==i]))
        count_vec_zero.append(len(zero_counts_sub[zero_counts_sub['n_der']==i]))
    
    sfs_file = out_pre + ".sfs_for_dfe.txt"
    f = open(sfs_file, "w")
    f.write('1')
    f.write('\n')
    f.write(str(n_alleles))
    f.write('\n')
    f.write(' '.join(str(d) for d in count_vec_zero))
    f.write('\n')
    f.write(' '.join(str(d) for d in count_vec_four))
    f.close()

    sel_div = ['1', dn_norm_obs, dn_count_obs]
    neut_div = ['0', ds_norm_obs, ds_count_obs]

    div_file = out_pre + ".div_for_alpha.txt"
    f = open(div_file, "w")
    f.write(' '.join(str(d) for d in sel_div))
    f.write('\n')
    f.write(' '.join(str(d) for d in neut_div))
    f.write('\n')
    f.close()

    site0_cmd = '{0}/est_dfe -c {1}'.format(dfe_path,site0_conf)
    site1_cmd = '{0}/est_dfe -c {1}'.format(dfe_path,site1_conf)
    alpha_cmd = '{0}/est_alpha_omega -c {1}'.format(dfe_path,alpha_conf)

    os.system(site0_cmd)
    os.system(site1_cmd)
    os.system(alpha_cmd)

    rename_site0 = 'mv site_0_2_epoch/est_dfe.out dfe_alpha_out_files/{}.obs_dfe_site_0.txt'.format(out_pre)
    rename_site1 = 'mv site_1_2_epoch/est_dfe.out dfe_alpha_out_files/{}.obs_dfe_site_1.txt'.format(out_pre)
    rename_alpha = 'mv alph_2_epoch/alpha_omega.txt {}.obs_alph_omega.txt'.format(out_pre)
    move_div = 'mv {}.div_for_alpha.txt dfe_alpha_in_files/'.format(out_pre)
    move_sfs = 'mv {}.sfs_for_dfe.txt dfe_alpha_in_files/'.format(out_pre)

    os.system(rename_site0)
    os.system(rename_site1)
    os.system(rename_alpha)
    os.system(move_div)
    os.system(move_sfs)


pn_ps_boot_list = []
dn_ds_boot_list = []

gene_indices = gene_df.index
for iter in range(0,boot):
    rand_index = random.choices(gene_indices, k = len(gene_indices))

    boot_genes_d = list(gene_df.loc[rand_index,'gene'])
    boot_genes_p = list(gene_df.loc[rand_index,'gene'])

    dn_count = np.sum(dn_ds_df['dn_tot'][dn_ds_df['gene'].isin(boot_genes_d)])
    dn_norm = np.sum(dn_ds_df['dn_norm_tot'][dn_ds_df['gene'].isin(boot_genes_d)])
    
    ds_count = np.sum(dn_ds_df['ds_tot'][dn_ds_df['gene'].isin(boot_genes_d)])
    ds_norm = np.sum(dn_ds_df['ds_norm_tot'][dn_ds_df['gene'].isin(boot_genes_d)])
    
    ds_jack = ds_count/ds_norm
    dn_jack = dn_count/dn_norm
    dn_ds_jack = dn_jack/ds_jack
    
    dn_ds_boot_list.append([dn_count,dn_norm,ds_count,ds_norm,dn_jack,ds_jack,dn_ds_jack])

    ps_jack = np.sum(pn_ps_df['het_four'][pn_ps_df['gene'].isin(boot_genes_p)])/np.sum(pn_ps_df['len_four'][pn_ps_df['gene'].isin(boot_genes_p)])
    pn_jack = np.sum(pn_ps_df['het_zero'][pn_ps_df['gene'].isin(boot_genes_p)])/np.sum(pn_ps_df['len_zero'][pn_ps_df['gene'].isin(boot_genes_p)])
    pn_ps_jack = pn_jack/ps_jack
    
    pn_ps_boot_list.append([pn_jack,ps_jack,pn_ps_jack])

    if freq_file_four is not None:
        count_vec_four = []
        count_vec_zero = []

        four_counts_sub = freq_counts_four[freq_counts_four['gene'].isin(boot_genes_p)]
        zero_counts_sub = freq_counts_zero[freq_counts_zero['gene'].isin(boot_genes_p)]
        for i in range(0,n_alleles+1):
            count_vec_four.append(len(four_counts_sub[four_counts_sub['n_der']==i]))
            count_vec_zero.append(len(zero_counts_sub[zero_counts_sub['n_der']==i]))
        
        sfs_file = out_pre + ".sfs_for_dfe.txt"
        f = open(sfs_file, "w")
        f.write('1')
        f.write('\n')
        f.write(str(n_alleles))
        f.write('\n')
        f.write(' '.join(str(d) for d in count_vec_zero))
        f.write('\n')
        f.write(' '.join(str(d) for d in count_vec_four))
        f.close()

        sel_div = ['1', dn_norm, dn_count]
        neut_div = ['0', ds_norm, ds_count]

        div_file = out_pre + ".div_for_alpha.txt"
        f = open(div_file, "w")
        f.write(' '.join(str(d) for d in sel_div))
        f.write('\n')
        f.write(' '.join(str(d) for d in neut_div))
        f.write('\n')
        f.close()
        

        site0_cmd = '{0}/est_dfe -c {1}'.format(dfe_path,site0_conf)
        site1_cmd = '{0}/est_dfe -c {1}'.format(dfe_path,site1_conf)
        alpha_cmd = '{0}/est_alpha_omega -c {1}'.format(dfe_path,alpha_conf)

        os.system(site0_cmd)
        os.system(site1_cmd)
        os.system(alpha_cmd)

        rename_site0 = 'mv site_0_2_epoch/est_dfe.out {0}.{1}_boot_dfe_site_0.txt'.format(out_pre, iter)
        rename_site1 = 'mv site_1_2_epoch/est_dfe.out {0}.{1}_boot_dfe_site_1.txt'.format(out_pre, iter)
        rename_alpha = 'mv alph_2_epoch/alpha_omega.txt {0}.{1}_boot_alph_omega.txt'.format(out_pre, iter)

        os.system(rename_site0)
        os.system(rename_site1)
        os.system(rename_alpha)


dn_ds_boot_df = pd.DataFrame(dn_ds_boot_list, columns = ['dn_count','dn_norm','ds_count','ds_norm','dn_boot','ds_boot','dn_ds_boot'])
pn_ps_boot_df = pd.DataFrame(pn_ps_boot_list, columns = ['pn_boot','ps_boot','pn_ps_boot'])


lower_cut = round(len(dn_ds_boot_df)*0.025)-1
upper_cut = round(len(pn_ps_boot_df)*0.975)-1

dn_upper = dn_ds_boot_df['dn_boot'].sort_values(ignore_index=True)[upper_cut]
dn_lower = dn_ds_boot_df['dn_boot'].sort_values(ignore_index=True)[lower_cut]
ds_upper  = dn_ds_boot_df['ds_boot'].sort_values(ignore_index=True)[upper_cut]
ds_lower = dn_ds_boot_df['ds_boot'].sort_values(ignore_index=True)[lower_cut]
dn_ds_upper = dn_ds_boot_df['dn_ds_boot'].sort_values(ignore_index=True)[upper_cut]
dn_ds_lower = dn_ds_boot_df['dn_ds_boot'].sort_values(ignore_index=True)[lower_cut]


pn_upper = pn_ps_boot_df['pn_boot'].sort_values(ignore_index=True)[upper_cut]
pn_lower = pn_ps_boot_df['pn_boot'].sort_values(ignore_index=True)[lower_cut]
ps_upper = pn_ps_boot_df['ps_boot'].sort_values(ignore_index=True)[upper_cut]
ps_lower = pn_ps_boot_df['ps_boot'].sort_values(ignore_index=True)[lower_cut]
pn_ps_upper = pn_ps_boot_df['pn_ps_boot'].sort_values(ignore_index=True)[upper_cut]
pn_ps_lower = pn_ps_boot_df['pn_ps_boot'].sort_values(ignore_index=True)[lower_cut]

if freq_file_four is not None:


    join_boot_alphas = 'cat *boot_alph_omega.txt > {}.boot_alph_omega.all_concat.txt'.format(out_pre)
    join_boot_dfes = 'cat *boot_dfe_site_1.txt > {}.boot_dfe_site_1.all_concat.txt'.format(out_pre)
    os.system(join_boot_alphas)
    os.system(join_boot_dfes)

    boot_alpha_file = out_pre + '.boot_alph_omega.all_concat.txt'
    boot_dfe_file = out_pre + '.boot_dfe_site_1.all_concat.txt'

    boot_alphas = pd.read_csv(boot_alpha_file, sep = ' ', names =['l1','l2','s1','s2','a1','a2','o1','o2'])
    boot_dfes = pd.read_csv(boot_dfe_file, sep = ' ', names = ['N1','n1_val', 'n2', 'n2_val', 't2', 't2_val','nw','nw_val','b','bval','es','es_val','f0','f0_val','L','L_val'])
    obs_alpha_file = out_pre + '.obs_alph_omega.txt'
    obs_dfe_file = 'dfe_alpha_out_files/' + out_pre + '.obs_dfe_site_1.txt'

    obs_alpha = pd.read_csv(obs_alpha_file, sep = ' ', names =['l1','l2','s1','s2','a1','a2','o1','o2'])
    obs_dfe = pd.read_csv(obs_dfe_file, sep = ' ', names = ['N1','n1_val', 'n2', 'n2_val', 't2', 't2_val','nw','nw_val','b','bval','es','es_val','f0','f0_val','L','L_val'])
    obs_es = obs_dfe['es_val'][0]
    n2_est = obs_dfe['n2_val'][0]
    
    obs_alpha_val = obs_alpha['a2'][0]
    obs_omega_val = obs_alpha['o2'][0]

    obs_es = obs_es*n2_est

    alpha_upper = boot_alphas['a2'].sort_values(ignore_index=True)[upper_cut]
    alpha_lower = boot_alphas['a2'].sort_values(ignore_index=True)[lower_cut]

    omega_upper = boot_alphas['o2'].sort_values(ignore_index=True)[upper_cut]
    omega_lower = boot_alphas['o2'].sort_values(ignore_index=True)[lower_cut]

    Es_upper = boot_dfes['es_val'].sort_values(ignore_index=True)[upper_cut]*n2_est
    Es_lower = boot_dfes['es_val'].sort_values(ignore_index=True)[lower_cut]*n2_est


    clean_cmd = 'rm *boot_alph_omega.txt'
    os.system(clean_cmd)

    clean_cmd = 'rm *_boot_dfe_site_*.txt'
    os.system(clean_cmd)
    sel_out_file = out_pre + ".dir_sel_obs.txt"
    sel_out_head = ['dN','dN_lower','dN_upper', 'dS','dS_lower', 'dS_upper', 'dN_dS', 'dN_dS_lower', 'dN_dS_upper', 'pN','pN_lower', 'pN_upper', 'pS','pS_lower', 'pS_upper', 'pN_pS', 'pN_pS_lower', 'pN_pS_upper', 'alpha', 'alpha_lower', 'alpha_upper', 'omega_a',  'omega_a_lower', 'omega_a_upper', 'NeS',  'NeS_lower','NeS_upper']
    sel_out_vals = [dn_obs, dn_lower, dn_upper, ds_obs, ds_lower, ds_upper, dn_ds_obs,  dn_ds_lower, dn_ds_upper, pn_obs, pn_lower, pn_upper, ps_obs, ps_lower, ps_upper, pn_ps_obs, pn_ps_lower, pn_ps_upper,obs_alpha_val, alpha_lower, alpha_upper, obs_omega_val, omega_lower, omega_upper, obs_es,  Es_lower, Es_upper]
    f = open(sel_out_file, "w")
    f.write('\t'.join(sel_out_head))
    f.write('\n')
    f.write('\t'.join(str(d) for d in sel_out_vals))
    f.close()


else:

    sel_out_file = out_pre + ".dir_sel_obs.txt"
    sel_out_head = ['dN','dN_lower','dN_upper', 'dS','dS_lower', 'dS_upper', 'dN_dS', 'dN_dS_lower', 'dN_dS_upper', 'pN','pN_lower', 'pN_upper', 'pS','pS_lower', 'pS_upper', 'pN_pS', 'pN_pS_lower', 'pN_pS_upper']
    sel_out_vals = [dn_obs, dn_lower, dn_upper, ds_obs, ds_lower, ds_upper, dn_ds_obs,  dn_ds_lower, dn_ds_upper, pn_obs,  pn_lower, pn_upper, ps_obs,  ps_lower, ps_upper, pn_ps_obs,  pn_ps_lower, pn_ps_upper]
    f = open(sel_out_file, "w")
    f.write('\t'.join(sel_out_head))
    f.write('\n')
    f.write('\t'.join(str(d) for d in sel_out_vals))
    f.close()
