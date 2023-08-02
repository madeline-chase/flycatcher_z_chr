import pandas as pd
import random
import numpy as np
import argparse
from sys import argv

parser = argparse.ArgumentParser(description = 'Estimate pZ/pA and resample sites for CIs')
parser.add_argument('-a', '--count_a', dest = 'count_a', help = 'Fasta formatted file containing coding sequences. CDS must be correct orientation.')
parser.add_argument('-z', '--count_z', dest = 'count_z', help = 'Fasta formatted file containing coding sequences. CDS must be correct orientation.')
parser.add_argument('-tA', '--tot_a', dest = 'tot_a', type = int)
parser.add_argument('-tZ', '--tot_z', dest = 'tot_z', type = int)
parser.add_argument('-n', dest='n_iter', type = int)
args = parser.parse_args()

count_a_f = args.count_a
count_z_f = args.count_z
tot_a = args.tot_a
tot_z = args.tot_z
n_iter = args.n_iter

def calc_pi(p_count,q_count,chr_n,tot):
    p = np.divide(p_count,chr_n)
    q = np.divide(q_count, chr_n)
    pi_snp = np.multiply(p,q)*2
    pi = np.sum(pi_snp)/tot
    return pi

count_a = pd.read_csv(count_a_f, names = ['scaff','start','stop','n_alleles','n_chr','allele1','count1','allele2','count2'], sep = '\t')
count_z = pd.read_csv(count_z_f, names = ['scaff','start','stop','n_alleles','n_chr','allele1','count1','allele2','count2'], sep = '\t')

pi_a_obs = calc_pi(count_a['count1'], count_a['count2'], count_a['n_chr'], tot_a) 
pi_z_obs = calc_pi(count_z['count1'], count_z['count2'], count_z['n_chr'], tot_z)
pi_z_a_obs = pi_z_obs/pi_a_obs

a_tot_snps = len(count_a)
z_tot_snps = len(count_z)

a_index = count_a.index
z_index = count_z.index

pi_a_resamp_list = []
pi_z_resamp_list = []
pi_z_a_resamp_list = []

for iter in range(0, n_iter):
    a_index_res = random.choices(a_index, k = a_tot_snps)
    z_index_res = random.choices(z_index, k = z_tot_snps)

    count_a_res = count_a.loc[a_index_res]
    count_z_res = count_z.loc[z_index_res]

    a_tot_res = tot_a - (a_tot_snps - len(count_a_res))
    z_tot_res = tot_z - (z_tot_snps - len(count_z_res))

    pi_a_res = calc_pi(count_a_res['count1'], count_a_res['count2'], count_a_res['n_chr'],a_tot_res)
    pi_z_res = calc_pi(count_z_res['count1'], count_z_res['count2'], count_z_res['n_chr'], z_tot_res)

    pi_a_resamp_list.append(pi_a_res)
    pi_z_resamp_list.append(pi_z_res)
    pi_z_a_resamp_list.append(pi_z_res/pi_a_res)

pi_a_resamp_list.sort()
pi_z_resamp_list.sort()
pi_z_a_resamp_list.sort()

lower_cut = round(len(pi_a_resamp_list)*0.025)-1
upper_cut = round(len(pi_a_resamp_list)*0.975)-1

pi_a_lower = pi_a_resamp_list[lower_cut]
pi_a_upper = pi_a_resamp_list[upper_cut]

pi_z_lower = pi_z_resamp_list[lower_cut]
pi_z_upper = pi_z_resamp_list[upper_cut]

pi_z_a_lower = pi_z_a_resamp_list[lower_cut]
pi_z_a_upper = pi_z_a_resamp_list[upper_cut]

print('piA', 'piA_lower' , 'piA_upper', 'piZ', 'piZ_lower', 'piZ_upper','piZ/piA', 'piZ/piA_lower','piZ/piA_upper', sep = '\t')
print(pi_a_obs, pi_a_lower, pi_a_upper, pi_z_obs, pi_z_lower, pi_z_upper, pi_z_a_obs, pi_z_a_lower, pi_z_a_upper, sep = '\t')
