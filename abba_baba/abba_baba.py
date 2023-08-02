import pandas as pd
import numpy as np
from sys import argv
import argparse

parser = argparse.ArgumentParser(description = 'Perform ABBA/BABA test and jackknife resampling')
parser.add_argument('-p1', dest = 'pop1', help = 'Bed formatted file with allele counts for pop1')
parser.add_argument('-p2', dest = 'pop2', help = 'Bed formatted file with allele counts for pop2')
parser.add_argument('-p3', dest = 'pop3', help = 'Bed formatted file with allele counts for pop3')
parser.add_argument('-w', dest = 'wins', help = 'Bed formatted file with window coordinates for resampling')
args = parser.parse_args()

p1_file = args.pop1
p2_file = args.pop2
p3_file = args.pop3
win_file = args.wins

pop1 = pd.read_csv(p1_file, names = ['Scaff','Start','Stop', 'Alleles','Chr', 'Allele1','Count1','Allele2','Count2'], sep = '\t')
pop2 = pd.read_csv(p2_file, names = ['Scaff','Start','Stop', 'Alleles','Chr', 'Allele1','Count1','Allele2','Count2'], sep = '\t')
pop3 = pd.read_csv(p3_file, names = ['Scaff','Start','Stop', 'Alleles','Chr', 'Allele1','Count1','Allele2','Count2'], sep = '\t')

def abba(p1,p2,p3):
    f = (1-p1)*p2*p3
    return f

def baba(p1,p2,p3):
    f = p1*(1-p2)*p3
    return f

abba_sites = abba((pop1['Count2']/pop1['Chr']),(pop2['Count2']/pop2['Chr']),(pop3['Count2']/pop3['Chr']))

baba_sites = baba((pop1['Count2']/pop1['Chr']),(pop2['Count2']/pop2['Chr']),(pop3['Count2']/pop3['Chr']))

d_obs = (np.sum(abba_sites)-np.sum(baba_sites))/(np.sum(abba_sites)+np.sum(baba_sites))

auto_wins = pd.read_csv(win_file, sep = '\t', names = ['scaff','start','stop'])

d_stats_jack = []

for win in range(0,len(auto_wins)):
    scaff = auto_wins.loc[win, 'scaff']
    start = auto_wins.loc[win, 'start']
    stop = auto_wins.loc[win, 'stop']
    
    pop1_jack = pop1[~((pop1['Scaff']==scaff) & (pop1['Stop'] >= start) & (pop1['Stop'] < stop))]
    pop2_jack = pop2[~((pop2['Scaff']==scaff) & (pop2['Stop'] >= start) & (pop2['Stop'] < stop))]
    pop3_jack = pop3[~((pop3['Scaff']==scaff) & (pop3['Stop'] >= start) & (pop3['Stop'] < stop))]
    
    abba_sites_jack = abba((pop1_jack['Count2']/pop1_jack['Chr']),(pop2_jack['Count2']/pop2_jack['Chr']),(pop3_jack['Count2']/pop3_jack['Chr']))
    baba_sites_jack = baba((pop1_jack['Count2']/pop1_jack['Chr']),(pop2_jack['Count2']/pop2_jack['Chr']),(pop3_jack['Count2']/pop3_jack['Chr']))
    
    d_jack = (np.sum(abba_sites_jack)-np.sum(baba_sites_jack))/(np.sum(abba_sites_jack)+np.sum(baba_sites_jack))
    
    d_stats_jack.append(d_jack)


jack_mean = np.mean(d_stats_jack)
jack_se =  np.sqrt(((len(auto_wins)-1)/len(auto_wins)) * np.sum((d_stats_jack - jack_mean)**2))
jack_z = d_obs/jack_se

print(d_obs, jack_mean, jack_se, jack_z, sep = '\t')