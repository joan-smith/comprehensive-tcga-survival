#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 25 21:28:40 2020

@author: Joan Smith
"""
import pandas as pd
import os
from scipy.stats import percentileofscore
import glob

import biomarker_survival as surv

from comprehensive_tcga_survival import mutations
#%%

dropbox_dir = '~/Dropbox/comprehensive-tcga-survival/'
outdir = dropbox_dir + 'Driver_genes/'
platforms = {'rppa': dropbox_dir + 'raw-data/TCGA-RPPA-pancan-clean.txt',
             'cn': dropbox_dir + 'cn/cn_by_gene.csv',
             #'mirna': dropbox_dir + 'comprehensive-tcga-survival/raw-data/pancanMiRs_EBadjOnProtocolPlatformWithoutRepsWithUnCorrectMiRs_08_04_16.csv',
             'rnaseq': dropbox_dir + 'raw-data/EBPlusPlusAdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.tsv',
             'methylation': (dropbox_dir + 'comprehensive-tcga-survival/raw-data/jhu-usc.edu_PANCAN_merged_HumanMethylation27_HumanMethylation450.betaValue_whitelisted.tsv',
                            dropbox_dir + 'comprehensive-tcga-survival/raw-data/HumanMethylation450_15017482_v1-2.csv'),
             'mutations-non-synonymous': dropbox_dir + 'raw-data/mc3.v0.2.8.PUBLIC.maf'}
clinical = dropbox_dir + 'raw-data/TCGA-CDR-SupplementalTableS1-2019-05-27.xlsx'


#%% Driver gene zscores, pancan driver list, individual cytpes

driver_genes = pd.read_excel(os.path.join(dropbox_dir, 'Driver_genes', 'oncogenes pancan.xlsx'), engine='openpyxl')
pancan_drivers = '\'' + driver_genes['Gene']

for p in platforms.keys():
    pancan_zscores = pd.read_csv(glob.glob(os.path.join(dropbox_dir, p, '*pancan.csv'))[0], index_col=0)
    present_drivers = list(set(pancan_drivers) & set(pancan_zscores.index.values))
    pancan_zscores.loc[present_drivers].to_csv(outdir + 'pancan_driver_gene_zscores/' + p + '_zscores.csv')

#%% pancan driver genes permutations, mutations

num_genes = len(pancan_drivers)
mutation_zscores = pd.read_csv(glob.glob(os.path.join(dropbox_dir, 'mutations-non-synonymous', '*pancan.csv'))[0], index_col=0)
mutation_zscores = mutation_zscores.drop('stouffer unweighted', axis=1)

vals = []
sig_ctypes = []
for i in range(1000):
  perm = mutation_zscores.sample(num_genes)
  vals.extend(perm.melt().dropna()['value'].values)
  significant_ctypes = (perm > 1.96).sum(axis=1).values
  sig_ctypes.extend(significant_ctypes)

#pd.DataFrame(vals).to_csv(outdir + 'mutations_1000_permutations.csv', index_label=None, header=None)
pd.DataFrame(sig_ctypes).to_csv(outdir + 'mutations_1000_permutations_significant_ctype_count.csv', index_label=None, header=None)


#%%% pancan driver genes permutations, rnaseq

num_genes = len(pancan_drivers)
rnaseq_zscores = pd.read_csv(glob.glob(os.path.join(dropbox_dir, 'rnaseq', '*pancan.csv'))[0], index_col=0)
rnaseq_zscores = rnaseq_zscores.drop('stouffer unweighted', axis=1)

vals = []
sig_ctypes = []
for i in range(1000):
  perm = mutation_zscores.sample(num_genes)
  vals.extend(perm.melt().dropna()['value'].values)
  significant_ctypes = (perm > 1.96).sum(axis=1).values
  sig_ctypes.extend(significant_ctypes)

#pd.DataFrame(vals).to_csv(outdir + 'rnaseq_1000_permutations.csv', index_label=None, header=None)
pd.DataFrame(sig_ctypes).to_csv(outdir + 'rnaseq_1000_permutations_significant_ctype_count.csv', index_label=None, header=None)


#%% Driver gene zscores per-cancer-type list

driver_genes =  pd.read_excel(os.path.join(dropbox_dir, 'Driver_genes', 'oncogenes - cancer specific.xlsx'), engine='openpyxl')
driver_genes['Gene'] = '\'' + driver_genes['Gene']

def select_driver_genes(dg, pancan_zscores=pd.DataFrame):
  present_drivers = list(set(dg['Gene']) & set(pancan_zscores.index.values))
  pancan_zscores.loc
  if dg.name in pancan_zscores:
    return pancan_zscores.loc[present_drivers, dg.name]

for p in platforms.keys():
   pancan_zscores = pd.read_csv(glob.glob(os.path.join(dropbox_dir, p, '*pancan.csv'))[0], index_col=0)
   ctype_drivers = driver_genes.groupby('Cancer').apply(select_driver_genes, pancan_zscores=pancan_zscores)
   ctype_drivers.to_csv(outdir + 'cancer_type_specific_driver_gene_zscores/' + p + '_zscores.csv')


#%% Pancan driver genes

driver_genes = pd.read_excel(os.path.join(dropbox_dir, 'Driver_genes', 'oncogenes.xlsx'), engine='openpyxl')
pancan_drivers = '\'' + driver_genes[driver_genes['Cancer'] == 'PANCAN']['Gene']


for p in platforms.keys():
  ranks = pd.DataFrame(index=pancan_drivers)
  pancan_zscores = pd.read_csv(glob.glob(os.path.join(dropbox_dir, p, '*pancan.csv'))[0], index_col=0)
  print(pancan_zscores)
  present_driver_genes = list(set(pancan_drivers) & set(pancan_zscores.index.values))
  print(present_driver_genes)

  all_type_ranks = {}
  for c in pancan_zscores:
    ctype_ranks = {}
    for gene in present_driver_genes:
      rank = percentileofscore(pancan_zscores[c].dropna().values, pancan_zscores.loc[gene, c], kind='weak')
      ctype_ranks[gene] = rank
      print(p, c, gene, rank)
      ranks.loc[gene, c] = rank
  ranks.to_csv(os.path.join(dropbox_dir, 'Driver_genes', p + '_pancan_driver_genes.csv'))

#%% Cancer type specific driver genes

driver_genes = pd.read_excel(os.path.join(dropbox_dir, 'Driver_genes', 'oncogenes.xlsx'), engine='openpyxl')

for p in platforms.keys():
  ranks = pd.DataFrame(index=pancan_drivers)
  pancan_zscores = pd.read_csv(glob.glob(os.path.join(dropbox_dir, p, '*pancan.csv'))[0], index_col=0)

  all_type_ranks = {}
  for c in pancan_zscores:
    ctype_ranks = {}
    present_driver_genes = list(set('\'' + driver_genes[driver_genes['Cancer'] == c]['Gene']) & set(pancan_zscores[c].dropna().index))
    print(present_driver_genes)
    for gene in present_driver_genes:
      rank = percentileofscore(pancan_zscores[c].dropna().values, pancan_zscores.loc[gene, c], kind='weak')
      ctype_ranks[gene] = rank
      print(p, c, gene, rank)
      ranks.loc[gene, c] = rank
  ranks.to_csv(os.path.join(dropbox_dir, 'Driver_genes', p + '_ctype_driver_genes.csv'))

#%% Count driver gene mutations

tcga_cdr = surv.TCGA_CDR_util(clinical)
driver_genes = pd.read_excel(os.path.join(dropbox_dir, 'Driver_genes', 'oncogenes and tumor suppressors pancan.xlsx'), engine='openpyxl')
mut_data = mutations.prep_data(platforms['mutations-non-synonymous'])
#%%
outdir = dropbox_dir + 'Driver_genes'
out_df = {}
for c in tcga_cdr.cancer_types():
  mut_ctype = mutations.ctype_cleaning(mut_data, c, tcga_cdr.cancer_type_data(c))
  ctype_drivers = driver_genes[(driver_genes['Cancer'] == c) | (driver_genes['Cancer'] == 'PANCAN')]
  ctype_drivers.loc[:,'Gene'] = '\'' + ctype_drivers['Gene']
  present_drivers = set(ctype_drivers['Gene'].values) & set(mut_ctype.columns.values)
  ctype_counts = mut_ctype.loc[:,present_drivers].sum()
  ctype_counts.name = c
  ctype_counts.loc['patient count'] = mut_ctype.shape[0]
  out_df[c] = ctype_counts

pd.DataFrame(out_df).to_csv(os.path.join(outdir, 'mutation_driver_gene_counts_oncogenes_and_tumor_suppressors.csv'))


