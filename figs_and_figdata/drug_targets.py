#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 26 10:12:04 2021

@author: Joan Smith
"""

import pandas as pd
import os
import glob
import sys
import numpy as np
#%%

dropbox_dir = '~/Dropbox/comprehensive-tcga-survival/'
outdir = dropbox_dir + 'Drug targets/'
platforms = {'rppa': dropbox_dir + 'raw-data/TCGA-RPPA-pancan-clean.txt',
             'cn': dropbox_dir + 'cn/cn_by_gene.csv',
             #'mirna': dropbox_dir + 'comprehensive-tcga-survival/raw-data/pancanMiRs_EBadjOnProtocolPlatformWithoutRepsWithUnCorrectMiRs_08_04_16.csv',
             'rnaseq': dropbox_dir + 'raw-data/EBPlusPlusAdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.tsv',
             'methylation': (dropbox_dir + 'comprehensive-tcga-survival/raw-data/jhu-usc.edu_PANCAN_merged_HumanMethylation27_HumanMethylation450.betaValue_whitelisted.tsv',
                            dropbox_dir + 'comprehensive-tcga-survival/raw-data/HumanMethylation450_15017482_v1-2.csv'),
             'mutations-non-synonymous': dropbox_dir + 'raw-data/mc3.v0.2.8.PUBLIC.maf'}
clinical = dropbox_dir + 'raw-data/TCGA-CDR-SupplementalTableS1-2019-05-27.xlsx'

#%% table of cancer types and drug targets from Drug targets.xlsx

drug_targets = pd.read_excel(dropbox_dir + 'Drug targets/Drug targets - post 2018.xlsx', engine='openpyxl', index_col=0)
#drug_targets = pd.read_excel(dropbox_dir + 'Drug targets/Drug targets.xlsx', engine='openpyxl', index_col=0)

drug_targets = drug_targets.dropna(how='any', subset=['Main target'])
cancer_types = drug_targets['2018 or later?'].str.split(',', expand=True)
cancer_types = cancer_types.applymap(lambda x: x.strip() if isinstance(x, str) else x)

df_or = lambda x, y: x | y

combined = pd.get_dummies(cancer_types[0]).astype(int)
for i in cancer_types.columns:
  combined = combined.combine(pd.get_dummies(cancer_types[i]), df_or, fill_value=0, overwrite=False).astype(int)
#%%

main_targets = drug_targets['Main target'].str.split(',', expand=True)
main_targets = main_targets.applymap(lambda x: x.strip() if isinstance(x, str) else x)
main_targets = main_targets.melt(ignore_index=False).dropna()

#%%

targets_and_types_with_drugs = main_targets.join(combined)
targets_and_types = targets_and_types_with_drugs.drop('variable', axis=1)
targets_and_types = targets_and_types.replace(0, np.nan)

def collapse(df):
  gene_row = (df.sum() >= 1).astype(int)
  gene_row.name= df.name
  gene_row = gene_row.replace(0, np.nan)
  gene_row = gene_row.replace(1, gene_row.name)
  return gene_row

targets_and_types = targets_and_types.groupby('value').apply(collapse)
#targets_and_types.to_csv(dropbox_dir + 'Drug targets/targets_and_types.csv')
targets_and_types.to_csv(dropbox_dir + 'Drug targets/targets_and_types-post-2018.csv')


#%%

def select_genes(drugs, pancan_zscores=pd.DataFrame, targets_and_types_with_drugs=pd.DataFrame):
  genes = '\'' + targets_and_types_with_drugs.loc[drugs.index.values, 'value']
  present_targets = set(list(genes.values)) & set(pancan_zscores.index.values)
  if len(present_targets) > 0 and drugs.name in pancan_zscores.columns:
    return pancan_zscores.loc[present_targets, drugs.name]

drugs_and_types = cancer_types.melt(ignore_index=False)

for p in platforms.keys():
   pancan_zscores = pd.read_csv(glob.glob(os.path.join(dropbox_dir, p, '*pancan.csv'))[0], index_col=0)
   ctype_targets = drugs_and_types.groupby('value').apply(select_genes, pancan_zscores=pancan_zscores, targets_and_types_with_drugs=targets_and_types_with_drugs)
   print(ctype_targets)
   ctype_targets = ctype_targets.reset_index().pivot(columns='gene', index='value').T
   ctype_targets.to_csv(outdir + 'cancer_type_specific_zscores_post_2018/' + p + '_zscores.csv')
   ctype_targets.columns.name = 'cancer type'
   ctype_targets.melt(ignore_index=False).dropna().to_csv(outdir + 'cancer_type_specific_zscores_post_2018/' + p + '_melted_zscores.csv')

#%% drug target permutations

sample_counts =  targets_and_types.count()
for p in  platforms.keys():
  pancan_zscores = pd.read_csv(glob.glob(os.path.join(dropbox_dir, p, '*pancan.csv'))[0], index_col=0).drop('stouffer unweighted', axis=1)
  perms = []
  for i in range(1000):
    perm = pancan_zscores.apply(lambda x: x.sample(n=sample_counts[x.name]) if x.name in sample_counts else pd.Series(dtype=np.float64), axis=0)
    perms.extend(perm.melt().dropna().value.values)
  pd.DataFrame(perms).to_csv(outdir + 'permutations/' + p + '1000_perms.csv')





