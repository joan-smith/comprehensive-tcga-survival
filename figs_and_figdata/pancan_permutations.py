#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 30 12:59:06 2020

@author: Joan Smith
"""
import pandas as pd
import os
import numpy as np

import glob
#%%

dropbox_dir = '~/Dropbox/comprehensive-tcga-survival'
permutation_dir = os.path.join(dropbox_dir, 'permutations')
#%%
platforms = ['mirna', 'rnaseq', 'cn', 'mutations-non-synonymous', 'methylation', 'rppa']
#%%
n_permutations = 1000

for p in platforms:
  df = pd.read_csv(glob.glob(os.path.join(dropbox_dir, p, '*pancan.csv'))[0], index_col=0)
  df = df.drop('stouffer unweighted', axis=1)

  permuted_counts = pd.DataFrame({'GT 1.96': [0]*df.shape[1], 'LT -1.96': [0]*df.shape[1]}, index=range(df.shape[1]))

  #start permute
  for _ in range(n_permutations):
    df_p = df.apply(np.random.permutation)
    gt_counts = (df_p > 1.96).sum(axis=1).value_counts().reindex(range(df.shape[1]), fill_value=0)
    lt_counts = (df_p < -1.96).sum(axis=1).value_counts().reindex(range(df.shape[1]), fill_value=0)
    permuted_counts['GT 1.96'] += gt_counts
    permuted_counts['LT -1.96'] += lt_counts

  permuted_counts.to_csv(os.path.join(permutation_dir, p + '_' + str(n_permutations) + '_permutation_counts.csv'))

#%%

data_dir = '~/Dropbox/comprehensive-tcga-survival/age-stage-grade-sex'
#data_dir = '~/Dropbox/comprehensive-tcga-survival'
#data_dir = '~/Dropbox/comprehensive-tcga-survival/p53-mutant-multivariate-zscores'

#%%
dfs = []

for p in platforms:
  df = pd.read_csv(glob.glob(os.path.join(data_dir, p, '*pancan.csv'))[0])
  df['platform'] = p
  dfs.append(df)

combined = pd.concat(dfs, axis=0)
combined = combined.set_index(['platform', 'gene'])
combined.to_csv(os.path.join(data_dir, 'combo_pancan.csv'))

#%% Univaraite Pancan Correlations per platform
combined = pd.read_csv(dropbox_dir + '/combo_pancan.csv')

for p in combined['platform'].unique():
  pdata = combined[combined['platform'] == p]
  print(pdata.shape[0])
  pdata.corr().to_csv(dropbox_dir + '/Heatmaps/' + p +'_pancan_correlations_pearson.csv')

#%%
combined = pd.read_csv(os.path.join(data_dir, 'combo_pancan.csv'), index_col=[0,1])
#%%
def largest(c):
  nlargest = c.abs().nlargest(100)
  df = c[nlargest.index].reset_index()
  a = df.platform + ' ' + df.gene
  return a

def platform_counts(c):
  nlargest = c.abs().nlargest(100)
  df = c[nlargest.index].reset_index()
  return df.value_counts('platform')

#%%
df = combined.apply(largest, result_type='expand', axis=0)
df.to_csv(os.path.join(data_dir, 'top-100.csv'), index_label=None)
#%%
platform_counts = combined.apply(platform_counts, result_type='expand', axis=0)
platform_counts.to_csv(os.path.join(data_dir,  'platform_counts.csv'), index_label=None)