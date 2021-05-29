#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 30 10:53:32 2020

@author: Joan Smith
"""
import pandas as pd
import os
import numpy as np


import seaborn as sns
from matplotlib import pyplot as plt
#%%
dropbox_dir = '~/Dropbox/comprehensive-tcga-survival'
cross_platform_correlations_dir = os.path.join(dropbox_dir, 'cross-platform correlations')

df = pd.read_csv(os.path.join(dropbox_dir, 'combo_pancan.csv'))
rppa_key = pd.read_excel(os.path.join(dropbox_dir, 'RPPA conversion.xlsx'), engine='openpyxl')

#%% Cross platform stouffer's zscore correlations

all_genes = df[df['platform'] == 'cn']['gene'].str.upper()
corr_df = pd.DataFrame(index=all_genes)

platforms = ['rnaseq', 'cn', 'mutations-non-synonymous', 'methylation']
for p in platforms:

  stouffer = df.loc[df['platform'] == p, ['stouffer unweighted', 'gene']]
  stouffer['gene'] = stouffer['gene'].str.upper()
  stouffer = stouffer.set_index('gene')
  corr_df[p] = stouffer

rppa = df.loc[df['platform'] == 'rppa',['stouffer unweighted', 'gene']]
rppa['gene'] = rppa['gene'].str.upper()
rppa = rppa.set_index('gene')

rppa_key = pd.read_excel(os.path.join(cross_platform_correlations_dir, 'RPPA conversion.xlsx'), engine='openpyxl')
for c in rppa_key:
  rppa_key[c] = '\'' + rppa_key[c].str.upper()
rppa_key = rppa_key.set_index('RPPA')

rppa = rppa.join(rppa_key, how='inner').set_index('RNA-Seq')
corr_df['rppa'] = rppa

correlates = corr_df.corr()
correlates.to_csv(os.path.join(cross_platform_correlations_dir, 'cross_platform_correlations.csv'))

#%%
sns.set_style({'font_family': 'sans-serif', 'font.sans-serif': 'Arial'})
sns.set(font_scale=2)

correlates.columns = ['Gene expression', 'Copy number alterations', 'Mutations', 'DNA methylation', 'Protein expression']
correlates.index = ['Gene expression', 'Copy number alterations', 'Mutations', 'DNA methylation', 'Protein expression']
ax = sns.clustermap(correlates, vmin=-1, vmax=1, cmap='twilight', cbar_pos=None, xticklabels=True, yticklabels=True)
ax.savefig(os.path.join(cross_platform_correlations_dir, 'correlates_heatmap.pdf'))

#%%
a = np.array([[-1,1]])
plt.figure(figsize=(1.5, 9))
img = plt.imshow(a, cmap='twilight')
plt.gca().set_visible(False)
cax = plt.axes([0.2, 0.6, 0.1, 0.2])
plt.colorbar(orientation="vertical", cax=cax, ticks=[-1, -0.5, 0, 0.5, 1])
plt.savefig(os.path.join(cross_platform_correlations_dir, "colorbar.pdf"))

#%%
platforms = ['rnaseq', 'cn', 'mutations-non-synonymous', 'methylation', 'rppa']

ctype_corrs = {}
for p in platforms:
  ctype_df = df.loc[df['platform'] == p].drop(['stouffer unweighted', 'platform'], axis=1)
  ctype_df['gene'] = ctype_df['gene'].str.upper()
  ctype_df = ctype_df.set_index('gene')
  if p == 'rppa':
    ctype_df = ctype_df.join(rppa_key, how='inner').set_index('RNA-Seq')
    ctype_df.index.name = 'gene'
    print(ctype_df)
  melt = ctype_df.melt(ignore_index=False).reset_index().set_index(['gene', 'variable'])
  ctype_corrs[p] = melt


#%%
ctype_corrs_df = pd.concat(ctype_corrs, axis=1)
#%%
ctype_correlates = ctype_corrs_df.corr()
ctype_correlates.to_csv(os.path.join(cross_platform_correlations_dir, 'cancer_type_cross_platform_correlates.csv'))

ctype_correlates.columns = ['Gene expression', 'Copy number alterations', 'Mutations', 'DNA methylation', 'Protein expression']
ctype_correlates.index = ['Gene expression', 'Copy number alterations', 'Mutations', 'DNA methylation', 'Protein expression']
ax = sns.clustermap(ctype_correlates, vmin=-1, vmax=1, cmap='twilight', cbar_pos=None, xticklabels=True, yticklabels=True)
ax.savefig(os.path.join(cross_platform_correlations_dir, 'cancer_type_correlates_heatmap.pdf'))

