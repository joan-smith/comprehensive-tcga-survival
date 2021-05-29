#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 29 12:21:34 2020

@author: Joan Smith
"""
import scanpy as sc
import pandas as pd
import os

from matplotlib import pyplot as plt
import numpy as np
import seaborn as sns
import glob

sns.set_style({'font_family': 'sans-serif', 'font.sans-serif': 'Arial'})
sns.set(font_scale=1.25)


#%%
dropbox_dir = '~/comprehensive-tcga-survival'
sc.settings.fidir = '/home/joan/code/comprehensive-tcga-survival/figs_and_figdata'

#%% per-platform correlation heatmaps with clustering

platforms = ['rppa', 'rnaseq', 'methylation', 'mutations-non-synonymous', 'cn', 'mirna']
#%%
for p in platforms:
  print(p)
  infile = glob.glob(os.path.join(dropbox_dir, 'Heatmaps', p + '*pancan_correlations_pearson.csv'))[0]
  pearson_corr = pd.read_csv(infile, index_col=0)
  empty_columns =  list(pearson_corr[pearson_corr.isna().sum() == pearson_corr.shape[0]].index.values)
  pearson_corr = pearson_corr.drop(['stouffer unweighted'] + empty_columns, axis=0).drop(['stouffer unweighted'] + empty_columns, axis=1)
  pearson_corr = pearson_corr.dropna(how='any', axis=1)
  ax = sns.clustermap(pearson_corr, cmap='magma', vmin=-0.5, vmax=0.5, cbar_pos=None, xticklabels=True, yticklabels=True)
  ax.savefig(os.path.join(dropbox_dir, 'Heatmaps', p + '_clustermap.pdf'))

#%%
a = np.array([[-0.5,0.5]])
plt.figure(figsize=(1.5, 9))
img = plt.imshow(a, cmap='magma')
plt.gca().set_visible(False)
cax = plt.axes([0.2, 0.6, 0.1, 0.2])
plt.colorbar(orientation="vertical", cax=cax, ticks=[-0.5, -0.25, 0, 0.25, 0.5])
plt.savefig(os.path.join(dropbox_dir, 'Heatmaps', "colorbar.pdf"))

#%%
df = pd.read_csv(os.path.join(dropbox_dir, 'combo_pancan.csv'))

sns.set(font_scale=1.4)

for i, p in enumerate(platforms):
  print(p)
  p_df = df[df['platform'] == p]
  p_df = p_df.set_index('gene').drop(['stouffer unweighted', 'platform'], axis=1)
  empty_cols = list(p_df.loc[:,p_df.isna().sum() == p_df.shape[0]].columns.values)
  do_clustering=True
  if p == 'mutations-non-synonymous':
    p_df = p_df.drop(empty_cols, axis=1).T
    p_df = p_df[p_df.isna().sum().sort_values().head(500).index.values]
    do_clustering=False

  else:
    p_df = p_df.drop(empty_cols, axis=1).dropna(how='any').T
    n = 500
    if p_df.shape[1] < 500:
      n = p_df.shape[1]
    p_df = p_df.sample(n=n, axis=1)
  p_df.columns.name = ''
  show_yticks = False
  if p == 'cn':
    show_yticks = True
  print(p_df.shape)
  ax = sns.clustermap(p_df, cmap='bwr', row_cluster=False, col_cluster=do_clustering, cbar_pos=None, vmin=-4, vmax=4,
                      xticklabels=False, yticklabels=show_yticks)
  if p == 'cn':
    ax.ax_heatmap.tick_params(right=False, left=True, labelleft=True, labelright=False)
  ax.ax_col_dendrogram.set_visible(False)
  ax.savefig(os.path.join(dropbox_dir, 'Heatmaps', p + '_zscore_clustermap_bwr.tiff'))

#%%
a = np.array([[-4,4]])
plt.figure(figsize=(1.5, 9))
img = plt.imshow(a, cmap='bwr')
plt.gca().set_visible(False)
cax = plt.axes([0.2, 0.6, 0.1, 0.2])
plt.colorbar(orientation="vertical", cax=cax, ticks=[-4, -2, 0, 2, 4])
plt.savefig(os.path.join(dropbox_dir, 'Heatmaps', "zscore_colorbar_bwr.tiff"))
