#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan  1 20:16:22 2021

@author: Joan Smith
"""

#%%
import pandas as pd
import os
import scipy
import numpy as np

from matplotlib import pyplot as plt

import matplotlib
import scipy.signal

matplotlib.use('Qt5Agg')
#%%

dropbox_dir = '~/Dropbox'

cn_in =  dropbox_dir + '/cn/cn_by_gene.csv'
#%%

cn = pd.read_csv(cn_in, index_col=0, usecols=[0, 21853, 21854])
cn.index = '\'' + cn.index

#%%

pancan_all = pd.read_csv(os.path.join(dropbox_dir, 'cn', 'cn-pancan.csv'), index_col=0)
pancan_all = pancan_all.join(cn).sort_values(['chr', 'start']).drop(['chr', 'start'], axis=1)

#%%
pancan = pancan_all[['stouffer unweighted']].join(cn).sort_values(['chr', 'start'])

#%%

plt.figure()

peaks, properties = scipy.signal.find_peaks(pancan['stouffer unweighted'], distance=500,
                                            width=4, height=3)

low_peaks, low_properties = scipy.signal.find_peaks(pancan['stouffer unweighted']*-1,
                                                     distance=500, width=4, height=3)

plt.scatter(low_peaks, low_properties['peak_heights']*-1, marker='x', c='blue')
plt.scatter(peaks, properties['peak_heights'], marker='x', c='red')

plt.scatter(range(pancan.shape[0]), pancan['stouffer unweighted'], s=.5, c='gray')

widths = scipy.signal.peak_widths(pancan['stouffer unweighted'], peaks, rel_height=.20)
low_widths = scipy.signal.peak_widths(pancan['stouffer unweighted']*-1, low_peaks, rel_height=.20)


pancan_gene_list_peaks = []
peak_indicies = [(int(np.floor(i[0])), int(np.ceil(i[1]))) for i in zip(widths[2], widths[3])]
for j,i in enumerate(peak_indicies):
    pancan_gene_list_peaks.extend(pancan.iloc[i[0]:i[1]].index.values)
    xs = [pancan.index.get_loc(x) for x in pancan.iloc[i[0]:i[1]].index]
    plt.scatter(xs, pancan.iloc[i[0]:i[1]]['stouffer unweighted'], s=.5, c='red')
pd.Series(pancan_gene_list_peaks).to_csv(dropbox_dir + '/Chromosome coordinates/peak details.csv')

pancan_gene_list_valleys = []
low_peak_indicies = [(int(np.floor(i[0])), int(np.ceil(i[1]))) for i in zip(low_widths[2], low_widths[3])]
for j,i in enumerate(low_peak_indicies):
    pancan_gene_list_valleys.extend(pancan.iloc[i[0]:i[1]].index.values)
    xs = [pancan.index.get_loc(x) for x in pancan.iloc[i[0]:i[1]].index]
    plt.scatter(xs, pancan.iloc[i[0]:i[1]]['stouffer unweighted'], s=.5, c='blue')
pd.Series(pancan_gene_list_valleys).to_csv(dropbox_dir + '/Chromosome coordinates/valley details.csv')

plt.show()
plt.savefig(dropbox_dir + '/Chromosome coordinates/colored peaks.png')


#%% per-cancer type peaks

peak_writer = pd.ExcelWriter(os.path.join(dropbox_dir, 'Chromosome coordinates', 'peaks_and_valleys_gene_list.xlsx'))

for c in pancan_all:
  print(c)
  peaks, properties = scipy.signal.find_peaks(pancan_all[c], distance=500,
                                            width=4, height=3)
  low_peaks, low_properties = scipy.signal.find_peaks(pancan_all[c]*-1,
                                                     distance=500, width=4, height=3)

  widths = scipy.signal.peak_widths(pancan_all[c], peaks, rel_height=.20)
  low_widths = scipy.signal.peak_widths(pancan_all[c]*-1, low_peaks, rel_height=.20)

  gene_list_peaks = []
  peak_indicies = [(int(np.ceil(i[0])), int(np.ceil(i[1]))) for i in zip(widths[2], widths[3])]
  ps = {}

  for j,i in enumerate(peak_indicies):
    gs = pancan_all.iloc[i[0]:i[1]].index.values
    p_df = pd.DataFrame(gs, dtype=str, index=gs).join(cn).sort_values(['chr', 'start'])
    start = str(int(p_df['start'].min()))
    end = str(int(p_df['start'].max()))
    assert len(p_df['chr'].value_counts()) == 1
    chrm = str(int(p_df['chr'].value_counts().index[0]))
    gene_names = p_df[0].values
    ps['chr ' + chrm + ': ' + start + '-' + end] = pd.Series(gene_names)
  pd.DataFrame(ps).to_excel(peak_writer, sheet_name=c + ' - peaks', index=False)

  gene_list_valleys = []
  low_peak_indicies = [(int(np.ceil(i[0])), int(np.ceil(i[1]))) for i in zip(low_widths[2], low_widths[3])]
  vs = {}
  for j,i in enumerate(low_peak_indicies):
    gs = pancan_all.iloc[i[0]:i[1]].index.values
    v_df = pd.DataFrame(gs, dtype=str, index=gs).join(cn).sort_values(['chr', 'start'])
    start = str(int(v_df['start'].min()))
    end = str(int(v_df['start'].max()))
    assert len(v_df['chr'].value_counts()) == 1
    chrm = str(int(v_df['chr'].value_counts().index[0]))
    gene_names = v_df[0].values
    vs['chr ' + chrm + ': ' + start + '-' + end] = pd.Series(gene_names)
  pd.DataFrame(vs).to_excel(peak_writer, sheet_name=c + ' - valleys', index=False)

peak_writer.close()


