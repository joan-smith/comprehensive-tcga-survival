#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 22 22:58:07 2020

@author: Joan Smith
"""

import os

import pandas as pd
import anndata as ad

import biomarker_survival as surv

from comprehensive_tcga_survival import rnaseq
from comprehensive_tcga_survival import rppa
from comprehensive_tcga_survival import cn
from comprehensive_tcga_survival import mutations
from comprehensive_tcga_survival import methylation
from comprehensive_tcga_survival import mirna


#%%
dropbox_dir = '~/Dropbox/comprehensive-tcga-survival/'

platforms = {'rppa': dropbox_dir + 'raw-data/TCGA-RPPA-pancan-clean.txt',
             'cn': dropbox_dir + 'cn/cn_by_gene.csv',
             'mirna': dropbox_dir + 'raw-data/pancanMiRs_EBadjOnProtocolPlatformWithoutRepsWithUnCorrectMiRs_08_04_16.csv',
             'rnaseq': dropbox_dir + 'raw-data/EBPlusPlusAdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.tsv',
             'methylation-large': (dropbox_dir + 'raw-data/jhu-usc.edu_PANCAN_HumanMethylation450.betaValue_whitelisted.tsv',
                                   dropbox_dir + 'raw-data/HumanMethylation450_15017482_v1-2.csv'),
             'methylation': (dropbox_dir + 'raw-data/jhu-usc.edu_PANCAN_merged_HumanMethylation27_HumanMethylation450.betaValue_whitelisted.tsv',
                             dropbox_dir +'raw-data/HumanMethylation450_15017482_v1-2.csv'),
             'mutations': dropbox_dir + 'raw-data/mc3.v0.2.8.PUBLIC.maf'}
clinical = dropbox_dir + 'raw-data/TCGA-CDR-SupplementalTableS1-2019-05-27.xlsx'
#%%
platforms_modules = {'rnaseq': rnaseq, 'rppa': rppa, 'cn': cn, 'mutations': mutations, 'methylation': methylation, 'mirna': mirna}

#%%
rppa_key = pd.read_excel(os.path.join(dropbox_dir, 'RPPA conversion.xlsx'), engine='openpyxl')
for c in rppa_key:
  rppa_key[c] = rppa_key[c].str.upper()
rppa_key = rppa_key.set_index('RPPA')

#%%

platforms_list = ['rppa', 'methylation', 'rnaseq', 'cn', 'mutations', 'mirna']


dfs = {}
for i, p in enumerate(platforms_list):
  print(p)
  prep_data = platforms_modules[p].prep_data
  if p == 'cn':
    prep_data = platforms_modules[p].prep_processed_data
  if p == 'methylation':
    zsc = platforms_modules[p].ZscoreCommon(prep_data, platforms_modules[p].ctype_cleaning, extra_data=platforms[p][1])
    data = zsc.metadata(platforms[p][0], clinical, return_data=True)
  else:
    zsc = platforms_modules[p].ZscoreCommon(prep_data, platforms_modules[p].ctype_cleaning)
    data = zsc.metadata(platforms[p], clinical, return_data=True)

  data.columns = [i[1:] if '\'' in i else i for i in data.columns]
  data = data.reset_index()
  data['level_1'] = ad.utils.make_index_unique(data['level_1'])
  data = data.set_index(['level_0', 'level_1', 'time', 'censor'])
  if p == 'rppa':
    rppa_genes = set(list(data.columns)) & set(list(rppa_key.index.values))
    relevant_rppa_key = rppa_key.loc[rppa_genes]
    protein_mapped_rppa = data.loc[:,rppa_genes]
    protein_mapped_rppa.columns = list(relevant_rppa_key['RNA-Seq'].values)
    data.loc[:,rppa_genes] = protein_mapped_rppa
    print(data)

  dfs[p] = data


#%% KMPlots for all genes and platforms
gene_list = list(dfs['cn'].columns)

km = surv.KMPlot()

for g in gene_list:
  one_gene_dfs = {}
  for p in platforms_list:
    if g in dfs[p].columns:
      one_gene_dfs[p] = dfs[p].loc[:,[g]]
      one_gene_dfs[p].columns = one_gene_dfs[p].columns + '_' + p
      g_df = dfs[p].loc[:,[g]].reset_index()
      for c, ctype_g in g_df.groupby('level_0'):
        json = km.km_plot_data(g, ctype_g['time'] , ctype_g['censor'], ctype_g[g])
        print(json)
    else:
      one_gene_dfs[p] = pd.DataFrame(index=['level_0', 'level_1', 'time', 'censor'])

  one_gene_dfs = {k:v for k,v in one_gene_dfs.items() if not v.empty}
  df = pd.concat(one_gene_dfs, axis=1)
  df.to_csv(os.path.join(dropbox_dir, 'KM plot data', 'all_genes', g + '.csv'))

