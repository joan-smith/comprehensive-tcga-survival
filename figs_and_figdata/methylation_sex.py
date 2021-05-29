#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 24 16:48:55 2020

@author: Joan Smith
"""

import pandas as pd
import os

import scipy.stats as stats
import biomarker_survival as surv

from comprehensive_tcga_survival import methylation
#%%
dropbox_dir = '~/Dropbox/'
methylation_data = dropbox_dir + 'comprehensive-tcga-survival/raw-data/jhu-usc.edu_PANCAN_merged_HumanMethylation27_HumanMethylation450.betaValue_whitelisted.tsv'
methylation_key = dropbox_dir + 'comprehensive-tcga-survival/raw-data/HumanMethylation450_15017482_v1-2.csv'
clinical = dropbox_dir + 'comprehensive-tcga-survival/raw-data/TCGA-CDR-SupplementalTableS1-2019-05-27.xlsx'

#%%
df = methylation.prep_data(methylation_data, extra_data=methylation_key)
tcga_cdr = surv.TCGA_CDR_util(clinical)
c_df = tcga_cdr.cancer_type_data('*', extra_cols=['gender'])[['gender']]


dfs = [methylation.ctype_cleaning(df, ctype, None) for ctype in tcga_cdr.cancer_types()]
dfs = pd.concat(dfs)
joined = dfs.join(c_df)

def ttest(j):
  site = j.name
  j = j.reset_index().dropna(how='any')
  m = j[j['gender'] == 'MALE'][site].array
  f = j[j['gender'] == 'FEMALE'][site].array
  t, p = stats.ttest_ind(m, f)
  return pd.Series({'p': p, 't': t}, name=site)

joined = joined.reset_index().set_index(['index', 'gender'])

results = joined.apply(ttest)
results.T.to_csv(os.path.join(dropbox_dir, 'comprehensive-tcga-survival', 'figs', 'methylation_sex_ttest.csv'))
