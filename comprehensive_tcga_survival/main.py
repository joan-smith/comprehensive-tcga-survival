#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Dec 20 14:43:34 2020

@author: Joan Smith
"""
import pandas as pd
import os
import pathlib
import scipy
import glob

import biomarker_survival as surv

from comprehensive_tcga_survival import rppa
from comprehensive_tcga_survival import cn
from comprehensive_tcga_survival import mirna
from comprehensive_tcga_survival import rnaseq
from comprehensive_tcga_survival import methylation
from comprehensive_tcga_survival import mutations

from comprehensive_tcga_survival import Multivariate

from statsmodels.stats.multitest import fdrcorrection


#%%
dropbox_dir = '~/Dropbox/'
stage_key = dropbox_dir + 'comprehensive-tcga-survival/Multivariate/Stage_key.xlsx'
grade_key = dropbox_dir + 'comprehensive-tcga-survival/Multivariate/Grade_key.xlsx'
clinical = dropbox_dir + 'comprehensive-tcga-survival/raw-data/TCGA-CDR-SupplementalTableS1-2019-05-27.xlsx'
raw_data = dropbox_dir + 'comprehensive-tcga-survival/raw-data/'
#%%

platforms = {'rppa': dropbox_dir + 'comprehensive-tcga-survival/raw-data/TCGA-RPPA-pancan-clean.txt',
             'cn': dropbox_dir + 'comprehensive-tcga-survival/cn/cn_by_gene.csv',
             'mirna': dropbox_dir + 'comprehensive-tcga-survival/raw-data/pancanMiRs_EBadjOnProtocolPlatformWithoutRepsWithUnCorrectMiRs_08_04_16.csv',
             'rnaseq': dropbox_dir + 'comprehensive-tcga-survival/raw-data/EBPlusPlusAdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.tsv',
             'methylation': (dropbox_dir + 'comprehensive-tcga-survival/raw-data/jhu-usc.edu_PANCAN_merged_HumanMethylation27_HumanMethylation450.betaValue_whitelisted.tsv',
                             dropbox_dir + 'comprehensive-tcga-survival/raw-data/HumanMethylation450_15017482_v1-2.csv'),
             'mutations-non-synonymous': dropbox_dir + 'comprehensive-tcga-survival/raw-data/mc3.v0.2.8.PUBLIC.maf'}
#%%

tcga_cdr_local = surv.TCGA_CDR_util(clinical)
cancer_types = tcga_cdr_local.cancer_types()

#%% Prep univariate zscores
platform_outdir = dropbox_dir + 'comprehensive-tcga-survival/'
parallel = 35
for p, rd in platforms.items():
  pathlib.Path(platform_outdir).mkdir(parents=True, exist_ok=True)
  if p == 'rppa':
    rppa.zscores(rd, clinical, platform_outdir, parallel)
  if p == 'cn':
    cn.zscores(rd, clinical, platform_outdir, parallel)
  if p == 'rnaseq':
    rnaseq.zscores(rd, clinical, platform_outdir, parallel)
  if p == 'mirna':
    mirna.zscores(rd, clinical, platform_outdir, parallel)
  if p == 'methylation':
    methylation_data = rd[0]
    methylation_key = rd[1]
    methylation.zscores(methylation_data, methylation_key, clinical, platform_outdir, parallel)
  if p == 'mutations':
    mutations.zscores(rd, clinical, platform_outdir, parallel)

#%% Prep Metadata
outdir = dropbox_dir + 'comprehensive-tcga-survival/'
metadata = {}
for p, rd in platforms.items():
  if p == 'rppa':
    metadata['rppa'] = rppa.metadata(rd, clinical)
  if p == 'cn':
    metadata['cn'] = cn.metadata(rd, clinical)
  if p == 'mirna':
    metadata['mirna'] = mirna.metadata(rd, clinical)
  if p == 'rnaseq':
    metadata['rnaseq'] = rnaseq.metadata(rd, clinical)
  if p == 'methylation-large':
    methylation_large_data = rd[0]
    methylation_key = rd[1]
    metadata['methylation-large'] = methylation.metadata(methylation_large_data, methylation_key, clinical)
  if p == 'methylation':
    methylation_data = rd[0]
    methylation_key = rd[1]
    metadata['methylation'] = methylation.metadata(methylation_data, methylation_key, clinical)
  if p == 'mutations':
    metadata['mutations'] = mutations.metadata(rd, clinical)

  pd.DataFrame(metadata).to_csv(os.path.join(outdir, 'patient_counts_cloud.csv'))

#%% Prep Multivariate
parallel = 35
outdir = dropbox_dir + 'comprehensive-tcga-survival/age-stage-grade-sex'

multivar = Multivariate(tcga_cdr_local, stage_key, grade_key)
ctype_multivars = multivar.prep_all_multivar()

for p, rd in platforms.items():
  platform_outdir = os.path.join(outdir, p)
  if p == 'rppa':
    rppa.zscores(rd, clinical, platform_outdir, parallel, additional_vars=ctype_multivars)
  if p == 'cn':
    cn.zscores(rd, clinical, platform_outdir, 14, additional_vars=ctype_multivars)
  if p == 'rnaseq':
    rnaseq.zscores(rd, clinical, platform_outdir, parallel, additional_vars=ctype_multivars)
  if p == 'mirna':
    mirna.zscores(rd, clinical, platform_outdir, parallel, additional_vars=ctype_multivars)
  if p == 'methylation':
    methylation_data = rd[0]
    methylation_key = rd[1]
    methylation.zscores(methylation_data, methylation_key, clinical, platform_outdir, parallel, additional_vars=ctype_multivars)
  if p == 'mutations':
    mutations.zscores(rd, clinical, platform_outdir, parallel, additional_vars=ctype_multivars)

#%% RNASEQ: Sex, RPSY41, and XIST

rnaseq_data = rnaseq.prep_data(platforms['rnaseq'])

dfs = []
for c in cancer_types:
  df = rnaseq.ctype_cleaning(rnaseq_data, c, None)
  df = df.join(tcga_cdr_local.cancer_type_data(c, extra_cols=['gender']), how='inner')
  out_df = df[['XIST', 'RPS4Y1', 'gender']].copy()
  print(c)
  print(out_df)

  out_df['cancer_type'] = c
  dfs.append(out_df)

dfs = pd.concat(dfs, axis=0)
dfs.to_csv(os.path.join(outdir, 'XIST_RPS4Y1.csv'), index_label='patient', columns=['cancer_type', 'gender', 'XIST', 'RPS4Y1'])

#%% Methylation: Sex + methylation sites

methylation_data = methylation.prep_data(platforms['methylation'][0], platforms['methylation'][1])

genes = ['MOSPD1', 'SLC9A7', 'MTMR1', 'APEX2', 'OTUD5','PHF8', 'OCRL', 'FAM50A',
         'IRAK1','FTSJ1', 'PRDX4', 'EMD', 'TMEM185A', 'NDUFB11', 'RBM10','TSPYL2',
         'ELK1', 'HTATSF1', 'BCOR','CLCN5']

dfs = []
for c in cancer_types:
  df = methylation.ctype_cleaning(methylation_data, c, None)
  df = df.join(tcga_cdr_local.cancer_type_data(c, extra_cols=['gender']), how='inner')
  out_df = df[genes + ['gender']].copy()

  out_df['cancer_type'] = c
  dfs.append(out_df)

dfs = pd.concat(dfs, axis=0)
dfs.to_csv(os.path.join(outdir, 'methylation_sex.csv'), index_label='patient')

#%% FDR correction -- group all platforms together for cancer type

indir = dropbox_dir + '/comprehensive-tcga-survival'

def count_sig(g):
  print(g.name)
  print('uncorrected sig:', (g['value'] <= 0.05).sum())
  significant_corrected = g['corrected-p'] <= 0.05
  print(significant_corrected.sum(), sum(significant_corrected))

  return sum(g['corrected-p'] <= 0.05)

significant_corrected_pvals = {}
for c in cancer_types:
  print(c)
  pvals = {}
  for p in platforms.keys():
    platform_dir = os.path.join(indir, p)
    ctype_file = c + '.zscores.out.csv'
    if p == 'mutations':
      ctype_file = c + '_mutation-fraction-0.02.zscores.out.csv'
    df = pd.read_csv(os.path.join(platform_dir, ctype_file), index_col=0)
    if('p' in df.columns):
      pvals[p] = df.p
      pval_df = pd.DataFrame(df.p).melt().dropna()
      rejected, corrected_p = fdrcorrection(pval_df['value'].values)
      pval_df['corrected-p'] = corrected_p

  significant_corrected_pvals[c] = pval_df.groupby('variable').apply(count_sig)

print(pd.DataFrame(significant_corrected_pvals))
outdf = pd.DataFrame(significant_corrected_pvals).T
outdf.to_csv(os.path.join(outdir, 'significant_corrected_pvals.csv'))


#%% FDR correction -- within cancer type+platform fdr

outdir = dropbox_dir + 'comprehensive-tcga-survival/univariate-fdr'
indir = dropbox_dir + 'comprehensive-tcga-survival'

corrected_p_df = {}
for c in cancer_types:
  platforms_corrected_counts = {}
  for p in platforms.keys():
    platform_dir = os.path.join(indir, p)
    ctype_file = c + '.zscores.out.csv'
    if 'mutations' in p:
      ctype_file = c + '_0.02.zscores.out.csv'
    df = pd.read_csv(os.path.join(platform_dir, ctype_file), index_col=0)
    if('p' in df.columns):
      rejected, corrected_p = fdrcorrection(df['p'].dropna(), alpha=0.05)
      platforms_corrected_counts[p] = rejected.sum()

  corrected_p_df[c] = platforms_corrected_counts

outdf = pd.DataFrame(corrected_p_df).T
outdf.to_csv(os.path.join(outdir, 'significant_corrected_pvals_ctype_0.05.csv'))

#%% stouffer fdr

outdir = dropbox_dir + 'comprehensive-tcga-survival/univariate-fdr'

for p in platforms.keys():
  pancan = pd.read_csv(glob.glob(os.path.join(indir, p, '*pancan.csv'))[0], index_col=0)
  pancan = pancan.dropna(subset=['stouffer unweighted'])
  pancan['stouffer-p'] = scipy.stats.norm.sf(abs(pancan['stouffer unweighted'].values))*2
  print(pancan[['stouffer unweighted', 'stouffer-p']])
  rejected05, corrected_p = fdrcorrection(pancan['stouffer-p'], alpha = 0.05)
  rejected01, corrected_p = fdrcorrection(pancan['stouffer-p'], alpha = 0.01)

  pancan.loc[:, 'rejected-0.05'] = rejected05
  pancan.loc[:, 'rejected-0.01'] = rejected01
  pancan.loc[:, 'corrected-p'] = corrected_p

  pancan.to_csv(os.path.join(outdir, p + '_univariate_fdr.csv'))


#%% P53-mutant multivariate zscores

parallel = 35
outdir = dropbox_dir + 'comprehensive-tcga-survival/p53-mutant-multivariate-zscores'

mutations_df = mutations.prep_data(platforms['mutations-non-synonymous'])
p53_muts = {}

for ctype in tcga_cdr_local.cancer_types():
  ctype_clinical = tcga_cdr_local.cancer_type_data(ctype)
  ctype_muts = mutations.ctype_cleaning(mutations_df, ctype, ctype_clinical)
  if '\'TP53' in list(ctype_muts.columns.values):
    p53_muts[ctype] = ctype_muts[['\'TP53']]
  else:
    p53_muts[ctype] = pd.DataFrame({'\'TP53': [0]*ctype_muts.shape[0]}, index=ctype_muts.index)
  p53_muts[ctype].columns =  [i[1:] if '\'' in i else i for i in p53_muts[ctype].columns]
  print(p53_muts[ctype].columns)
print('mutant selection complete')


for k in p53_muts.keys():
  p53_muts[k].columns = p53_muts[k].columns = [i + '_mut' for i in p53_muts[k].columns]

print(p53_muts['ACC'].columns)

for p, rd in platforms.items():
  platform_outdir = os.path.join(outdir, p)
  pathlib.Path(platform_outdir).mkdir(parents=True, exist_ok=True)
  if p == 'rppa':
    rppa.zscores(rd, clinical, platform_outdir, parallel, additional_vars=p53_muts)
  if p == 'cn':
    cn.zscores(rd, clinical, platform_outdir, 14, additional_vars=p53_muts)
  if p == 'rnaseq':
    rnaseq.zscores(rd, clinical, platform_outdir, parallel, additional_vars=p53_muts)
  if p == 'mirna':
    mirna.zscores(rd, clinical, platform_outdir, parallel, additional_vars=p53_muts)
  if p == 'methylation':
    methylation_data = rd[0]
    methylation_key = rd[1]
    methylation.zscores(methylation_data, methylation_key, clinical, platform_outdir, parallel, additional_vars=p53_muts)
  if p == 'mutations-non-synonymous':
    mutations.zscores(rd, clinical, platform_outdir, parallel, additional_vars=p53_muts)

