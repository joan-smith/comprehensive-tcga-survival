#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Dec 26 17:56:46 2020

@author: Joan Smith
"""
import pandas as pd
import os

import biomarker_survival as surv

from comprehensive_tcga_survival import ZscoreCommon
from comprehensive_tcga_survival import mutations


#%%
dropbox_dir = '~/Dropbox/'

mutations_file = dropbox_dir + 'comprehensive-tcga-survival/raw-data/mc3.v0.2.8.PUBLIC.maf'
clinical = dropbox_dir + 'comprehensive-tcga-survival/raw-data/TCGA-CDR-SupplementalTableS1-2019-05-27.xlsx'
point_mutations = dropbox_dir + 'comprehensive-tcga-survival/Driver_genes/point_mutations.xlsx'
outdir = dropbox_dir + 'comprehensive-tcga-survival/Driver_genes/point-mutations-gene'
MUTATION_PERCENT = 0.02
#%%
point_mut_data = pd.read_excel(point_mutations, engine='openpyxl', sheet_name=1)

#%% per-gene driver mutation zscores (e.g. collapse all the relevant driver mutations into their appropriate gene, then do zscores)
outdir = dropbox_dir + 'comprehensive-tcga-survival/Driver_genes/point-mutations-gene'

def prep_data(mutations_file, extra_data):
  mut_data = mutations.prep_data(mutations_file)
  relevant_mut_data = mut_data[mut_data['Hugo_Symbol'].isin(extra_data['Gene'])]
  relevant_mut_data = relevant_mut_data[relevant_mut_data['HGVSp_Short'].isin(extra_data['Mutation'])].copy(deep=True)

  #slightly gross hack to keep all the patients in the data, without including any of the uninteresting mutations
  all_patients = mut_data[['identifier', 'Tumor_Sample_Barcode']].copy(deep=True)
  all_patients.loc[:,'Variant_Classification'] = 'Silent'
  all_patients = all_patients.drop_duplicates(subset='identifier', keep='first')
  relevant_mut_data = pd.concat([relevant_mut_data, all_patients])

  return relevant_mut_data.copy(deep=True)

def ctype_cleaning(mutation_df, ctype, ctype_clinical):
  print(ctype)
  ctype_cleaned = mutations.ctype_cleaning(mutation_df, ctype, ctype_clinical).copy(deep=True)
  return ctype_cleaned
#%%
mutation_zscores = ZscoreCommon(prep_data, ctype_cleaning,
                                extra_data=point_mut_data,
                                use_presence_percentage_of_patients_threshold=True,
                                threshold_percent=MUTATION_PERCENT)
mutation_zscores.zscores(mutations_file, clinical, outdir, parallel_workers=0)

pancan_df = surv.pancan(outdir, multivariate=False)
pancan_df.to_csv(os.path.join(outdir, 'mutations-pancan.csv'), index_label='gene')

#%% per-driver mutation zscores

outdir = dropbox_dir + 'comprehensive-tcga-survival/Driver_genes/point-mutations'

def prep_data_uncollapsed_by_gene(mutations, extra_data):
  relevant_mut_data = prep_data(mutations, extra_data)
  data_by_mutation = relevant_mut_data.copy(deep=True)
  data_by_mutation.loc[:,'Hugo_Symbol'] = data_by_mutation['Hugo_Symbol'] + '_' + data_by_mutation['HGVSp_Short']
  return data_by_mutation
#%%

mutation_zscores = ZscoreCommon(prep_data_uncollapsed_by_gene, ctype_cleaning,
                                extra_data=point_mut_data,
                                use_presence_percentage_of_patients_threshold=True,
                                threshold_percent=MUTATION_PERCENT)

mutation_zscores.zscores(mutations_file, clinical, outdir, parallel_workers=0)
pancan_df = surv.pancan(outdir, multivariate=False)
pancan_df.to_csv(os.path.join(outdir, 'mutations-pancan.csv'), index_label='gene')

#%%

outdir = dropbox_dir + 'comprehensive-tcga-survival/Driver_genes/point-mutations-no-cutoff'

mutation_zscores = ZscoreCommon(prep_data_uncollapsed_by_gene, ctype_cleaning,
                                extra_data=point_mut_data,
                                use_presence_percentage_of_patients_threshold=True,
                                threshold_percent=0.0)

mutation_zscores.zscores(mutations_file, clinical, outdir, parallel_workers=0)
pancan_df = surv.pancan(outdir, multivariate=False)
pancan_df.to_csv(os.path.join(outdir, 'mutations-pancan.csv'), index_label='gene')

#%% KMPlot outputs

tcga_cdr = surv.TCGA_CDR_util(clinical)

data_by_gene = prep_data(mutations_file, point_mut_data)
print('data prepped')
data_by_mutation= data_by_gene.copy(deep=True)
print('by_mutation copy')
data_by_mutation.loc[:,'Hugo_Symbol'] = data_by_mutation['Hugo_Symbol'] + '_' + data_by_mutation['HGVSp_Short']
print('by mutation renamed')

outdir_by_mut = dropbox_dir + 'comprehensive-tcga-survival/Driver_genes/point-mutation-kmplot-data'
outdir_by_gene = dropbox_dir + 'comprehensive-tcga-survival/Driver_genes/point-mutation-kmplot-data-hotspot-mutations-per-gene'
#%%
for c in tcga_cdr.cancer_types():

  df_by_mut = ctype_cleaning(data_by_mutation, c, tcga_cdr.cancer_type_data(c))
  df_by_mut.to_csv(os.path.join(outdir_by_mut, c + '.csv'))

  df_by_gene = ctype_cleaning(data_by_gene, c, tcga_cdr.cancer_type_data(c))
  df_by_gene.to_csv(os.path.join(outdir_by_gene, c +'_by_gene.csv'))


#%% % of KRAS Mutations in 12, 13, 61

prepped_data = mutation_zscores.metadata(mutations_file, clinical, return_data=True)
kras = prepped_data.loc[:,prepped_data.columns.str.contains('KRAS')]

total_kras_mutations = kras.sum().sum()
kras_12 = kras.sum().loc[kras.sum().index.str.contains('12')].sum()
kras_13 = kras.sum().loc[kras.sum().index.str.contains('13')].sum()
kras_61 = kras.sum().loc[kras.sum().index.str.contains('61')].sum()

print('Percent KRAS mutations in 12, 13, 61: ', (kras_12+kras_13+kras_61)/total_kras_mutations)

