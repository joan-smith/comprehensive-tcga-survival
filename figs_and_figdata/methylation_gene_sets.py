#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan  1 15:07:42 2021

@author: Joan Smith
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan  1 11:12:02 2021

@author: Joan Smith
"""
import pandas as pd
import os
import pathlib

import biomarker_survival as surv
from comprehensive_tcga_survival import methylation
from comprehensive_tcga_survival import StageGradeParser


#%%

dropbox_dir = '~/Dropbox/comprehensive-tcga-survival/'
clinical = dropbox_dir + 'raw-data/TCGA-CDR-SupplementalTableS1-2019-05-27.xlsx'
methylation_input = dropbox_dir + 'raw-data/jhu-usc.edu_PANCAN_merged_HumanMethylation27_HumanMethylation450.betaValue_whitelisted.tsv'
methylation_key = dropbox_dir + 'raw-data/HumanMethylation450_15017482_v1-2.csv'

#%%

tcga_cdr = surv.TCGA_CDR_util(clinical)
cancer_types = tcga_cdr.cancer_types()

#%%
gene_list = pd.read_excel(os.path.join(dropbox_dir, 'Gene Ontology', 'Methylation events', 'Meth gene sets.xlsx'), engine='openpyxl')
gene_list = gene_list.apply(lambda x: x.str.replace('\'', '')) # remove quotes to facilitate matching. Will add them back in the end
df = methylation.prep_data(methylation_input, extra_data=methylation_key)

#%%
outdir = os.path.join(dropbox_dir, 'Gene Ontology', 'Methylation events', 'mean normalization')
for c in cancer_types:
  pathlib.Path(outdir).mkdir(parents=True, exist_ok=True)
  ctype_data = methylation.ctype_cleaning(df, c, None)
  print(ctype_data.shape)
  ctype_data = ctype_data.join(tcga_cdr.cancer_type_data(c), how='inner')
  for l in gene_list:
    print(l)
    present_genes = set(gene_list[l].dropna().values) & set(ctype_data.columns.values)
    ctype_genes = ctype_data.loc[:, present_genes]
    mean_centered = ctype_genes - ctype_genes.mean()
    mean_centered[['time', 'censor']] = ctype_data[['time', 'censor']]
    mean_centered.columns = '\'' + mean_centered.columns
    mean_centered.to_csv(os.path.join(outdir, c+ '_' + l + '_methylation_genes.csv'))
#%%
grade_key = dropbox_dir + 'Multivariate/Grade_key.xlsx'
grades = StageGradeParser(grade_key, tcga_cdr)
#%%
outdir = os.path.join(dropbox_dir, 'Gene Ontology', 'Methylation events', 'Grade normalization')

indir = dropbox_dir + 'Gene Ontology/Methylation events/'
gene_list = pd.read_excel(indir + 'SUZ12 adverse meth sites.xlsx', engine='openpyxl')
gene_list[gene_list.columns[0]] = [i[1:] if '\'' in i else i for i in gene_list[gene_list.columns[0]]]
#%%

outdir = indir + '/SUZ12'
pathlib.Path(outdir).mkdir(parents=True, exist_ok=True)
for c in cancer_types:
  grade_df = grades.parse_for_ctype(c)
  if grade_df.empty:
    continue

  ctype_data = methylation.ctype_cleaning(df, c, None)
  ctype_data = ctype_data.join(tcga_cdr.cancer_type_data(c), how='inner')

  gene_lists = {}
  for l in gene_list:
    print(l)
    present_genes = set(list(gene_list[l].values)) & set(ctype_data.columns.values)
    ctype_genes = ctype_data.loc[:, present_genes]
    mean_centered = ctype_genes - ctype_genes.mean()
    per_patient_mean = mean_centered.mean(axis=1)
    gene_lists[l] = per_patient_mean

  ctype_out = pd.concat(gene_lists, axis=1)
  ctype_out[['time', 'censor']] = ctype_data[['time', 'censor']]
  ctype_out = ctype_out.join(grade_df)
  ctype_out.to_csv(os.path.join(outdir, c + '_suz12_adverse_meth_sites_grade.csv'))


