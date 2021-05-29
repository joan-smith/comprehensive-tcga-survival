#!/usr/bin/env python
# encoding: utf-8
'''
mutations.py

Created by Joan Smith
on 2019-6-03.

Copyright (c) 2019 All rights reserved.
'''

import pandas as pd
import argparse
import sys
import os

import biomarker_survival as surv
from .zscore_common import ZscoreCommon



MUTATION_PERCENT = .02

def ctype_cleaning(mutation_df, ctype, ctype_clinical):
  non_synonymous_classifications = ['Missense_Mutation', 'Nonsense_Mutation', 'Frame_Shift_Del',
                                    'Splice_Site', 'Frame_Shift_Ins', 'In_Frame_Del',
                                    'Translation_Start_Site','Nonstop_Mutation',
                                    'In_Frame_Ins']

  df, patients_w_sequence_and_clinical, num_patients = surv.prep_data(mutation_df, ctype_clinical, ctype)
  df_non_synonymous = df.where(df['Variant_Classification'].isin(non_synonymous_classifications)).copy(deep=True)
  df_non_synonymous = df_non_synonymous.dropna(subset=['Variant_Classification'])

  df_non_synonymous = df_non_synonymous.reset_index()
  df_non_synonymous.loc[:,'Hugo_Symbol'] = '\'' + df_non_synonymous['Hugo_Symbol'].astype(str)

  pivoted = surv.pivot_mutation_list_to_mutations_per_gene(df_non_synonymous)
  final = pd.concat([pivoted, patients_w_sequence_and_clinical[['time', 'censor']]], axis=1).fillna(0).drop(['time', 'censor'], axis=1)
  return final

def prep_data(mutation, **kwargs):
  df = surv.prep_mutation_data_alone(mutation)
  df = df[['Hugo_Symbol', 'Variant_Classification', 'Chromosome', 'Start_Position',
           'End_Position', 'Strand', 'HGVSp_Short', 'Tumor_Sample_Barcode', 'identifier']]
  return df.copy(deep=True)

def zscores(mutation, clinical, outdir, parallel, additional_vars={}):
  mutation_zscores = ZscoreCommon(prep_data, ctype_cleaning,
                                  use_presence_percentage_of_patients_threshold=True,
                                  threshold_percent=MUTATION_PERCENT)
  mutation_zscores.zscores(mutation, clinical, outdir, parallel_workers=parallel, additional_vars=additional_vars)

  pancan_df = surv.pancan(outdir, multivariate=(len(additional_vars) > 0))
  pancan_df.to_csv(os.path.join(outdir, 'mutations-pancan.csv'), index_label='gene')

def get_options(argv):
  parser = argparse.ArgumentParser(description='Get mutation file, clinical file, optional output dir')
  parser.add_argument('-m', action='store', dest='mutation')
  parser.add_argument('-c', action='store', dest='tcga_cdr')
  parser.add_argument('-p', action='store', dest='parallel_workers', type=int)
  parser.add_argument('-o', action='store', dest='output_directory', default='.')
  ns = parser.parse_args()

  return ns.mutation, ns.tcga_cdr, ns.output_directory, ns.parallel_workers


def metadata(mutation, clinical):
  tcga_cdr = surv.TCGA_CDR_util(clinical)
  muts = surv.prep_data(mutation)

  metadata = {}
  for ctype in tcga_cdr.cancer_types():
    ctype_clinical = tcga_cdr.cancer_type_data(ctype)
    ctype_muts = surv.prep_mutation_data(muts, ctype_clinical, ctype)
    df = ctype_muts.join(ctype_clinical, how='inner')
    metadata[ctype] = df.shape[0]
  return metadata

def main(argv=None):
  mutation, clinical, outdir, p = get_options(argv)
  zscores(mutation, clinical, outdir, p)

if __name__ == "__main__":
  main()
