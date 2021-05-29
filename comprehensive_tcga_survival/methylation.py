#!/usr/bin/env python
# encoding: utf-8
'''
methylation.py

Created by Joan Smith
on 2019-10-10.

Copyright (c) 2019 All rights reserved.
'''

import pandas as pd
import argparse
import sys
import os

import biomarker_survival as surv

from . import ZscoreCommon

def get_options(argv):
  parser = argparse.ArgumentParser(description='Get methylation file, clinical file, optional output dir')
  parser.add_argument('-m', action='store', dest='methylation')
  parser.add_argument('-k', action='store', dest='methylation_key')
  parser.add_argument('-c', action='store', dest='tcga_cdr')
  parser.add_argument('-p', action='store', dest='parallel', type=int)
  parser.add_argument('-o', action='store', dest='output_directory', default='.')
  parser.add_argument('-no_zscores', action='store_false', dest='zscores')
  parser.set_defaults(zscores=True)
  ns = parser.parse_args()

  return ns.methylation, ns.methylation_key, ns.tcga_cdr, ns.output_directory, ns.parallel, ns.zscores


def ctype_cleaning(df, ctype, ctype_clinical): #ctype_clinical is unused
  df = df.reset_index()
  df = surv.maybe_clear_non_01s(df, 'index', ctype)
  df = surv.add_identifier_column(df, 'index')
  df = df.set_index('identifier')
  df = df.drop('index', axis=1)
  return df


# Extra data = the methylation to genes key path
def prep_data(methylation, extra_data=None, outdir=None):
  methylation_data = pd.read_csv(methylation, sep='\t', header=0,
                                 na_values='???', index_col=0)

  methylation_key = pd.read_csv(extra_data, header=7,
                                usecols=[0,21])
  # We expect that the key file does not contain any probe more than one time.
  # If it does, we'll need to modify the prep below to account for the same probe
  # appearing twice in the file.
  assert(methylation_key[methylation_key['IlmnID'].duplicated()].shape == (0,2))

  methylation_key = methylation_key.set_index('IlmnID')
  key = methylation_key['UCSC_RefGene_Name'].str.strip().str.split(';', expand=True)
  key = key.reset_index()
  key = pd.melt(key, id_vars='IlmnID', value_name='gene')
  key = key.drop('variable', axis=1)
  key = key.drop_duplicates().dropna(subset=['gene'])

  methylation_with_genes = key.join(methylation_data, how='inner', on='IlmnID')
  methylation_mean_by_gene = methylation_with_genes.groupby('gene')

  # Intermediate output for debugging/verification
  # methylation_mean_by_gene.count().to_csv('methylation_probe_counts.csv')
  # print methylation_mean_by_gene.get_group('HOPX')
  # methylation_mean_by_gene.get_group('HOPX').to_csv('methylation_HOPX.csv')

  methylation_mean = methylation_mean_by_gene.mean().T
  if outdir:
    methylation_mean.to_csv(os.path.join(outdir, 'methylation_mean.csv'))
  return methylation_mean

def metadata(methylation, methylation_key, clinical):
  methylation_zscores = ZscoreCommon(prep_data, ctype_cleaning, extra_data=methylation_key)
  return methylation_zscores.metadata(methylation, clinical)

def zscores(methylation, methylation_key, clinical, outdir, parallel, additional_vars):
  methylation_zscores = ZscoreCommon(prep_data, ctype_cleaning, extra_data=methylation_key)
  methylation_zscores.zscores(methylation, clinical, outdir, parallel_workers=parallel, additional_vars=additional_vars)

  pancan_df = surv.pancan(outdir, multivariate=(len(additional_vars) > 0))
  pancan_df.to_csv(os.path.join(outdir, 'pancan.csv'), index_label='gene')

def main(argv=None):
  methylation, methylation_key, clinical, outdir, parallel, zscores = get_options(argv)
  if zscores:
    zscores(methylation, methylation_key, clinical, outdir, parallel)
  else:
    # Don't calculate zscores, just transform the methylation data with the key file to produce methylation mean for all cancer types
    prep_data(methylation, extra_data=methylation_key, outdir=outdir)

if __name__ == "__main__":
  main()
