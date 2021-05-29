#!/usr/bin/env python
# encoding: utf-8
'''
rnaseq.py

Created by Joan Smith
on 2019-7-07.

Copyright (c) 2019 All rights reserved.
'''

import pandas as pd
import numpy as np
import argparse
import sys
import os

import biomarker_survival as surv

from .zscore_common import ZscoreCommon

def get_options(argv):
  parser = argparse.ArgumentParser(description='Get rnaseq file, clinical file, optional output dir')
  parser.add_argument('-r', action='store', dest='rnaseq')
  parser.add_argument('-c', action='store', dest='tcga_cdr')
  parser.add_argument('-p', action='store', dest='parallel', type=int)
  parser.add_argument('-o', action='store', dest='output_directory', default='.')
  ns = parser.parse_args()

  return ns.rnaseq, ns.tcga_cdr, ns.output_directory, ns.parallel

def ctype_cleaning(df, ctype, ctype_clinical): #ctype_clinical is unused
  df = surv.maybe_clear_non_01s(df, 'index', ctype)
  df = surv.add_identifier_column(df, 'index')
  df = df.set_index('identifier')
  df = df.drop('index', axis=1)
  return df

def prep_data(rnaseq, extra_data=None):
  rnaseq = pd.read_csv(rnaseq, sep='\t', header=0,
                      na_values='???')

  genes = rnaseq['gene_id'].str.split('|', expand=True)
  dupes = genes[genes.duplicated(subset=0, keep=False)]
  genes.loc[genes.duplicated(subset=0, keep=False),0] = dupes[0] + '|' + dupes[1]
  rnaseq['gene_id'] = genes[0]
  rnaseq = rnaseq.set_index('gene_id')

  rnaseq_log2 = rnaseq.apply(np.log2)
  rnaseq_clipped_log2 = rnaseq_log2.clip(lower=0)

  return rnaseq_clipped_log2.T.reset_index()

def metadata(rnaseq, clinical):
  rnaseq_zscores = ZscoreCommon(prep_data, ctype_cleaning)
  return rnaseq_zscores.metadata(rnaseq, clinical)

def zscores(rnaseq, clinical, outdir, parallel, additional_vars={}):
  rnaseq_zscores = ZscoreCommon(prep_data, ctype_cleaning)
  rnaseq_zscores.zscores(rnaseq, clinical, outdir, parallel_workers=parallel, additional_vars=additional_vars)

  pancan_df = surv.pancan(outdir, multivariate=(len(additional_vars) > 0))
  pancan_df.to_csv(os.path.join(outdir, 'pancan.csv'), index_label='gene')


def main(argv=None):
  rnaseq, clinical, outdir, parallel = get_options(argv)
  zscores(rnaseq, clinical, outdir, parallel)

if __name__ == "__main__":
  main()
