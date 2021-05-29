#!/usr/bin/env python
# encoding: utf-8
'''
mirna.py

Created by Joan Smith
on 2019-8-29.

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
  parser = argparse.ArgumentParser(description='Get mirna file, clinical file, optional output dir')
  parser.add_argument('-m', action='store', dest='mirna')
  parser.add_argument('-c', action='store', dest='tcga_cdr')
  parser.add_argument('-p', action='store', dest='parallel', type=int)
  parser.add_argument('-o', action='store', dest='output_directory', default='.')
  ns = parser.parse_args()

  return ns.mirna, ns.tcga_cdr, ns.output_directory, ns.parallel

def prep_data(mirna_path, extra_data=None):
  mirna = pd.read_csv(mirna_path, header=0, na_values='???', index_col=0)
  mirna = mirna.drop('Correction', axis=1)

  mirna_log2 = mirna.apply(np.log2)
  mirna_clipped_log2 = mirna_log2.clip(lower=0)

  return mirna_clipped_log2.T.reset_index()

def ctype_cleaning(df, ctype, ctype_ctype_clinical): #ctype_clinical unused
  df = surv.maybe_clear_non_01s(df, 'index', ctype)
  df = surv.add_identifier_column(df, 'index')
  df = df.set_index('identifier')
  df = df.drop('index', axis=1)
  return df

def metadata(mirna, clinical):
  mirna_zscores = ZscoreCommon(prep_data, ctype_cleaning)
  return mirna_zscores.metadata(mirna, clinical)

def zscores(mirna, clinical, outdir, parallel, additional_vars={}):
  mirna_zscores = ZscoreCommon(prep_data, ctype_cleaning)
  mirna_zscores.zscores(mirna, clinical, outdir, parallel_workers=parallel, additional_vars=additional_vars)

  pancan_df = surv.pancan(outdir, multivariate=(len(additional_vars) > 0))
  pancan_df.to_csv(os.path.join(outdir, 'pancan.csv'), index_label='gene')



def main(argv=None):
  mirna, clinical, outdir, parallel = get_options(argv)
  zscores(mirna, clinical, outdir, parallel)

if __name__ == "__main__":
  main()
