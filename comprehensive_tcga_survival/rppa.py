#!/usr/bin/env python
# encoding: utf-8
'''
rppa.py

Created by Joan Smith
on 2019-8-29.

Copyright (c) 2019 All rights reserved.
'''

import pandas as pd
import argparse
import sys
import os

import biomarker_survival as surv

from .zscore_common import ZscoreCommon

def get_options(argv):
  parser = argparse.ArgumentParser(description='Get rppa file, clinical file, optional output dir')
  parser.add_argument('-r', action='store', dest='rppa')
  parser.add_argument('-c', action='store', dest='tcga_cdr')
  parser.add_argument('-p', action='store', dest='parallel', type=int)
  parser.add_argument('-o', action='store', dest='output_directory', default='.')
  ns = parser.parse_args()

  return ns.rppa, ns.tcga_cdr, ns.output_directory, ns.parallel


def ctype_cleaning(df, ctype, ctype_clinical): #ctype_clinical is unused
  ctype_filter = ctype
  if ctype in ['COAD', 'READ']:
    ctype_filter = 'CORE'
  df = df[df['TumorType'] == ctype_filter]
  df = surv.maybe_clear_non_01s(df, 'SampleID', ctype)
  df = surv.add_identifier_column(df, 'SampleID')
  df = df.set_index('identifier')
  df = df.drop(['SampleID', 'TumorType'], axis=1)
  return df


def prep_data(rppa, extra_data=None):
  rppa_data = pd.read_csv(rppa, sep='\t', header=0, na_values='???')
  return rppa_data


def zscores(rppa, clinical, outdir, parallel, additional_vars={}):
  rppa_zscores = ZscoreCommon(prep_data, ctype_cleaning)
  rppa_zscores.zscores(rppa, clinical, outdir, parallel_workers=parallel, additional_vars=additional_vars)

  pancan_df = surv.pancan(outdir, multivariate=(len(additional_vars) > 0))
  pancan_df.to_csv(os.path.join(outdir, 'pancan.csv'), index_label='gene')

def metadata(rppa, clinical):
  rppa_zscores = ZscoreCommon(prep_data, ctype_cleaning)
  return rppa_zscores.metadata(rppa, clinical)


def main(argv=None):
  rppa, clinical, outdir, parallel = get_options(argv)
  zscores(rppa, clinical, outdir, parallel)

if __name__ == "__main__":
  main()
