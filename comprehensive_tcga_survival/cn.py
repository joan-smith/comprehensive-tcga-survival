#!/usr/bin/env python
# encoding: utf-8
'''
cn.py

Created by Joan Smith
on 2019-10-10.

Copyright (c) 2019 All rights reserved.
'''

import pandas as pd
import numpy as np
import dask.dataframe as dd
import argparse
import sys
import os

import biomarker_survival as surv

from . import ZscoreCommon

from intervaltree import IntervalTree

def get_options(argv):
  parser = argparse.ArgumentParser(description='Get CN file, clinical file, optional output dir')
  parser.add_argument('-d', action='store', dest='cn')
  parser.add_argument('-k', action='store', dest='annotation_file')
  parser.add_argument('-c', action='store', dest='tcga_cdr')
  parser.add_argument('-p', action='store', dest='parallel', type=int)
  parser.add_argument('-o', action='store', dest='output_directory', default='.')
  parser.add_argument('--no_zscores', dest='zscores', action='store_false')
  parser.add_argument('--no_from_raw', dest='from_raw', action='store_false')
  parser.set_defaults(zscores=True, from_raw=True, parallel=1)
  ns = parser.parse_args()

  return ns.cn, ns.annotation_file, ns.tcga_cdr, ns.output_directory, ns.parallel, ns.zscores, ns.from_raw


def ctype_cleaning(df, ctype, ctype_clinical): #ctype_clinical is unused
  df = df.reset_index()
  df = surv.maybe_clear_non_01s(df, 'index', ctype)
  df = surv.add_identifier_column(df, 'index')
  if ctype == 'LAML':
    PRIMARY_BLOOD_DERIVED_CANCER_PATIENT_ID_REGEX = '^.{4}-.{2}-.{4}-03.*'
    non_primary_blood_barcodes = df[~df['index'].str.contains(PRIMARY_BLOOD_DERIVED_CANCER_PATIENT_ID_REGEX)]['index']
    df = df.drop(non_primary_blood_barcodes.index)
  df.set_index('identifier', inplace=True, drop=True)
  df = df.drop('index', axis=1)
  return df


def prep_annotation_file(annotation_path):
  gtf = pd.read_csv(annotation_path, sep='\t', header=None, comment="#")
  gtf.columns = ['chr', 'annotation_src', 'feature_type', 'start', 'end', 'score', 'genomic_strand', 'genomic_phase', 'extra']
  gtf['gene_name'] = gtf['extra'].str.extract(pat='gene_name "(.*?)";')
  gtf['transcript_type'] = gtf['extra'].str.extract(pat='transcript_type "(.*?)";')
  gtf = gtf[gtf['feature_type'] == 'transcript']

  gtf = gtf.groupby('gene_name').apply(lambda x: x[x['start'] == x['start'].min()].iloc[0])
  return gtf[['chr', 'start', 'end', 'feature_type', 'transcript_type']]


# Takes a raw segment file, returns a dictionary of lists of interval trees.
# Dictionary is keyed by patient id
# List is "chromosome copy number intervals". List is 24 long, has one slot for each chromosome
# Interval tree has copy number value at each interval on a particular chromosome, for each patient.
def prep_interval_trees(cn_data):
  patient_data = {}
  for index, row in cn_data.iterrows():
    chromosome = int(row['Chromosome'])
    # Add the new patient to the data dict, initializing the chromosome list
    if not index in patient_data:
      # note the length of the list is 24, so we can use chromosome number (and not worry about index 0)
      patient_data[index] = [0] * 24
      # print('New id: ', index)

    # Initialize the interval tree for the new chromosome for this patient
    if patient_data[index][chromosome] == 0:
      patient_data[index][chromosome] = IntervalTree()

    # Add the range and copy number in this row to the correct patient_data/chromosome location
    # Note the interval tree implementation uses half close intervals, but copy number data
    # uses closed intervals, so we add 1 to the end to ensure our intervaltree matches the data.
    start, end, copy_number = row['Start'], row['End'], row['Segment_Mean']
    patient_data[index][chromosome][start:end+1] = copy_number
  return patient_data

def get_copy_number_for_gene_start(chr_cn, annotation_start, patient, gene):
  if chr_cn == 0:
    return np.nan

  intervals = chr_cn[annotation_start]
  if len(intervals) > 1:
    print('Err: more than one interval found for patient:', patient, 'for gene:', gene)
    print('Got:', intervals)
    return np.nan
  if len(intervals) == 0:
    # print('Err: no interval found for patient:', patient, 'for gene:', gene)
    return np.nan
  else:
    interval = intervals.pop()
    return interval.data

def parse_chrom(chrom):
  if not 'chr' in chrom:
    return None
  if chrom[3:] in ['X', 'Y']:
    return 23
  elif chrom[3:] in ['M']:
    return 0
  else:
    return int(chrom[3:])

# Extra data = the annotation file
def prep_data(cn_path, extra_data=None, parallel=1, outdir=None, **kwargs):
  annotations = prep_annotation_file(extra_data)
  print('annotations prepped')
  cn_data = pd.read_csv(cn_path, sep='\t', header=0,
                                 na_values='???', index_col=0)

  cn_interval_trees = prep_interval_trees(cn_data)
  print('interval trees created')
  ddata = dd.from_pandas(annotations, npartitions=parallel, sort=False)
  patients = sorted(tuple(cn_interval_trees.keys()))
  cols = np.append(patients, ['chr', 'start'])
  meta_df = pd.DataFrame(index=annotations.index, columns=cols).astype(float)

  def make_cn_row(annotation_row):
    chrom = parse_chrom(annotation_row['chr'])
    cn_row = pd.Series(index=patients, name=annotation_row.name)
    for i, patient in enumerate(patients):
      cn_row.loc[patient] = get_copy_number_for_gene_start(cn_interval_trees[patient][chrom], annotation_row['start'], patient, annotation_row.name)
    cn_row['chr'] = chrom
    cn_row['start'] = annotation_row['start']
    return cn_row

  gene_cn_data = ddata.apply(make_cn_row, axis=1, meta=meta_df).compute(scheduler='processes')
  if outdir:
    annotations.to_csv(os.path.join(outdir, 'prepped_annotations.csv'))
    gene_cn_data.to_csv(os.path.join(outdir, 'cn_by_gene.csv'))
  else:
    gene_cn_data = gene_cn_data.drop(['chr', 'start'], axis=1)
  return gene_cn_data.T


def prep_processed_data(cn_path, **kwargs):
  return pd.read_csv(cn_path, header=0, index_col=0).T


def metadata(cn, clinical):
  cn_zscores = ZscoreCommon(prep_processed_data, ctype_cleaning)
  return cn_zscores.metadata(cn, clinical)

def zscores(cn, clinical, outdir, parallel, additional_vars={}):
  '''
  Run zscores on demand _not_ from raw.

  Parameters
  ----------
  cn : str
    pre-processed input file.
  clinical : str
    TCGA CDR.
  outdir : str
  parallel : int.
  additional_vars : dict
    Dictionary of additional variables to include in the cox proportional hazards
    Keyed by cancer type.

  Returns
  -------
  None.

  '''
  cn_zscores = ZscoreCommon(prep_processed_data, ctype_cleaning)
  cn_zscores.zscores(cn, clinical, outdir, parallel_workers=parallel, additional_vars=additional_vars)

  pancan_df = surv.pancan(outdir, multivariate=(len(additional_vars) > 0))
  pancan_df.to_csv(os.path.join(outdir, 'cn-pancan.csv'), index_label='gene')

def main(argv=None):
  cn, annotation_file, clinical, outdir, parallel, zscores, from_raw = get_options(argv)

  if not zscores and not from_raw:
    print('Error: --no_zscores and --no_from_raw are not compatible')
    print('To preprocess data, use --no_zscores alone. To run zscores from preprocessed data, use --no_from_raw')
  if zscores and from_raw:
    print('Zscores from raw data')
    cn_zscores = ZscoreCommon(prep_data, ctype_cleaning, extra_data=annotation_file)
    cn_zscores.zscores(cn, clinical, outdir, parallel_workers=parallel)

    pancan_df = surv.pancan(outdir)
    pancan_df.to_csv(os.path.join(outdir, 'pancan.csv'), index_label='gene')
  elif zscores and not from_raw:
    print('Zscores from pre-processed data')
    cn_zscores = ZscoreCommon(prep_processed_data, ctype_cleaning)
    cn_zscores.zscores(cn, clinical, outdir, parallel_workers=parallel, additional_vars={})
  else:
    # Don't calculate zscores, just transform the CN data and annotation file to produce a copy number by gene file, for all cancer types
    print('Process data only, no zscore calculations')
    prep_data(cn, extra_data=annotation_file, outdir=outdir, parallel=parallel)

if __name__ == "__main__":
  main()
