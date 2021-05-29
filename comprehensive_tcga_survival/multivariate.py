#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Dec 20 09:34:57 2020

@author: Joan Smith
"""
import pandas as pd
import argparse
import sys
import os

from comprehensive_tcga_survival import StageGradeParser
import biomarker_survival as surv

class Multivariate():
  def __init__(self, tcga_cdr, stage_key, grade_key):
    self.tcga_cdr = tcga_cdr
    self.stage = StageGradeParser(stage_key, tcga_cdr)
    self.grade = StageGradeParser(grade_key, tcga_cdr)

  def prep_ctype_multivar(self, ctype):
    gender_age =  self.tcga_cdr.cancer_type_data(ctype,
                                               extra_cols=['gender', 'age_at_initial_pathologic_diagnosis'])
    gender_age = gender_age[['gender', 'age_at_initial_pathologic_diagnosis']]

    gender_age = gender_age.replace({'FEMALE': 0, 'MALE': 1})
    stage_cols = self.stage.parse_for_ctype(ctype)
    grade_cols = self.grade.parse_for_ctype(ctype)

    multivar_cols = gender_age
    if not stage_cols.empty:
      multivar_cols = multivar_cols.join(stage_cols, how='inner')
    if not grade_cols.empty:
      multivar_cols = multivar_cols.join(grade_cols, how='inner')
    return multivar_cols


  # returns a dictionary of age, gender, and appropriate stage/grade features,
  # keyed by cancer type
  def prep_all_multivar(self):
    multivars = {}
    for c in self.tcga_cdr.cancer_types():
      multivars[c] = self.prep_ctype_multivar(c)

    return multivars
#%%%

def get_options(argv):
  parser = argparse.ArgumentParser(description='Get rppa file, clinical file, optional output dir')
  parser.add_argument('-c', action='store', dest='tcga_cdr')
  parser.add_argument('-s', action='store', dest='stage_key')
  parser.add_argument('-g', action='store', dest='grade_key')
  parser.add_argument('-o', action='store', dest='output_directory', default='.')
  ns = parser.parse_args()

  return ns.tcga_cdr, ns.stage_key, ns.grade_key, ns.output_directory


def multivar_main(argv=None):
  clinical, stage_key, grade_key, outdir = get_options(argv)

  tcga_cdr = surv.TCGA_CDR_util(clinical)

  m = Multivariate(tcga_cdr, stage_key, grade_key)

  cox_dicts = {}
  for ctype in tcga_cdr.cancer_types():
    print(ctype)
    data_cols = m.prep_multivar(ctype)
    time_censor = tcga_cdr.cancer_type_data(ctype)
    ctype_patient_count = time_censor.shape[0]
    df = time_censor.join(data_cols, how='inner').dropna(how='any')

    cox_dict = surv.do_multivariate_cox(df['time'], df['censor'], df['gender'], df.drop(['time', 'censor', 'gender'], axis=1))
    cox_dict['total patient count'] = ctype_patient_count
    cox_dicts[ctype] = pd.Series(cox_dict)


  cox_df = pd.DataFrame(cox_dicts)
  cox_df.index = ([i.replace('var', 'gender') for i in list(cox_df.index)])


  cox_df.to_csv(os.path.join(outdir, 'clinical_multivariate.csv'))


if __name__ == "__main__":
  multivar_main()
