#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Dec 20 09:44:57 2020

@author: Joan Smith
"""

import pytest
import os
import pandas as pd
import biomarker_survival as surv

from comprehensive_tcga_survival import Multivariate

FIXTURE_DIR = os.path.join(
   os.path.dirname(os.path.realpath(__file__)),
      'test_data',
        )

def test_prep_multivar():
  stage_key_path = os.path.join(FIXTURE_DIR, 'Stage_key.xlsx')
  grade_key_path = os.path.join(FIXTURE_DIR, 'Grade_key.xlsx')
  path = os.path.join(FIXTURE_DIR, 'TCGA-CDR-SupplementalTableS1-2019-05-27.xlsx')
  clin = surv.TCGA_CDR_util(path)

  m = Multivariate(clin, stage_key_path, grade_key_path)
  cols = m.prep_ctype_multivar('COAD')

  print(cols)
  print(cols.columns)
  assert cols.shape == (457, 5)
  assert list(cols.columns).sort() == ['gender', 'race', 'age_at_initial_pathologic_diagnosis', 'GTE_Stage_IV', 'GTE_Stage_III', 'GTE_Stage_II'].sort()

def test_prep_ctype_multivar_esca():
  stage_key_path = os.path.join(FIXTURE_DIR, 'Stage_key.xlsx')
  grade_key_path = os.path.join(FIXTURE_DIR, 'Grade_key.xlsx')
  path = os.path.join(FIXTURE_DIR, 'TCGA-CDR-SupplementalTableS1-2019-05-27.xlsx')
  clin = surv.TCGA_CDR_util(path)

  m = Multivariate(clin, stage_key_path, grade_key_path)

  cols_esca = m.prep_ctype_multivar('ESCA')
  assert cols_esca.shape == (185, 4)
  assert list(cols_esca.columns).sort() == ['gender', 'race', 'age_at_initial_pathologic_diagnosis', 'GTE_G3', 'GTE_Stage_III'].sort()

def test_prep_all_multivar():
  stage_key_path = os.path.join(FIXTURE_DIR, 'Stage_key.xlsx')
  grade_key_path = os.path.join(FIXTURE_DIR, 'Grade_key.xlsx')
  path = os.path.join(FIXTURE_DIR, 'TCGA-CDR-SupplementalTableS1-2019-05-27.xlsx')
  clin = surv.TCGA_CDR_util(path)

  m = Multivariate(clin, stage_key_path, grade_key_path)
  multivars = m.prep_all_multivar()

  assert len(multivars) == 33
  assert 'ACC' in multivars.keys()
  assert 'BLCA' in multivars.keys()
