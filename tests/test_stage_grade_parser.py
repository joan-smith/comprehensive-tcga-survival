#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Dec 19 21:52:15 2020

@author: Joan Smith
"""

import pytest
import os
import pandas as pd
import biomarker_survival as surv

from comprehensive_tcga_survival import StageGradeParser

FIXTURE_DIR = os.path.join(
   os.path.dirname(os.path.realpath(__file__)),
      'test_data',
        )

def test_generate_comparison_dict():
    key_path = os.path.join(FIXTURE_DIR, 'Stage_key.xlsx')
    path = os.path.join(FIXTURE_DIR, 'TCGA-CDR-SupplementalTableS1-2019-05-27.xlsx')
    clin = surv.TCGA_CDR_util(path)

    key = StageGradeParser(key_path, clin)

    comparison_dict = key.generate_comparison_dict('COAD')
    assert len(comparison_dict.keys()) == 3
    assert list(comparison_dict.keys()).sort() == ['GTE_Stage_II', 'GTE_Stage_III', 'GTE_Stage_IV'].sort()
    assert len(comparison_dict['GTE_Stage_III']) > len(comparison_dict['GTE_Stage_IV'])
    assert len(comparison_dict['GTE_Stage_III']) == 7

    comparison_dict = key.generate_comparison_dict('ESCA')
    assert len(comparison_dict.keys()) == 1
    assert list(comparison_dict.keys()).sort() == ['GTE_Stage_III'].sort()

def test_missing_cancer_type():
  key_path = os.path.join(FIXTURE_DIR, 'Stage_key.xlsx')
  path = os.path.join(FIXTURE_DIR, 'TCGA-CDR-SupplementalTableS1-2019-05-27.xlsx')
  clin = surv.TCGA_CDR_util(path)

  key = StageGradeParser(key_path, clin)

  df = key.parse_for_ctype('GBM')
  assert df.empty

def test_case_insensitive():
  key_path = os.path.join(FIXTURE_DIR, 'Stage_key.xlsx')
  path = os.path.join(FIXTURE_DIR, 'TCGA-CDR-SupplementalTableS1-2019-05-27.xlsx')
  clin = surv.TCGA_CDR_util(path)

  key = StageGradeParser(key_path, clin)
  df = key.parse_for_ctype('THYM')
  assert df.shape == (123, 1)



def test_read_file():
  key_path = os.path.join(FIXTURE_DIR, 'Stage_key.xlsx')
  path = os.path.join(FIXTURE_DIR, 'TCGA-CDR-SupplementalTableS1-2019-05-27.xlsx')
  clin = surv.TCGA_CDR_util(path)

  key = StageGradeParser(key_path, clin)
  df = key.parse_for_ctype('COAD')

  print( df['GTE_Stage_III'].value_counts())
  assert df.shape == (457, 3)
  assert df['GTE_Stage_III'].value_counts()[0] == 263
  assert df['GTE_Stage_III'].value_counts()[1] == 194
