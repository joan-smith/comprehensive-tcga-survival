import pytest
import os
import pandas as pd
import biomarker_survival as surv

from comprehensive_tcga_survival import rppa
from comprehensive_tcga_survival import ZscoreCommon

FIXTURE_DIR = os.path.join(
   os.path.dirname(os.path.realpath(__file__)),
      'test_data',
        )


def test_prep_and_clean():
  rppa_path = os.path.join(FIXTURE_DIR, 'rppa_small.txt')
  rppa_data = rppa.prep_data(rppa_path)
  rppa_ACC = rppa.ctype_cleaning(rppa_data, 'ACC', pd.DataFrame())

  assert rppa_ACC.shape == (46, 198)
  assert rppa_data.shape == (46, 200) # includes columns for "SampleID" and "TumorType"
  assert 'TCGA-PA-A5YG' in list(rppa_ACC.index)


def test_zscores(tmpdir):
  path = os.path.join(FIXTURE_DIR, 'TCGA-CDR-SupplementalTableS1-2019-05-27.xlsx')
  clin = surv.TCGA_CDR_util(path)

  rppa_path = os.path.join(FIXTURE_DIR, 'rppa_small.txt')
  rppa_data = rppa.prep_data(rppa_path)

  # 0th is the "index" col which is dropped during cleanup
  rppa_subset = rppa_data[rppa_data.columns[0:52]]


  rppa_zscores = ZscoreCommon(rppa.prep_data, rppa.ctype_cleaning)
  cox = rppa_zscores.cancer_type_cox(rppa_subset, clin.cancer_type_data('ACC'), 'ACC', str(tmpdir))

  assert cox.shape == (50, 6)
  assert 'AKT' in list(cox.index)
  assert cox.loc['AKT'].round(6).to_dict() == {
    'hazard_ratio':       0.609382,
    'lower_conf':         0.212732,
    'n':                  46.000000,
    'p':                  0.356296,
    'upper_conf':         1.745601,
    'z':                  -0.922446}


def test_multivar_zscore(tmpdir):
  path = os.path.join(FIXTURE_DIR, 'TCGA-CDR-SupplementalTableS1-2019-05-27.xlsx')
  clin = surv.TCGA_CDR_util(path)

  rppa_path = os.path.join(FIXTURE_DIR, 'rppa_small.txt')
  rppa_data = rppa.prep_data(rppa_path)

  # 0th is the "index" col which is dropped during cleanup
  rppa_subset = rppa_data[rppa_data.columns[0:52]]


  extra_cols = clin.cancer_type_data('ACC', extra_cols=['gender'])[['gender']]
  extra_cols = extra_cols.replace({'MALE': 0, 'FEMALE': 1})

  rppa_zscores = ZscoreCommon(rppa.prep_data, rppa.ctype_cleaning)
  cox = rppa_zscores.cancer_type_cox_multivar(rppa_subset, clin.cancer_type_data('ACC'), 'ACC', str(tmpdir), extra_cols)

  print(tmpdir)
  assert cox.shape == (50, 11)
  assert 'AKT' in list(cox.index)
  assert cox.loc['AKT'].round(6).to_dict() == {
    'gender-p': 0.104235,
    'gender-z': 1.624659,
    'gender_hazard_ratio': 2.893846,
    'gender_lower_conf': 0.803079,
    'gender_upper_conf': 10.427789,
    'var-n': 46.0,
    'var-p': 0.320603,
    'var-z': -0.99322,
    'var_hazard_ratio': 0.598101,
    'var_lower_conf': 0.216907,
    'var_upper_conf': 1.649208
   }
