import pytest
import os
import pandas as pd
import numpy as np
import biomarker_survival as surv

from comprehensive_tcga_survival import mirna
from comprehensive_tcga_survival import ZscoreCommon

FIXTURE_DIR = os.path.join(
   os.path.dirname(os.path.realpath(__file__)),
      'test_data',
        )


def test_prep_data():
  mirna_path = os.path.join(FIXTURE_DIR, 'mirna_small.csv')
  mirna_data = mirna.prep_data(mirna_path).set_index('index')
  assert mirna_data.shape == (508, 743)
  assert np.sum((mirna_data < 0).values.ravel()) == 0 # log2 norm and clip means no negative values

def test_prep_and_clean_data():
  mirna_path = os.path.join(FIXTURE_DIR, 'mirna_small.csv')
  mirna_data = mirna.prep_data(mirna_path)
  mirna_ACC = mirna.ctype_cleaning(mirna_data, 'ACC', pd.DataFrame()) #ctype cleaning is unused, so empty is fine

  assert mirna_ACC.shape == (488, 743)
  assert 'TCGA-C4-A0F6' in list(mirna_ACC.index)
  assert np.sum((mirna_ACC < 0).values.ravel()) == 0 # log2 norm and clip means no negative values


def test_zscores(tmpdir):
  path = os.path.join(FIXTURE_DIR, 'TCGA-CDR-SupplementalTableS1-2019-05-27.xlsx')
  clin = surv.TCGA_CDR_util(path)

  mirna_path = os.path.join(FIXTURE_DIR, 'mirna_small.csv')
  mirna_data = mirna.prep_data(mirna_path)

  # 0th is the "index" col which is dropped during cleanup
  mirna_subset = mirna_data[mirna_data.columns[0:51]]
  assert mirna_subset.shape == (508, 51)

  mirna_zscores = ZscoreCommon(mirna.prep_data, mirna.ctype_cleaning)

  cox_ACC = mirna_zscores.cancer_type_cox(mirna_subset, clin.cancer_type_data('ACC'), 'ACC', str(tmpdir))
  assert cox_ACC.shape == (49, 6)
  assert 'hsa-let-7a-2-3p' in list(cox_ACC.index)
  print(cox_ACC.loc['hsa-let-7a-2-3p'])
  assert cox_ACC.loc['hsa-let-7a-2-3p'].round(6).to_dict() == {
    'hazard_ratio':       1.449653,
    'lower_conf':         1.046805,
    'n':                  79.000000,
    'p':                  0.025396,
    'upper_conf':         2.007532,
    'z':                  2.235327}

  cox_BLCA = mirna_zscores.cancer_type_cox(mirna_subset, clin.cancer_type_data('BLCA'), 'BLCA', str(tmpdir))
  assert cox_BLCA.shape == (50, 6)
  assert 'hsa-let-7a-2-3p' in list(cox_BLCA.index)
  print(cox_BLCA.loc['hsa-let-7a-2-3p'])
  assert cox_BLCA.loc['hsa-let-7a-2-3p'].round(6).to_dict() == {
    'hazard_ratio':       1.161701,
    'lower_conf':         1.050871,
    'n':                  405.000000,
    'p':                  0.003391,
    'upper_conf':         1.284219,
    'z':                  2.929917}

