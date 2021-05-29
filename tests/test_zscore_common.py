import pytest
import os
import pandas as pd
import rpy2
import biomarker_survival as surv

from comprehensive_tcga_survival import rppa
from comprehensive_tcga_survival import ZscoreCommon

FIXTURE_DIR = os.path.join(
   os.path.dirname(os.path.realpath(__file__)),
      'test_data',
        )

rpy2.robjects.r['options'](warn=-1)

def test_zscores_serial(tmpdir):
  clin_path = os.path.join(FIXTURE_DIR, 'TCGA-CDR-SupplementalTableS1-2019-05-27.xlsx')
  rppa_path = os.path.join(FIXTURE_DIR, 'rppa_small.txt')

  rppa_zscores = ZscoreCommon(rppa.prep_data, rppa.ctype_cleaning)
  cox_all = rppa_zscores.zscores(rppa_path, clin_path, str(tmpdir))

  print(cox_all)
  print(len(cox_all))
  assert len(cox_all) == 33
  cox_ACC = cox_all[0]
  assert cox_ACC.shape == (194, 6)
  assert 'AKT' in list(cox_ACC.index)
  assert cox_ACC.loc['AKT'].round(6).to_dict() == {
    'hazard_ratio':       0.609382,
    'lower_conf':         0.212732,
    'n':                  46.000000,
    'p':                  0.356296,
    'upper_conf':         1.745601,
    'z':                  -0.922446}


def test_zscores_parallel(tmpdir):
  clin_path = os.path.join(FIXTURE_DIR, 'TCGA-CDR-SupplementalTableS1-2019-05-27.xlsx')
  rppa_path = os.path.join(FIXTURE_DIR, 'rppa_small.txt')

  rppa_zscores = ZscoreCommon(rppa.prep_data, rppa.ctype_cleaning)
  cox_all = rppa_zscores.zscores(rppa_path, clin_path, str(tmpdir), parallel_workers=2)

  print(cox_all)
  print(len(cox_all))
  assert len(cox_all) == 33
  cox_ACC = cox_all[0]
  assert cox_ACC.shape == (194, 6)
  assert 'AKT' in list(cox_ACC.index)
  assert cox_ACC.loc['AKT'].round(6).to_dict() == {
    'hazard_ratio':       0.609382,
    'lower_conf':         0.212732,
    'n':                  46.000000,
    'p':                  0.356296,
    'upper_conf':         1.745601,
    'z':                  -0.922446}
