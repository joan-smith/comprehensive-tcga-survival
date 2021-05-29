import pytest
import os
import pandas as pd
import biomarker_survival as surv

from comprehensive_tcga_survival import methylation
from comprehensive_tcga_survival import ZscoreCommon

FIXTURE_DIR = os.path.join(
   os.path.dirname(os.path.realpath(__file__)),
      'test_data',
        )


def test_prep_and_clean():
  methylation_path = os.path.join(FIXTURE_DIR, 'methylation_small.tsv')

  methylation_key_path = os.path.join(FIXTURE_DIR,
                                'methylation_to_genes_key_450_15017482_v1-2_small.csv')
  methylation_data = methylation.prep_data(methylation_path,
                                           extra_data=methylation_key_path)

  print(methylation_data.columns)
  assert methylation_data.shape == (12039, 20)
  assert 'TCGA-ZS-A9CD-01A-11D-A36Y-05' in list(methylation_data.index)
  assert 'HOPX' in methylation_data.columns

def test_prep_and_clean_with_output(tmpdir):
  methylation_path = os.path.join(FIXTURE_DIR, 'methylation_small.tsv')

  methylation_key_path = os.path.join(FIXTURE_DIR,
                                'methylation_to_genes_key_450_15017482_v1-2_small.csv')
  methylation_data = methylation.prep_data(methylation_path,
                                           extra_data=methylation_key_path,
                                           outdir=tmpdir)

  print(methylation_data.columns)
  assert methylation_data.shape == (12039, 20)
  assert len(os.listdir(tmpdir)) == 1
  assert os.listdir(tmpdir) == ['methylation_mean.csv']
  assert 'TCGA-ZS-A9CD-01A-11D-A36Y-05' in list(methylation_data.index)
  assert 'HOPX' in methylation_data.columns

def test_zscores(tmpdir):
  methylation_key_path = os.path.join(FIXTURE_DIR,
                                'methylation_to_genes_key_450_15017482_v1-2_small.csv')
  path = os.path.join(FIXTURE_DIR, 'TCGA-CDR-SupplementalTableS1-2019-05-27.xlsx')
  clin = surv.TCGA_CDR_util(path)

  methylation_path = os.path.join(FIXTURE_DIR, 'methylation_small.tsv')
  methylation_data = methylation.prep_data(methylation_path,
                                           extra_data=methylation_key_path)
  # methylation_subset = methylation_data[methylation_data.columns[0:50]]

  methylation_zscores = ZscoreCommon(methylation.prep_data,
                                     methylation.ctype_cleaning,
                                     extra_data=methylation_key_path)

  cox = methylation_zscores.cancer_type_cox(methylation_data,
                                            clin.cancer_type_data('ACC'),
                                            'ACC', str(tmpdir))
  print(cox)
  assert cox.shape == (20, 6)
  assert 'HOPX' in list(cox.index)
  assert cox.loc['HOPX'].round(6).to_dict() == {
    'hazard_ratio':       0.62198399999999998,
    'lower_conf':         0.124599,
    'n':                  79.00000,
    'p':                  0.56269499999999995,
    'upper_conf':         3.104873,
    'z':                  -0.578843}


