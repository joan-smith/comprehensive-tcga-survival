import pytest
import os
import pandas as pd
import numpy as np
import biomarker_survival as surv

from comprehensive_tcga_survival import rnaseq
from comprehensive_tcga_survival import ZscoreCommon

FIXTURE_DIR = os.path.join(
   os.path.dirname(os.path.realpath(__file__)),
      'test_data',
        )


def test_prep_and_clean():
  rnaseq_path = os.path.join(FIXTURE_DIR, 'rnaseq_small.tsv')
  rnaseq_data = rnaseq.prep_data(rnaseq_path)
  rnaseq_ACC = rnaseq.ctype_cleaning(rnaseq_data, 'ACC', pd.DataFrame()) #clinical data is unneeded.

  assert rnaseq_ACC.shape == (323, 20531)
  assert rnaseq_data.shape == (337, 20532) # includes column for "index"
  assert 'TCGA-OR-A5J5' in list(rnaseq_ACC.index)
  assert np.sum((rnaseq_ACC < 0).values.ravel()) == 0 # log2 norm and clip means no negative values


def test_zscores(tmpdir):
  path = os.path.join(FIXTURE_DIR, 'TCGA-CDR-SupplementalTableS1-2019-05-27.xlsx')
  clin = surv.TCGA_CDR_util(path)

  rnaseq_path = os.path.join(FIXTURE_DIR, 'rnaseq_small.tsv')
  rnaseq_data = rnaseq.prep_data(rnaseq_path)

  # 0th is the "index" col which is dropped during cleanup
  rnaseq_subset = rnaseq_data[rnaseq_data.columns[0:51]]

  rnaseq_zscores = ZscoreCommon(rnaseq.prep_data, rnaseq.ctype_cleaning)
  cox = rnaseq_zscores.cancer_type_cox(rnaseq_subset, clin.cancer_type_data('BLCA'), 'BLCA', str(tmpdir))

  assert cox.shape == (50, 6)
  assert '?|100130426' in list(cox.index)
  print(cox.loc['?|100130426'])
  assert cox.loc['?|100130426'].round(6).to_dict() == {
    'hazard_ratio':       0.000002,
    'lower_conf':         0.000000,
    'n':                  241.000000,
    'p':                  0.995170,
    'upper_conf':         np.inf,
    'z':                  -0.006054}


