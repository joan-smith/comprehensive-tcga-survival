import pytest
import os
import pandas as pd

import biomarker_survival as surv
from comprehensive_tcga_survival import mutations

from comprehensive_tcga_survival import ZscoreCommon

FIXTURE_DIR = os.path.join(
   os.path.dirname(os.path.realpath(__file__)),
      'test_data',
        )


def test_cancer_type_cox(tmpdir):
  path = os.path.join(FIXTURE_DIR, 'TCGA-CDR-SupplementalTableS1-2019-05-27.xlsx')
  clin = surv.TCGA_CDR_util(path)

  muts_path = os.path.join(FIXTURE_DIR, 'OV_mutations.txt')
  muts_data = mutations.prep_data(muts_path)

  mutation_zscores = ZscoreCommon(mutations.prep_data, mutations.ctype_cleaning,
                                  use_presence_percentage_of_patients_threshold=True,
                                  threshold_percent=0.02)


  cox = mutation_zscores.cancer_type_cox(muts_data, clin.cancer_type_data('OV'), 'OV', str(tmpdir))

  assert len(tmpdir.listdir()) == 1
  assert tmpdir.listdir()[0].basename == 'OV_0.02.zscores.out.csv'


  reference_results = pd.read_csv(os.path.join(FIXTURE_DIR, 'OV_0.02.zscores.out.csv'), index_col=0).astype(float)

  cox_results = pd.read_csv(tmpdir.listdir()[0], index_col=0).astype(float)
  cox_results = cox_results[reference_results.columns]
  pd.testing.assert_frame_equal(cox_results, reference_results, check_less_precise=5)

