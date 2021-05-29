import pytest
import os
import pandas as pd
import biomarker_survival as surv

from comprehensive_tcga_survival import cn
from comprehensive_tcga_survival import ZscoreCommon

FIXTURE_DIR = os.path.join(
   os.path.dirname(os.path.realpath(__file__)),
      'test_data',
        )


def test_prep_annotation_file_structure():
  hgnc_annotation = os.path.join(FIXTURE_DIR,
                                'gencode.v32.annotation.small.gtf')
  annotation = cn.prep_annotation_file(hgnc_annotation)
  assert annotation.index.name == 'gene_name'
  assert annotation['feature_type'].value_counts().shape == (1,)

def test_prep_annotation_file_one_per_gene():
  hgnc_annotation = os.path.join(FIXTURE_DIR,
                                'gencode.v32.annotation.small.gtf')
  annotation = cn.prep_annotation_file(hgnc_annotation)
  # check that there's only one row for each gene
  assert annotation.reset_index().groupby('gene_name').count().max()['start'] == 1
  assert annotation.shape == (977, 5)

def test_prep_annotation_file_minimum():
  hgnc_annotation = os.path.join(FIXTURE_DIR,
                                'gencode.v32.annotation.small.gtf')
  annotation = cn.prep_annotation_file(hgnc_annotation)
  #  check that we successfully took the minimum
  assert annotation.loc['LINC01409']['start'] == 778747

def test_parse_chrom():
  assert cn.parse_chrom('chr1') == 1
  assert cn.parse_chrom('chr10') == 10
  assert cn.parse_chrom('chrX') == 23
  assert cn.parse_chrom('foo') == None
  assert cn.parse_chrom('chrM') == 0

def test_prep_and_clean():
  cn_path = os.path.join(FIXTURE_DIR, 'cn_small.seg')
  hgnc_annotation = os.path.join(FIXTURE_DIR,
                                'gencode.v32.annotation.small.gtf')
  cn_data = cn.prep_data(cn_path, extra_data=hgnc_annotation)

  assert cn_data.shape == (36, 977)
  assert 'TCGA-OR-A5J3-10A-01D-A29K-01' in list(cn_data.index)
  assert 'ZNF683' in cn_data.columns
  assert 'chr' not in cn_data.index
  assert 'chr' not in cn_data.columns

def test_prep_and_clean_with_output(tmpdir):
  cn_path = os.path.join(FIXTURE_DIR, 'cn_small.seg')
  hgnc_annotation = os.path.join(FIXTURE_DIR,
                                'gencode.v32.annotation.small.gtf')
  cn_data = cn.prep_data(cn_path, extra_data=hgnc_annotation, outdir=tmpdir)

  assert cn_data.shape == (38, 977) # includes chromosome and start location
  assert 'TCGA-OR-A5J3-10A-01D-A29K-01' in list(cn_data.index)
  assert 'ZNF683' in cn_data.columns
  assert 'chr' in cn_data.index
  assert len(os.listdir(tmpdir)) == 2
  assert sorted(os.listdir(tmpdir)) == sorted(['cn_by_gene.csv', 'prepped_annotations.csv'])

def test_zscores(tmpdir):
  cn_path = os.path.join(FIXTURE_DIR, 'cn_small.seg')
  hgnc_annotation = os.path.join(FIXTURE_DIR,
                                'gencode.v32.annotation.small.gtf')
  path = os.path.join(FIXTURE_DIR, 'TCGA-CDR-SupplementalTableS1-2019-05-27.xlsx')
  clin = surv.TCGA_CDR_util(path)

  cn_data = cn.prep_data(cn_path, extra_data=hgnc_annotation)
  small_cn_data = cn_data.iloc[:,0:20]

  cn_zscores = ZscoreCommon(cn.prep_data,
                            cn.ctype_cleaning,
                            extra_data=hgnc_annotation)

  cox = cn_zscores.cancer_type_cox(small_cn_data, clin.cancer_type_data('ACC'),
                                            'ACC', str(tmpdir))
  print(cox)
  assert cox.shape == (13, 6)
  assert 'AGTRAP' in list(cox.index)
  assert cox.loc['AGTRAP'].round(6).to_dict() == {
    'hazard_ratio':       1.641397,
    'lower_conf':         0.12899,
    'n':                  18.0,
    'p':                  0.702574,
    'upper_conf':         20.886813,
    'z':                  0.381848}


