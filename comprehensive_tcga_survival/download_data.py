import pandas as pd
import wget
import os
from pathlib import Path
import gzip
import shutil

home_dir = os.path.expanduser('~')
analysis_dir = Path(os.path.join(home_dir, 'comprehensive-tcga-survival-analysis'))
data_dir = Path(os.path.join(analysis_dir, 'raw-data'))
multivariate_dir = Path(os.path.join(analysis_dir, 'Multivariate'))

Path.mkdir(analysis_dir, exist_ok=True, parents=True)
Path.mkdir(data_dir, exist_ok=True, parents=True)
Path.mkdir(multivariate_dir, exist_ok=True, parents=True)

def download_file(url, name, kind=None, ddir=data_dir):
    rna = 'http://api.gdc.cancer.gov/data/3586c0da-64d0-4b74-a449-5ff4d9136611'
    path = os.path.join(ddir, name)
    wget.download(url, path)
    if kind:
      print(kind, 'download complete')

def ungzip(p):
    with gzip.open(p, 'rb') as f_in:
        with open(p[:-3], 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)

download_file(
    'http://api.gdc.cancer.gov/data/3586c0da-64d0-4b74-a449-5ff4d9136611',
    'EBPlusPlusAdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.tsv',
    'rna'
)

download_file(
    'http://api.gdc.cancer.gov/data/fcbb373e-28d4-4818-92f3-601ede3da5e1',
    'TCGA-RPPA-pancan-clean.txt',
    'rppa'
)

download_file(
    'http://api.gdc.cancer.gov/data/1c6174d9-8ffb-466e-b5ee-07b204c15cf8',
    'pancanMiRs_EBadjOnProtocolPlatformWithoutRepsWithUnCorrectMiRs_08_04_16.csv',
)

download_file(
    'http://api.gdc.cancer.gov/data/55d9bf6f-0712-4315-b588-e6f8e295018e',
    'PanCanAtlas_miRNA_sample_information_list.txt',
    'mirna'
)

download_file(
    'http://api.gdc.cancer.gov/data/00a32f7a-c85f-4f86-850d-be53973cbc4d',
    'broad.mit.edu_PANCAN_Genome_Wide_SNP_6_whitelisted.seg',
    'cn'
)

download_file(
    'http://api.gdc.cancer.gov/data/1c8cfe5f-e52d-41ba-94da-f15ea1337efc',
    'mc3.v0.2.8.PUBLIC.maf.gz',
)
ungzip(os.path.join(data_dir, 'mc3.v0.2.8.PUBLIC.maf.gz'))
print('mutations download complete')

download_file(
    'https://api.gdc.cancer.gov/data/1b5f413e-a8d1-4d10-92eb-7c4ae739ed81',
    'TCGA-CDR-SupplementalTableS1.xlsx',
    'clinical data'
)

download_file(
    'https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_32/gencode.v32.annotation.gtf.gz',
    'gencode.v32.annotation.gtf.gz'
)
ungzip(os.path.join(data_dir, 'gencode.v32.annotation.gtf.gz'))
print('gencode annotation download complete')


download_file(
    'http://api.gdc.cancer.gov/data/d82e2c44-89eb-43d9-b6d3-712732bf6a53',
    'jhu-usc.edu_PANCAN_merged_HumanMethylation27_HumanMethylation450.betaValue_whitelisted.tsv'
)

download_file(
    'https://webdata.illumina.com/downloads/productfiles/humanmethylation450/Humanmethylation450_15017482_v1-2.csv',
    'HumanMethylation450_15017482_v1-2.csv',
    'methylation'
)

download_file(
  'https://storage.googleapis.com/cancer-survival-km-2021/Grade_key.xlsx',
  'Grade_key.xlsx',
  'grade key',
  ddir=multivariate_dir)

download_file(
  'https://storage.googleapis.com/cancer-survival-km-2021/Grade_key.xlsx',
  'Grade_key.xlsx',
  'grade key',
  ddir=multivariate_dir
)

download_file(
  'https://storage.googleapis.com/cancer-survival-km-2021/Stage_key.xlsx',
  'Stage_key.xlsx',
  'Stage key',
  ddir=multivariate_dir
)

