from setuptools import setup, find_packages

setup(name='comprehensive_tcga_survival',
      version='0.2',
      description='Comprehensive TCGA Survival',
      url='http://github.com/joansmith/comprehensive-tcga-survival',
      author='Joan Smith',
      author_email='joans@alum.mit.edu',
      license='MIT',
      packages=find_packages(),
      zip_safe=False,
      setup_requires=['pytest-runner'],
      install_requires=[
        'numpy~=1.20.3',
        'pandas>=0.25.1',
        'cloudpickle~=1.2.2',
        'dask[dataframe]~=2.9.0',
        'rpy2==3.4.5',
        'anndata',
        'wget',
        'matplotlib',
        'scipy>=1.3.1',
        'intervaltree',
        'scanpy',
        'statsmodels',
        'tzlocal',
        'pathos',
        'openpyxl',
        'biomarker-survival>=0.2.4'
        ],
      tests_require=['pytest', 'pytest-datafiles'])

