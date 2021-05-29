import os
import pandas as pd
from pathos.multiprocessing import ProcessPool
import rpy2

import biomarker_survival as surv
from matplotlib import pyplot
pyplot.switch_backend('agg')

# R is a bit spammy with warnings (e.g. every non-converged cox)
rpy2.robjects.r['options'](warn=-1)

class ZscoreCommon():
  def __init__(self, prep_data, ctype_cleaning, extra_data=None,
               threshold_percent=0.02,
               use_presence_percentage_of_patients_threshold=False):
    # fn that takes the path of the raw data, and returns df
    #  (path) -> df
    self.prep_data = prep_data

    # if data prep requires additional data (e.g. for methylation)
    self.extra_data = extra_data

    # fn that takes the whole df, cancer type string, and clinical data (only needed for mutations
    # and returns an appropriately cleaned
    # (e.g. removal of non-01s)
    # df indexed by patients with site (gene, etc) names in the columns
    #  (df, ctype, ctype_clinical) -> ctype_df
    self.ctype_cleaning = ctype_cleaning

    # Mutation cutoffs work a bit differently from the other types.
    # instead of using total number of observations, we want a cutoff
    # based on the number of mutations in the patient population for a single gene.
    self.use_presence_percentage_of_patients_threshold = use_presence_percentage_of_patients_threshold
    self.threshold_percent = threshold_percent



  def cancer_type_cox_multivar(self, df, ctype_clinical, ctype, outdir, additional_vars):
    df = self.ctype_cleaning(df, ctype, ctype_clinical)
    df = df.join(ctype_clinical, how='inner')
    df = df.join(additional_vars, how='inner')
    num_patients = len(df)

    cox_df = pd.DataFrame()
    cox = pd.Series()
    for gene in df:
      if gene in (['time', 'censor'] + list(additional_vars.columns)):
        continue

      # don't do zscores if there aren't enough non-missing observations
      cox_in_df = df[[gene, 'time', 'censor'] + list(additional_vars.columns)]
      cox_in_df = cox_in_df.dropna(how='any', axis=0)

      if self.use_presence_percentage_of_patients_threshold:
        num_mutations = int(cox_in_df[gene].sum())
        if num_mutations >= self.threshold_percent * num_patients:
          cox = surv.do_multivariate_cox(cox_in_df['time'], cox_in_df['censor'], cox_in_df[gene], additional_vars=additional_vars)
          cox['num mutated patients'] = num_mutations
          cox = pd.DataFrame(cox, index=[gene])
          cox_df = cox_df.append(cox, sort=True)
      elif cox_in_df.shape[0] > 0 and len(cox_in_df[gene]) > 10:
        cox = surv.do_multivariate_cox(cox_in_df['time'], cox_in_df['censor'], cox_in_df[gene], additional_vars=additional_vars)
        cox = pd.DataFrame(cox, index=[gene])
        cox_df = cox_df.append(cox, sort=True)


    outfile = os.path.join(outdir, ctype + '.zscores.out.csv')
    if self.use_presence_percentage_of_patients_threshold:
        outfile = os.path.join(outdir, ctype + '_' + str(self.threshold_percent) + '.zscores.out.csv')
    cox_df.to_csv(outfile, index_label='gene')
    return cox_df

  def cancer_type_cox(self, df, ctype_clinical, ctype, outdir):
    df = self.ctype_cleaning(df, ctype, ctype_clinical)
    df = df.join(ctype_clinical, how='inner')
    num_patients = len(df)

    cox_df = pd.DataFrame()
    cox = pd.Series(dtype='float64')
    for gene in df:
      if gene in ['time', 'censor']:
        continue

      # don't do zscores if there aren't enough non-missing observations
      cox_in_df = df.loc[:,[gene, 'time', 'censor']]
      cox_in_df = cox_in_df.dropna(how='any', axis=0)

      if self.use_presence_percentage_of_patients_threshold:
        # if we're here, this is probably a mutations analysis, so make it easy to remember that.
        num_mutations = int(cox_in_df[gene].sum())
        if num_mutations >= self.threshold_percent * num_patients:
          cox = surv.do_cox(cox_in_df['time'], cox_in_df['censor'], cox_in_df[gene])
          cox['num mutated patients'] = num_mutations
          cox = pd.DataFrame(cox, index=[gene])
          cox_df = cox_df.append(cox, sort=True)
      elif cox_in_df.shape[0] > 0 and len(cox_in_df[gene]) > 10:
        cox = surv.do_cox(cox_in_df['time'], cox_in_df['censor'], cox_in_df[gene])
        cox = pd.DataFrame(cox, index=[gene])
        cox_df = cox_df.append(cox, sort=True)

    outfile = os.path.join(outdir, ctype + '.zscores.out.csv')
    if self.use_presence_percentage_of_patients_threshold:
        outfile = os.path.join(outdir, ctype + '_' + str(self.threshold_percent) + '.zscores.out.csv')
    cox_df.to_csv(outfile, index_label='gene')
    return cox_df

  def multiprocess_cox(self, args):
    return self.cancer_type_cox(*args)

  def multiprocess_multivar_cox(self, args):
    return self.cancer_type_cox_multivar(*args)

  def zscores(self, df_path, clinical, outdir, additional_vars={}, parallel_workers=0):
    '''
    Zscores for all cancer types, maybe in parallel, maybe with multivariates

    Parameters
    ----------
    df_path : str
      Path to raw input data.
    clinical : str
      Path to TCGA CDR clinical data.
    outdir : str
      destination directory for output files.
    additional_vars : dict, optional
      Dictionary of dataframes, keyed by cancer type variables to include in a multivariate analysis. The default is {}.
    parallel_workers : int, optional
      Number of processes to start in the multiprocess pool. The default is 0.

    Returns
    -------
    None.
    '''
    df = self.prep_data(df_path, extra_data=self.extra_data)
    tcga_cdr = surv.TCGA_CDR_util(clinical)

    args = []
    results = []
    for ctype in tcga_cdr.cancer_types():
      ctype_clinical = tcga_cdr.cancer_type_data(ctype)
      if parallel_workers == 0 and ctype not in additional_vars.keys():
        ctype_cox = self.cancer_type_cox(df, ctype_clinical, ctype, outdir)
        results.append(ctype_cox)
      elif parallel_workers == 0 and ctype in additional_vars.keys():
        ctype_cox = self.cancer_type_cox_multivar(df, ctype_clinical, ctype, outdir, additional_vars[ctype])
      elif parallel_workers > 0 and ctype not in additional_vars.keys():
        args.append([df, ctype_clinical, ctype, outdir])
      else:
        args.append([df, ctype_clinical, ctype, outdir, additional_vars[ctype]])

    if parallel_workers > 0:
      p = ProcessPool(nodes=parallel_workers)
      if ctype not in additional_vars.keys():
        results = p.map(self.multiprocess_cox, args)
      else:
        results = p.map(self.multiprocess_multivar_cox, args)

    print('Zscores complete')
    return results

  def metadata(self, df_path, clinical, return_data=False, extra_data=None):
    df = self.prep_data(df_path, extra_data=self.extra_data)
    tcga_cdr = surv.TCGA_CDR_util(clinical)
    metadata_df = {}
    prepped_data = {}
    for ctype in tcga_cdr.cancer_types():
      print(ctype)
      ctype_clinical = tcga_cdr.cancer_type_data(ctype)
      ctype_df = self.ctype_cleaning(df, ctype, ctype_clinical)
      ctype_df = ctype_df.join(ctype_clinical, how='inner')
      metadata_df[ctype] = ctype_df.shape[0]
      prepped_data[ctype] = ctype_df

    if return_data:
      return pd.concat(prepped_data)
    return metadata_df
