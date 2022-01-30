import os
from rpy2.robjects import r
import rpy2
import numpy as np
import rpy2.rinterface as ri

import pytest

# Explicit test to make sure R, rpy2, and survival are set up and functioning well.
# Basically duplicated from biomarker-survival
def test_setup():
  print(os.environ['R_HOME'])

  rpy2.robjects.packages.importr('survival')
  r.assign('time',rpy2.robjects.IntVector(np.array([1,2,3])))
  r.assign('censor', rpy2.robjects.IntVector(np.array([0,1,0])))
  r.assign('split', rpy2.robjects.FloatVector(np.array([1,0,1])))

  coxuh_output = r('summary( coxph(formula = Surv(time, censor) ~ split, model=FALSE, x=FALSE, y=FALSE))')

  print(coxuh_output)
  coef_ind = list(coxuh_output.names).index('coefficients')
  coeffs = coxuh_output[coef_ind]

  patient_count_ind = list(coxuh_output.names).index('n')
  patient_count = coxuh_output[patient_count_ind][0]
  zscore = get_zscore('split', coeffs)
  pvalue = get_pvalue('split', coeffs)

  conf_int_ind = list(coxuh_output.names).index('conf.int')
  conf_int = coxuh_output[conf_int_ind]
  hazard_ratio = conf_int.rx('split', 'exp(coef)')[0]
  lower_conf = conf_int.rx('split', 'lower .95')[0]
  upper_conf = conf_int.rx('split', 'upper .95')[0]

  cox_dict = {
      'n': patient_count,
      'z': zscore,
      'p': pvalue,
      'hazard_ratio': hazard_ratio,
      'lower_conf': lower_conf,
      'upper_conf': upper_conf,
      }
  cox_dict == {'n': 3, 'z': pytest.approx(-0.0005275274394041736), 'p': pytest.approx(0.9995790940202215),
               'hazard_ratio': pytest.approx(6.190130270042719e-10), 'lower_conf': 0.0, 'upper_conf': np.inf}

def get_zscore(variate, coeffs):
  zscore = coeffs.rx(variate, 'z')[0]
  if type(zscore) == rpy2.rinterface_lib.sexp.NARealType:
    zscore = np.nan
  return zscore

def get_pvalue(variate, coeffs):
  pvalue = coeffs.rx(variate, 'Pr(>|z|)')[0]
  if type(pvalue) == rpy2.rinterface_lib.sexp.NARealType:
    pvalue = np.nan
  return pvalue

def get_hazards(variate, cox_output):
  conf_int_ind = list(cox_output.names).index('conf.int')
  conf_int = cox_output[conf_int_ind]
  hazard_ratio = conf_int.rx(variate, 'exp(coef)')[0]
  lower_conf = conf_int.rx(variate, 'lower .95')[0]
  upper_conf = conf_int.rx(variate, 'upper .95')[0]

  return {
    variate + '_hazard_ratio': hazard_ratio,
    variate + '_lower_conf': lower_conf,
    variate + '_upper_conf': upper_conf,
  }


