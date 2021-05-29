#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Dec 19 21:51:54 2020

@author: Joan Smith
"""

import pandas as pd
import argparse
import sys
import os


class StageGradeParser():
  def __init__(self, key_file, clin):
    self.key_file = key_file
    self.ctype_key_dict = pd.read_excel(self.key_file, sheet_name=None, header=None, engine='openpyxl')
    self.clin = clin

  def generate_comparison_dict(self, ctype):
    '''
    The labeled/grouped stage/grades are cumulative.
    If an input sheet has N rows, there will be N-1 output columns.
    If a patient is e.g. stage 3, and the groupings are stage <2 and stage >=2, the patient will have a 1 for stage >=2.

    If a patient is stage 4, and the groupings are stage <2, stage >=3, stage >=4, the patient will have a 1 for stage >=3, and stage>=4
    '''
    comparison_dict = {}
    key = self.ctype_key_dict[ctype].reindex(index=self.ctype_key_dict[ctype].index[::-1])
    prev = []
    for i, g in key.iterrows():
      if i < 2:
        continue
      curr = g.dropna(how='all')
      prev = prev + list(curr)
      comparison_dict['GTE_' + g[0].replace(' ', '_')] = prev
    return comparison_dict

  def parse_for_ctype(self, ctype):
    if ctype not in self.ctype_key_dict:
      return pd.DataFrame()

    comparison_dict = self.generate_comparison_dict(ctype)

    data_column_name = self.ctype_key_dict[ctype].iloc[0,0]
    df = self.clin.cancer_type_data(ctype, [data_column_name])
    col = df[data_column_name].str.upper()
    comparison_dict = {k: [i.upper() for i in v] for (k,v) in comparison_dict.items()}

    stage_grade_df = pd.DataFrame()
    for k in comparison_dict.keys():
      stage_grade_df[k] = col.apply(lambda x: x in comparison_dict[k]).astype(int)

    return stage_grade_df
