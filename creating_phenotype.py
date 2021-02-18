# -*- coding: utf-8 -*-
"""
Created on Sun May 10 20:25:17 2020

@author: HO18971
"""

import pandas as pd
from utilities import load_task, plot_olr
df_threshold = pd.read_csv('thresholds_olr.csv').set_index('thresholds')

# ordered_LR_prediction
def f_associate(x):
    if x.gender == 1:
        sex = 'female'
    else:
        sex = 'male'
    if x.age <= df_threshold.loc['thr_1'][sex]:
        return 0
    elif x.age <= df_threshold.loc['thr_2'][sex]:
        return 1
    elif x.age <= df_threshold.loc['thr_3'][sex]:
        return 2
    elif x.age <= df_threshold.loc['thr_4'][sex]:
        return 3
    else:
        return 4

df_task = load_task('phenotype.csv')
df_task['ordered_LR_prediction'] = df_task[['age', 'gender']].apply(f_associate, axis=1)
df_task['delta_grading_olr'] = df_task['grading'] - df_task['ordered_LR_prediction']

# task 0/1 from grading
def f_(t):
    if t in grading_positive:
        return 1
    elif t in grading_negative:
        return 0
    else:
        return 'none'

grading_positive = [1, 2, 3, 4, 5]
grading_negative = [-1, -2, -3, -4, 5]
df_task['delta_grading'] = df_task['delta_grading_olr'].apply(f_)
df_task.to_csv('age_sex_adjusted.csv')

plot_olr(df_task, title='application_olr_model')

