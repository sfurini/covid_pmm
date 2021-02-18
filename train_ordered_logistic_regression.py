# -*- coding: utf-8 -*-
"""
Created on Sun May 10 20:25:17 2020

@author: HO18971
"""

from mord import LogisticAT
from utilities import load_task, plot_olr
import pandas as pd

df_task = load_task('phenotype.csv') # CHANGE TO THE NAME OF YOUR PHENOTYPE FILE

model_ordinal_m = LogisticAT(alpha=0)
df_task_original_m = df_task[df_task['gender']==0]
model_ordinal_m.fit(df_task_original_m[['age']].astype(int), df_task_original_m['grading'].astype(int))
y_pred_m = model_ordinal_m.predict(df_task_original_m[['age']])
df_task.loc[df_task_original_m.index, 'ordered_LR_prediction'] = y_pred_m

model_ordinal_f = LogisticAT(alpha=0)
df_task_original_f = df_task[df_task['gender']==1]
model_ordinal_f.fit(df_task_original_f[['age']].astype(int), df_task_original_f['grading'].astype(int))
y_pred_f = model_ordinal_f.predict(df_task_original_f[['age']])
df_task.loc[df_task_original_f.index, 'ordered_LR_prediction'] = y_pred_f

thresholds_m = model_ordinal_m.theta_/model_ordinal_m.coef_
thresholds_f = model_ordinal_f.theta_/model_ordinal_f.coef_

df_threshold = pd.DataFrame({'male':thresholds_m, 'female':thresholds_f}, index=['thr_1', 'thr_2', 'thr_3', 'thr_4', 'thr_5'])
df_threshold.index.name = 'thresholds'
df_threshold.to_csv('thresholds_olr.csv')

plot_olr(df_task, title='fitted_olr')
