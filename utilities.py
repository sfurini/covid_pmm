# -*- coding: utf-8 -*-
"""
Created on Sun May 10 20:25:17 2020

@author: HO18971
"""

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib as mpl
import matplotlib.pyplot as pl
from sklearn.metrics import confusion_matrix
from sklearn.model_selection import StratifiedKFold
from sklearn.model_selection import cross_val_predict
from sklearn.model_selection import cross_validate
import copy

mpl.rcParams['xtick.labelsize'] = 8
mpl.rcParams['ytick.labelsize'] = 8
mpl.rcParams['lines.linewidth'] = 1
mpl.rcParams['lines.markersize'] = 4
mpl.rcParams['font.size'] = 7
mpl.rcParams['axes.linewidth'] = 0.8
mpl.rcParams['axes.labelsize'] = 8


def metrics(model, X_input, y_input, fold, sample_weight=None, fit_params=None, random_state=0, n_jobs=4):
    list_score = ['accuracy', 'precision', 'recall', 'balanced_accuracy', 'roc_auc'] 
    cv = StratifiedKFold(n_splits=fold, shuffle=True, random_state=random_state)
    df_cv = pd.DataFrame(cross_validate(model, X_input, y_input, cv=cv, scoring=list_score, fit_params=fit_params))[['test_'+tt for tt in list_score]]
    y_pred = cross_val_predict(model, X_input, y_input, cv=cv, fit_params=fit_params, n_jobs=n_jobs)
    cm = confusion_matrix(y_input, y_pred)
    cm_df = pd.DataFrame(cm, index=['Negative', 'Positive'], columns=['Negative (predicted)', 'Positive (predicted)'])
    pl.figure(str(model.__class__) + '_confusion_matrix')
    sns.set(font_scale=2)
    sns.heatmap(cm_df, center=0, cmap=sns.diverging_palette(220, 15, as_cmap=True), annot=True, fmt='g')
    pl.xticks(rotation=0)
    pl.yticks(rotation=90, va="center")
    mng = pl.get_current_fig_manager()
    mng.window.showMaximized()
    pl.show()
    df_overall = df_cv.mean().to_frame(model.__class__)
    return df_cv, df_overall

def load_task(filename):
    #----------------------------------------- load df_task
    df_task = pd.read_csv(filename).set_index('sample')[['age', 'gender', 'grading']]    

    assert np.isin(df_task['grading'].unique(), [0, 1, 2, 3, 4, 5]).all()
    assert (df_task['age'] >= 15).all()
    assert np.isin(df_task['gender'].unique(), [0, 1]).all()
    return df_task

def load_datasets():
    df_original_rare_al1 = pd.read_csv('data_al1_rare.csv').set_index('Unnamed: 0')
    col_name = []
    for gene_ in df_original_rare_al1.columns:
        col_name.append(gene_ + '_rare_al1')
    df_original_rare_al1.columns = col_name
    df_original_rare_al2 = pd.read_csv('data_al2_rare.csv').set_index('Unnamed: 0')
    col_name = []
    for gene_ in df_original_rare_al2.columns:
        col_name.append(gene_ + '_rare_al2')
    df_original_rare_al2.columns = col_name  
    df_original_gc = pd.read_csv('data_gc_unique_hetero.csv').set_index('Unnamed: 0')
    df_original_gc_homo = pd.read_csv('data_gc_unique_homo.csv').set_index('Unnamed: 0')
    col_name = []
    for gene_ in df_original_gc_homo.columns:
        col_name.append(gene_ + '_homo')
    df_original_gc_homo.columns = col_name
    df = pd.concat((df_original_rare_al1, df_original_rare_al2, df_original_gc, df_original_gc_homo), 1)
    return df

def load_genes(sex):
    df_gene_al1 = pd.read_excel('genes_al1_' + sex + '.xlsx')
    df_gene_al2 = pd.read_excel('genes_al2_' + sex + '.xlsx')
    df_gene_gc = pd.read_excel('genes_GC_' + sex + '.xlsx')
    df_gene_gc_homo = pd.read_excel('genes_GC_HOMO_' + sex + '.xlsx')
    df_gene = pd.concat((df_gene_al1, df_gene_al2, df_gene_gc, df_gene_gc_homo), 0)
    if sex == 'male':
        df_gene_al1_x = pd.read_excel('genes_al1_' + 'X_' + sex + '.xlsx')
        df_gene_gc_x = pd.read_excel('genes_GC_' + 'X_' + sex + '.xlsx')
        df_gene = pd.concat((df_gene, df_gene_al1_x, df_gene_gc_x), 0)        
    return df_gene


def plot_olr(df_task, title):
    
    pl.figure(title)
    df_plot_m = df_task[df_task['gender']==0]
    pl.subplot(2,1,1)
    pl.title('male')
    df_plot_m['Patient distribution'] = 'actual distribution'
    df_task_original_copy = copy.deepcopy(df_plot_m)
    df_task_original_copy['grading'] = df_plot_m['ordered_LR_prediction']
    df_task_original_copy['Patient distribution'] = 'predicted distribution'
    df_task_ = pd.concat([df_plot_m[['age', 'grading', 'Patient distribution']], df_task_original_copy[['age', 'grading', 'Patient distribution']]], 0)
    sns.violinplot(x='age', y='grading', data=df_task_, orient='h', inner=None, hue='Patient distribution', split=True,
                   order = [5, 4, 3, 2, 1, 0], color='gray', scale='count')
    pl.plot(df_plot_m['age'], 5-df_plot_m['grading'], 'or')
    pl.plot(df_plot_m['age'], 5-df_plot_m['ordered_LR_prediction'], 'ok', label='predicted grading (OLR)')
    ok = df_plot_m[(df_task['grading'] - df_task['ordered_LR_prediction'])<0]
    pl.plot(ok['age'], 5-ok['grading'], 'og')
    pl.xlabel('')
    pl.ylabel('grading')
    pl.legend()
    
    df_plot_f = df_task[df_task['gender']==1]
    pl.subplot(2, 1, 2)
    pl.title('female')
    df_plot_f['Patient distribution'] = 'actual distribution'
    df_task_original_copy = copy.deepcopy(df_plot_f)
    df_task_original_copy['grading'] = df_plot_f['ordered_LR_prediction']
    df_task_original_copy['Patient distribution'] = 'predicted distribution'
    df_task_ = pd.concat([df_plot_f[['age', 'grading', 'Patient distribution']], df_task_original_copy[['age', 'grading', 'Patient distribution']]], 0)
    sns.violinplot(x='age', y='grading', data=df_task_, orient='h', inner=None, hue='Patient distribution', split=True,
                   order = [5, 4, 3, 2, 1, 0], color='gray', scale='count')
    pl.plot(df_plot_f['age'], 5-df_plot_f['grading'], 'or', label='actual')
    pl.plot(df_plot_f['age'], 5-df_plot_f['ordered_LR_prediction'], 'ok', label='predicted (OLR)')
    ok = df_plot_f[(df_task['grading'] - df_task['ordered_LR_prediction'])<0]
    pl.plot(ok['age'], 5-ok['grading'], 'og', label='actual')
    pl.ylabel('grading')
    pl.xlabel('age')
    pl.legend('', frameon=False)
    pl.savefig(title + '.pdf')
    pl.show()
    pl.close()
    
