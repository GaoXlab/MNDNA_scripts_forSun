## /Users/xingyun/anaconda3/bin/pip3 install pyreadr 
import pyreadr
import pandas as pd
import numpy as np
import json, pickle, argparse
import math, os, random
from collections import Counter
from sklearn.svm import SVC
from sklearn.model_selection import StratifiedKFold, KFold, RepeatedStratifiedKFold, cross_val_score, cross_val_predict, cross_validate
from sklearn.metrics import roc_auc_score
random.seed(1234)
np.random.seed(1234)

filepath = 'modelData/'
trnval_file = pyreadr.read_r(filepath+'Input_raw_nor.RData')
test_file = pyreadr.read_r(filepath+'Input_raw_nor.test.RData')

feature = pd.read_csv('results/DataMatrix_SignedFindMarkers_40regions.log',sep='\t').index.tolist()

trnval = trnval_file['nor'].loc[feature]
test = test_file['nor'].loc[feature]

#### training set
trnval_x = trnval[[i for i in trnval.columns.tolist() if 'A_Apc' in i or 'A_Ctrl' in i]].T
trnval_y = [1 if 'A_Apc' in i else 0 for i in trnval_x.index.tolist()]; trnval_y = pd.DataFrame(trnval_y); trnval_y.index = trnval_x.index.tolist(); trnval_y= trnval_y[0]

#### test set
test_x = test.T
test_x.columns.tolist()==trnval_x.columns.tolist()

# model
svc = SVC(kernel='linear', probability=True, random_state=1234)#gamma='auto', 

n_splits = 10
df_merge = pd.DataFrame()
cv = KFold(n_splits=n_splits, shuffle=True, random_state=1234)
df_train_proba_temp, df_train_pred_temp, df_train_label_temp = [], [], []
df_valid_proba_temp, df_valid_pred_temp, df_valid_label_temp, df_valid_name_temp = [], [], [], []
df_test_proba_temp, df_test_pred_temp = pd.DataFrame(), pd.DataFrame()
round_i = 0
classifier = svc
for i, (train, valid) in enumerate(cv.split(trnval_x)):
    print(valid)
    train_X, valid_X, train_y, valid_y = trnval_x.iloc[train],trnval_x.iloc[valid],trnval_y.iloc[train],trnval_y.iloc[valid]
    classifier.fit(train_X,train_y)
    df_train_label_temp.extend(np.array(train_y))
    df_train_pred_temp.extend(classifier.predict(train_X))
    df_train_proba_temp.extend(classifier.predict_proba(train_X)[:,1])
    df_valid_name_temp.extend(valid_X.index.tolist())
    df_valid_label_temp.extend(np.array(valid_y))
    df_valid_pred_temp.extend(classifier.predict(valid_X))
    df_valid_proba_temp.extend(classifier.predict_proba(valid_X)[:,1])
    df_test_pred_temp['test_pred_%i' % (round_i)] = classifier.predict(test_x)
    df_test_proba_temp['test_pred_%i' % (round_i)] = classifier.predict_proba(test_x)[:,1]
    round_i = round_i+1

df_valid_proba = pd.DataFrame(df_valid_proba_temp); df_valid_proba.index = df_valid_name_temp
df_valid_proba['testset_label'] = 'cross-validation'
df_test_proba_temp.index = test_x.index.tolist()
df_test_proba_mean = pd.DataFrame(df_test_proba_temp.mean(axis=1))

df_valid_proba.to_csv('results/svc_cv_predictivescore.log', sep='\t')
df_test_proba_mean.to_csv('results/svc_test_predictivescore.log', sep='\t')
pd.concat([pd.DataFrame(trnval_x.columns.tolist()), pd.DataFrame(classifier.coef_.tolist()[0])], axis=1).to_csv('results/feature_importance.log',sep='\t')