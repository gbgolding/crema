#!/usr/bin/python

from Bio import SeqIO
from sklearn.svm import SVC
from sklearn.model_selection import GridSearchCV, cross_val_score, KFold
import numpy as np
import csv
from featuresetup_module import transcript_info, transcript_info_dict
from sklearn.externals import joblib
from sklearn import preprocessing
import datetime
from sklearn.ensemble import GradientBoostingClassifier
from sklearn.feature_selection import RFECV
#from collections import Counter
from treeinterpreter import treeinterpreter as ti

hsapiens_info,hsapiens_dict,hsapiens_names = transcript_info_dict("../data/training_files/h_sapiens_random3000.fa", "../data/training_files/h_sapiens_random3000.cpat.txt", "../data/training_files/h_sapiens_random3000.fa.tab")
print("imported human info")

osat_info, osat_dict, osat_names = transcript_info_dict("../data/training_files/o_sativa_random3000.fa", "../data/training_files/o_sativa_random3000.cpat.txt", "../data/training_files/o_sativa_random3000.fa.tab")
print("imported rice info")

lncRNA_info, lncRNA_dict, lncRNA_names = trans_info_dict("../data/training_files/all_lncRNA_nodup.fa","../data/training_files/all_lncRNA_nodup.humantrained.cpat.txt", "../data/training_files/all_lncRNA_nodup.fa.tab")

print("imported lncRNA info")


wanted_keys = ["align_perc_len", "align_perc_ORF", "align_length", "ORF", "length", "GC", "fickett", "hexamer", "identity"] # try minimal number of features!


# this will be able to be removed because we are using all possible keys most likely?
hsapiens_sub = {gene:{feature:hsapiens_dict[gene][feature] for feature in wanted_keys} for gene in hsapiens_names} 

osat_sub = {gene:{feature:osat_dict[gene][feature] for feature in wanted_keys} for gene in osat_names} 
lncRNA_sub = {gene:{feature:lncRNA_dict[gene][feature] for feature in wanted_keys} for gene in lncRNA_names} 


hsapiens_info = np.array([[hsapiens_sub[gene][feature] for feature in sorted(hsapiens_sub[gene])] for gene in hsapiens_names], dtype=float) #to keep order!!!!!!!	y

osat_info = np.array([[osat_sub[gene][feature] for feature in sorted(osat_sub[gene])] for gene in osat_names], dtype=float) #to keep order!!!!!!!	y
lncRNA_info = np.array([[lncRNA_sub[gene][feature] for feature in sorted(lncRNA_sub[gene])] for gene in lncRNA_names], dtype=float) #to keep order!!!!!!!	y

proteins = np.concatenate((hsapiens_info, osat_info), axis=0)


y = []


for num in range(0,proteins.shape[0]):
        y.append(0)
for num in range(0,lncRNA_info.shape[0]):
        y.append(1)

y = np.asarray(y)

X = np.concatenate((proteins, lncRNA_info), axis=0)

print(y.shape)
print(X.shape)

X_normalized = preprocessing.normalize(X, norm='l2')

trans_info, trans_dict, subset_names = transcript_info_dict("ips1.fa", "ips1.cpat.txt", "ips1_blast.tab")
#trans_info, trans_dict, subset_names = transcript_info_dict("cool_cold.fa", "cool_cold.fa.cpat.txt", "cool_cold.fa.tab")

test_sub = {gene:{feature:trans_dict[gene][feature] for feature in wanted_keys} for gene in subset_names} 

test_info = np.array([[test_sub[gene][feature] for feature in sorted(test_sub[gene])] for gene in subset_names], dtype=float) #to keep order!!!!!!!	y

test_info_n = preprocessing.normalize(test_info, norm='l2') # this should be done WITH the !!!!


clf = GradientBoostingClassifier(n_estimators=100, learning_rate= 0.04, subsample=0.6, max_depth=10, random_state=339)
clf.fit(X_normalized, y)

# Uncomment to save classifer as pickle
#joblib.dump(clf, 'model2.pkl')


test_pred = clf.predict(test_info_n)
test_prob = clf.predict_proba(test_info_n)


