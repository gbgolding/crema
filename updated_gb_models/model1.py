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
import featuresetup_module


hsapiens_info,hsapiens_dict,hsapiens_names = transcript_info_dict("../data/training_files/h_sapiens_random3000.fa", "../data/training_files/h_sapiens_random3000.cpat.txt", "../data/training_files/h_sapiens_random3000.fa.tab")
print("imported human info")
mmus_info, mmus_dict, mmus_names = transcript_info_dict("../data/training_files/m_musculus_random1000.fa", "../data/training_files/m_musculus_random1000.cpat.txt", "../data/training_files/m_musculus_random1000.fa.tab")
print("imported mouse info")
osat_info, osat_dict, osat_names = transcript_info_dict("../data/training_files/o_sativa_random3000.fa", "../data/training_files/o_sativa_random3000.cpat.txt", "../data/training_files/o_sativa_random3000.fa.tab")
print("imported rice info")
#lncRNA_info, lncRNA_dict, lncRNA_names = transcript_info_dict("fasta_files/subset_lncrna_V2.fa","fasta_files/subset_lncrna_v2.cpat.txt", "fasta_files/subset_lncrna_V2.fa.tab")
#lncRNA_info, lncRNA_dict, lncRNA_names = trans_info_dict_cc("fasta_files/all_lncRNA_nodup.fa","fasta_files/all_lncRNA_nodup.humantrained.cpat.txt", "fasta_files/all_lncRNA_nodup.fa.tab")
lncRNA_info, lncRNA_dict, lncRNA_names = transcript_info_dict("../data/training_files/all_lncRNA_nodup.fa","../data/training_files/all_lncRNA_nodup.humantrained.cpat.txt", "../data/training_files/all_lncRNA_nodup.fa.tab")
print("imported lncRNA info")


wanted_keys = ["align_perc_len", "align_perc_ORF", "align_length", "ORF", "length", "GC", "fickett", "hexamer", "identity"] # try minimal number of features!


hsapiens_sub = {gene:{feature:hsapiens_dict[gene][feature] for feature in wanted_keys} for gene in hsapiens_names} 
mmus_sub = {gene:{feature:mmus_dict[gene][feature] for feature in wanted_keys} for gene in mmus_names} 
osat_sub = {gene:{feature:osat_dict[gene][feature] for feature in wanted_keys} for gene in osat_names} 
lncRNA_sub = {gene:{feature:lncRNA_dict[gene][feature] for feature in wanted_keys} for gene in lncRNA_names} 


hsapiens_info = np.array([[hsapiens_sub[gene][feature] for feature in sorted(hsapiens_sub[gene])] for gene in hsapiens_names], dtype=float) #to keep order!!!!!!!	y
mmus_info = np.array([[mmus_sub[gene][feature] for feature in sorted(mmus_sub[gene])] for gene in mmus_names], dtype=float) #to keep order!!!!!!!	y
osat_info = np.array([[osat_sub[gene][feature] for feature in sorted(osat_sub[gene])] for gene in osat_names], dtype=float) #to keep order!!!!!!!	y
lncRNA_info = np.array([[lncRNA_sub[gene][feature] for feature in sorted(lncRNA_sub[gene])] for gene in lncRNA_names], dtype=float) #to keep order!!!!!!!	y

proteins = np.concatenate((hsapiens_info, mmus_info, osat_info), axis=0)


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


clf = GradientBoostingClassifier(n_estimators=100, learning_rate= 0.04, subsample=0.6, max_depth=10, random_state=339)
clf.fit(X_normalized, y)

## Uncomment if you want to save this model as a pickle

joblib.dump(clf, 'model1.pkl')
