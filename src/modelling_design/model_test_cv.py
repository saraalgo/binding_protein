# THERE IS NOT AUTOMATIZED, MODELS WILL BE IMPLEMENT WITH OBJECT ORIENTED PROGRAMMING!!
# JUST USED TO OBTAIN FAST RESULTS!!
import os
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'
""" 
0 = all messages are logged (default behavior)
1 = INFO messages are not printed
2 = INFO and WARNING messages are not printed
3 = INFO, WARNING, and ERROR messages are not printed 
"""

from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
from xgboost import XGBClassifier
from sklearn.neighbors import KNeighborsClassifier
from sklearn.metrics import classification_report, confusion_matrix, accuracy_score
from sklearn.model_selection import cross_val_score
from sklearn.model_selection import RepeatedStratifiedKFold
import matplotlib.pyplot as plt 
import pandas as pd
import numpy as np

# Load data by hand, neccesary to implement as function!!! + Add Label column
Mono_freq = pd.read_csv('../../data/trial2/descriptors/monomers_descriptors.csv', index_col=0)
output = Mono_freq.index.str.replace('\d+', '')
Mono_freq["target"] = output.values
Di_freq = pd.read_csv('../../data/trial2/descriptors/dimers_descriptors.csv', index_col=0)
Di_freq["target"] = output.values
Tetra_freq = pd.read_csv('../../data/trial2/descriptors/tetramers_descriptors.csv', index_col=0)
Tetra_freq["target"] = output.values

MGW = pd.read_csv('../../data/trial2/descriptors/MGW_descriptors.csv', index_col=0)
MGW["target"] = output.values
HelT = pd.read_csv('../../data/trial2/descriptors/HelT_descriptors.csv', index_col=0)
HelT["target"] = output.values
ProT = pd.read_csv('../../data/trial2/descriptors/ProT_descriptors.csv', index_col=0)
ProT["target"] = output.values
Roll = pd.read_csv('../../data/trial2/descriptors/Roll_descriptors.csv', index_col=0)
Roll["target"] = output.values

MGW_pca = pd.read_csv('../../data/trial2/descriptors/MGW_descriptors_PCA.csv', index_col=0)
MGW_pca["target"] = output.values
HelT_pca = pd.read_csv('../../data/trial2/descriptors/HelT_descriptors_PCA.csv', index_col=0)
HelT_pca["target"] = output.values
ProT_pca = pd.read_csv('../../data/trial2/descriptors/ProT_descriptors_PCA.csv', index_col=0)
ProT_pca["target"] = output.values
Roll_pca = pd.read_csv('../../data/trial2/descriptors/Roll_descriptors_PCA.csv', index_col=0)
Roll_pca["target"] = output.values

MGW_kbest = pd.read_csv('../../data/trial2/descriptors/MGW_descriptors_kbest.csv', index_col=0)
MGW_kbest["target"] = output.values
HelT_kbest = pd.read_csv('../../data/trial2/descriptors/HelT_descriptors_kbest.csv', index_col=0)
HelT_kbest["target"] = output.values
ProT_kbest = pd.read_csv('../../data/trial2/descriptors/ProT_descriptors_kbest.csv', index_col=0)
ProT_kbest["target"] = output.values
Roll_kbest = pd.read_csv('../../data/trial2/descriptors/Roll_descriptors_kbest.csv', index_col=0)
Roll_kbest["target"] = output.values

Tetra_freq_pca = pd.read_csv('../../data/trial2/descriptors/tetramers_descriptors_PCA.csv', index_col=0)
Tetra_freq_pca["target"] = output.values
Tetra_freq_kbest = pd.read_csv('../../data/trial2/descriptors/tetramers_descriptors_kbest.csv', index_col=0)
Tetra_freq_kbest["target"] = output.values

list_features = [Mono_freq,Di_freq,Tetra_freq,MGW,MGW_pca,MGW_kbest,HelT,HelT_pca,HelT_kbest,ProT,ProT_pca,ProT_kbest,Roll,Roll_pca,Roll_kbest,Tetra_freq_pca,Tetra_freq_kbest]   
names_features = ["Mono_freq","Di_freq","Tetra_freq","MGW","MGW_pca","MGW_kbest","HelT","HelT_pca","HelT_kbest","ProT","ProT_pca","ProT_kbest","Roll","Roll_pca","Roll_kbest","Tetra_freq_pca","Tetra_freq_kbest"]


## Random Forest 
accuracy_RF_CV = []

for idx, feature in enumerate(list_features):
    accuracy_cv = []
    for cv in range(1,11):
        pos = feature[0:879]
        neg = feature.loc[feature.index.str.endswith("."+str(cv))]
        data = pos.append(neg)
        X_train, X_test, y_train, y_test = train_test_split(data.loc[:, data.columns != 'target'],data.loc[:, data.columns == 'target'], test_size=0.30, random_state=0)
        RF = RandomForestClassifier(n_estimators=1000,random_state=0)
        RF.fit(X_train, y_train.values.ravel())
        y_pred = RF.predict(X_test)
        accuracy_new = accuracy_score(y_test, y_pred)
        "Accuracy: %.2f%%" % (accuracy_new * 100.0)
        accuracy_cv.append(accuracy_new)
    print("Mean accuracy for {} CV: {}".format(names_features[idx], (np.mean(accuracy_cv) * 100.0)))
    accuracy_RF_CV.append(np.mean(accuracy_cv))

plt.figure(figsize=(15, 15))
plt.scatter(accuracy_RF_CV,names_features)
plt.title('Random Forest accuracy')
plt.ylabel('Number of features')
plt.xlabel('Accuracy')
plt.savefig('../../data/trial2/results/randomforest_CV.png')
plt.show()
plt.close()

## Random Forest RepeatedStratifiedKFold
aucroc_RF_cv = []

for idx, feature in enumerate(list_features):
    RF = RandomForestClassifier(n_estimators=1000,random_state=0)
    cv = RepeatedStratifiedKFold(n_splits=10, n_repeats=3, random_state=1)
    # evaluate model
    scores = cross_val_score(RF, feature.loc[:, feature.columns != 'target'], feature.loc[:, feature.columns == 'target'].values.ravel(), scoring='roc_auc', cv=cv)
    # summarize performance
    print('Mean ROC AUC para las features de {}: {:.2f}'.format(names_features[idx], np.mean(scores)))
    aucroc_RF_cv.append(np.mean(scores))

plt.figure(figsize=(15, 15))
plt.scatter(aucroc_RF_cv,names_features)
plt.title('Random Forest accuracy')
plt.ylabel('Number of features')
plt.xlabel('Accuracy')
plt.savefig('../../data/trial3/results/aucroc_RF_cv.png')
plt.show()
plt.close()

# ## XGboost
# accuracy_XGboost_CV = []

# for idx, feature in enumerate(list_features):
#     accuracy_cv = []
#     for cv in range(1,11):
#         pos = feature[0:879]
#         neg = feature.loc[feature.index.str.endswith("."+str(cv))]
#         data = pos.append(neg)
#         X_train, X_test, y_train, y_test = train_test_split(data.loc[:, data.columns != 'target'],data.loc[:, data.columns == 'target'], test_size=0.30, random_state=0)
#         Xgboost = XGBClassifier(verbosity=0)
#         Xgboost.fit(X_train, y_train.values.ravel())
#         y_pred = Xgboost.predict(X_test)
#         accuracy_cv.append(accuracy_score(y_test, y_pred))
#     print("Mean accuracy for {} CV: {}".format(names_features[idx], (np.mean(accuracy_cv) * 100.0)))
#     accuracy_XGboost_CV.append(np.mean(accuracy_cv))

# plt.scatter(accuracy_XGboost_CV,names_features)
# plt.title('XGBoost accuracy')
# plt.ylabel('Number of features')
# plt.xlabel('Accuracy')
# plt.savefig('../../data/trial2/results/XGboost_CV.png')
# plt.show()
# plt.close()

## XGboost RepeatedStratifiedKFold
aucroc_XGboost_cv = []

for idx, feature in enumerate(list_features):
    Xgboost = XGBClassifier(scale_pos_weight=10) # Having into account that there are 10 times more negative data
    cv = RepeatedStratifiedKFold(n_splits=10, n_repeats=3, random_state=1)
    # evaluate model
    scores = cross_val_score(Xgboost, feature.loc[:, feature.columns != 'target'], feature.loc[:, feature.columns == 'target'], scoring='roc_auc', cv=cv)
    # summarize performance
    print('Mean ROC AUC para las features de {}: {:.2f}'.format(names_features[idx], np.mean(scores)))
    aucroc_XGboost_cv.append(np.mean(scores))


plt.scatter(aucroc_XGboost_cv,names_features)
plt.title('XGBoost accuracy')
plt.ylabel('Number of features')
plt.xlabel('Accuracy')
plt.savefig('../../data/trial3/results/aucroc_XGboost_cv.png')
plt.show()
plt.close()


## KNN
k_range =range(1,26)
scores_list = []
accuracy_KNN = []

for idx, feature in enumerate(list_features):
    accuracy_cv = []
    for cv in range(1,11):
        pos = feature[0:879]
        neg = feature.loc[feature.index.str.endswith("."+str(cv))]
        data = pos.append(neg)
        X_train, X_test, y_train, y_test = train_test_split(data.loc[:, data.columns != 'target'],data.loc[:, data.columns == 'target'], test_size=0.30, random_state=0)
        for k in k_range:
            knn=KNeighborsClassifier(n_neighbors=k)
            knn.fit(X_train,y_train.values.ravel())
            y_pred = knn.predict(X_test)
            scores_list.append(accuracy_score(y_test,y_pred))
        accuracy_cv.append(max(scores_list))
    print("Mean accuracy for {} CV: {}".format(names_features[idx], (np.mean(accuracy_cv) * 100.0)))
    accuracy_KNN.append(np.mean(accuracy_cv))

plt.figure(figsize=(15, 15))
plt.scatter(accuracy_KNN,names_features)
plt.title('KNN accuracy')
plt.ylabel('Number of features')
plt.xlabel('Accuracy')
plt.savefig('../../data/trial2/results/KNN_CV.png')
plt.show()
plt.close()



# ML plot summary for all the features and all the models
Models_test = pd.DataFrame()
Models_test["Accuracy"] = accuracy_RF+accuracy_XGboost+accuracy_KNN
Models_test["Number of features"] = names_features*3
Models_test["ML Algorithms"] = ["RandomForest"]*17+["XGboost"]*17+["KNN"]*17
Models_test


import seaborn as sns
sns.set_theme(style="darkgrid")

sns.set(rc={'figure.figsize':(15,15)})
sns_plot = sns.relplot(x="Number of features",y="Accuracy",hue="ML Algorithms",data=Models_test)
sns_plot.set_xticklabels(rotation=90)
sns_plot.savefig("../../data/trial2/results/Models_test_CV.png")
plt.close()
