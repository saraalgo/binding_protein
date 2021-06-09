# THERE IS NOT AUTOMATIZED, MODELS WILL BE IMPLEMENT WITH OBJECT ORIENTED PROGRAMMING!!
# JUST USED TO OBTAIN FAST RESULTS!!

from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
from xgboost import XGBClassifier
from sklearn.neighbors import KNeighborsClassifier
from sklearn.metrics import classification_report, confusion_matrix, accuracy_score
import matplotlib.pyplot as plt 
import pandas as pd

# Load data by hand, neccesary to implement as function!!! + Add Label column
Mono_freq = pd.read_csv('../../data/trial3/descriptors/monomers_descriptors.csv', index_col=0)
output = Mono_freq.index.str.replace('\d+', '')
Mono_freq["target"] = output.values
Di_freq = pd.read_csv('../../data/trial3/descriptors/dimers_descriptors.csv', index_col=0)
Di_freq["target"] = output.values
Tetra_freq = pd.read_csv('../../data/trial3/descriptors/tetramers_descriptors.csv', index_col=0)
Tetra_freq["target"] = output.values

MGW = pd.read_csv('../../data/trial3/descriptors/MGW_descriptors.csv', index_col=0)
MGW["target"] = output.values
HelT = pd.read_csv('../../data/trial3/descriptors/HelT_descriptors.csv', index_col=0)
HelT["target"] = output.values
ProT = pd.read_csv('../../data/trial3/descriptors/ProT_descriptors.csv', index_col=0)
ProT["target"] = output.values
Roll = pd.read_csv('../../data/trial3/descriptors/Roll_descriptors.csv', index_col=0)
Roll["target"] = output.values

MGW_pca = pd.read_csv('../../data/trial3/descriptors/MGW_descriptors_PCA.csv', index_col=0)
MGW_pca["target"] = output.values
HelT_pca = pd.read_csv('../../data/trial3/descriptors/HelT_descriptors_PCA.csv', index_col=0)
HelT_pca["target"] = output.values
ProT_pca = pd.read_csv('../../data/trial3/descriptors/ProT_descriptors_PCA.csv', index_col=0)
ProT_pca["target"] = output.values
Roll_pca = pd.read_csv('../../data/trial3/descriptors/Roll_descriptors_PCA.csv', index_col=0)
Roll_pca["target"] = output.values

MGW_kbest = pd.read_csv('../../data/trial3/descriptors/MGW_descriptors_kbest.csv', index_col=0)
MGW_kbest["target"] = output.values
HelT_kbest = pd.read_csv('../../data/trial3/descriptors/HelT_descriptors_kbest.csv', index_col=0)
HelT_kbest["target"] = output.values
ProT_kbest = pd.read_csv('../../data/trial3/descriptors/ProT_descriptors_kbest.csv', index_col=0)
ProT_kbest["target"] = output.values
Roll_kbest = pd.read_csv('../../data/trial3/descriptors/Roll_descriptors_kbest.csv', index_col=0)
Roll_kbest["target"] = output.values

Tetra_freq_pca = pd.read_csv('../../data/trial3/descriptors/tetramers_descriptors_PCA.csv', index_col=0)
Tetra_freq_pca["target"] = output.values
Tetra_freq_kbest = pd.read_csv('../../data/trial3/descriptors/tetramers_descriptors_kbest.csv', index_col=0)
Tetra_freq_kbest["target"] = output.values

list_features = [Mono_freq,Di_freq,Tetra_freq,MGW,MGW_pca,MGW_kbest,HelT,HelT_pca,HelT_kbest,ProT,ProT_pca,ProT_kbest,Roll,Roll_pca,Roll_kbest,Tetra_freq_pca,Tetra_freq_kbest]   
names_features = ["Mono_freq","Di_freq","Tetra_freq","MGW","MGW_pca","MGW_kbest","HelT","HelT_pca","HelT_kbest","ProT","ProT_pca","ProT_kbest","Roll","Roll_pca","Roll_kbest","Tetra_freq_pca","Tetra_freq_kbest"]

# Random Forest 
accuracy_RF = []

for feature in list_features:
    X_train, X_test, y_train, y_test = train_test_split(feature.loc[:, feature.columns != 'target'],feature.loc[:, feature.columns == 'target'], test_size=0.30, random_state=0)
    RF = RandomForestClassifier(n_estimators=1000,random_state=0)
    RF.fit(X_train, y_train.values.ravel())
    y_pred = RF.predict(X_test)
    accuracy_new = accuracy_score(y_test, y_pred)
    "Accuracy: %.2f%%" % (accuracy_new * 100.0)
    accuracy_RF.append(accuracy_new)

plt.figure(figsize=(15, 15))
plt.scatter(accuracy_RF,names_features)
plt.title('Random Forest accuracy')
plt.ylabel('Number of features')
plt.xlabel('Accuracy')
plt.savefig('../../data/trial3/results/randomforest.png')
plt.show()
plt.close()

## XGboost
accuracy_XGboost = []

for feature in list_features:
    X_train, X_test, y_train, y_test = train_test_split(feature.loc[:, feature.columns != 'target'],feature.loc[:, feature.columns == 'target'], test_size=0.30, random_state=0)
    Xgboost = XGBClassifier()
    Xgboost.fit(X_train, y_train.values.ravel())
    y_pred = Xgboost.predict(X_test)
    accuracy_XGboost.append(accuracy_score(y_test, y_pred))

plt.scatter(accuracy_XGboost,names_features)
plt.title('XGBoost accuracy')
plt.ylabel('Number of features')
plt.xlabel('Accuracy')
plt.savefig('../../data/trial3/results/XGboost.png')
plt.show()
plt.close()


## KNN
k_range =range(1,26)
scores_list = []
accuracy_KNN = []

for feature in list_features:
    X_train, X_test, y_train, y_test = train_test_split(feature.loc[:, feature.columns != 'target'],feature.loc[:, feature.columns == 'target'], test_size=0.30, random_state=0)
    for k in k_range:
        knn=KNeighborsClassifier(n_neighbors=k)
        knn.fit(X_train,y_train.values.ravel())
        y_pred = knn.predict(X_test)
        scores_list.append(accuracy_score(y_test,y_pred))
    accuracy_KNN.append(max(scores_list))

plt.figure(figsize=(15, 15))
plt.scatter(accuracy_KNN,names_features)
plt.title('KNN accuracy')
plt.ylabel('Number of features')
plt.xlabel('Accuracy')
plt.savefig('../../data/trial3/results/KNN.png')
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
sns_plot.savefig("../../data/trial3/results/Models_test.png")
plt.close()
