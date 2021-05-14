# THERE IS NOT AUTOMATIZED, MODELS WILL BE IMPLEMENT WITH OBJECT ORIENTED PROGRAMMING!!
# JUST USED TO OBTAIN FAST RESULTS!!

from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
from xgboost import XGBClassifier
from sklearn.neighbors import KNeighborsClassifier
from sklearn.metrics import classification_report, confusion_matrix, accuracy_score
import matplotlib.pyplot as plt 
import pandas as pd

# Bring data fom data_preprocessing file
exec(compile(open('data_preprocessing.py').read(), 'data_preprocessing.py', 'exec'))

# Split data in dataframes to work with easily
MGW = d_hist_features['../../data/trial1/hist_features_MGW.csv']
HelT = d_hist_features['../../data/trial1/hist_features_HelT.csv']
ProT = d_hist_features['../../data/trial1/hist_features_ProT.csv']
Roll = d_hist_features['../../data/trial1/hist_features_Roll.csv']

Mono_freq = d_descrip_featu['../../data/trial1/dataSequence_1.csv']
Di_freq = d_descrip_featu['../../data/trial1/dataSequence_2.csv']
Tetra_freq = d_descrip_featu['../../data/trial1/dataSequence_4.csv']

list_features = [MGW,HelT,ProT,Roll,Mono_freq,Di_freq,Tetra_freq]
names_features = ["MGW","HelT","ProT","Roll","Mono_freq","Di_freq","Tetra_freq"]

# Load data by hand, neccesary to implement as function!!!
MGW_kruskal = pd.read_csv('../../data/trial2/FS/MGW_kruskal.csv')
MGW_kruskal = MGW_kruskal.drop(MGW_kruskal.columns[0], axis=1)
HelT_kruskal = pd.read_csv('../../data/trial2/FS/HelT_kruskal.csv')
HelT_kruskal = HelT_kruskal.drop(HelT_kruskal.columns[0], axis=1)
ProT_kruskal = pd.read_csv('../../data/trial2/FS/ProT_kruskal.csv')
ProT_kruskal = ProT_kruskal.drop(ProT_kruskal.columns[0], axis=1)
Roll_kruskal = pd.read_csv('../../data/trial2/FS/Roll_kruskal.csv')
Roll_kruskal = Roll_kruskal.drop(Roll_kruskal.columns[0], axis=1)
MGW_pca = pd.read_csv('../../data/trial2/FS/MGW_pca.csv')
MGW_pca = MGW_pca.drop(MGW_pca.columns[0], axis=1)
HelT_pca = pd.read_csv('../../data/trial2/FS/HelT_pca.csv')
HelT_pca = HelT_pca.drop(HelT_pca.columns[0], axis=1)
ProT_pca = pd.read_csv('../../data/trial2/FS/ProT_pca.csv')
ProT_pca = ProT_pca.drop(ProT_pca.columns[0], axis=1)
Roll_pca = pd.read_csv('../../data/trial2/FS/Roll_pca.csv')
Roll_pca = Roll_pca.drop(Roll_pca.columns[0], axis=1)

dnashape_fcbf = pd.read_csv('../../data/trial2/FS/dnashape_fcbf.csv')
dnashape_fcbf = dnashape_fcbf.drop(dnashape_fcbf.columns[0], axis=1)

Mono_freq = pd.read_csv('../../data/trial2/FS/nt1.csv')
Mono_freq = Mono_freq.drop(Mono_freq.columns[0], axis=1)
Di_freq = pd.read_csv('../../data/trial2/FS/nt2.csv')
Di_freq = Di_freq.drop(Di_freq.columns[0], axis=1)
Tetra_freq_kruskal = pd.read_csv('../../data/trial2/FS/nt4_kruskal.csv')
Tetra_freq_kruskal = Tetra_freq_kruskal.drop(Tetra_freq_kruskal.columns[0], axis=1)
Tetra_freq_pca = pd.read_csv('../../data/trial2/FS/nt4_pca.csv')
Tetra_freq_pca = Tetra_freq_pca.drop(Tetra_freq_pca.columns[0], axis=1)
Tetra_freq_fcbf = pd.read_csv('../../data/trial2/FS/nt4_fcbf.csv')
Tetra_freq_fcbf = Tetra_freq_fcbf.drop(Tetra_freq_fcbf.columns[0], axis=1)

list_features = [MGW_kruskal,HelT_kruskal,ProT_kruskal,Roll_kruskal,MGW_pca,HelT_pca,ProT_pca,Roll_pca,dnashape_fcbf,Mono_freq,Di_freq,Tetra_freq_kruskal,Tetra_freq_pca,Tetra_freq_fcbf]   
names_features = ["MGW_kruskal","HelT_kruskal","ProT_kruskal","Roll_kruskal","MGW_pca","HelT_pca","ProT_pca","Roll_pca","dnashape_fcbf","Mono_freq","Di_freq","Tetra_freq_kruskal","Tetra_freq_pca","Tetra_freq_fcbf"]

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
plt.savefig('../../data/trial2/results/randomforest.png')
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
plt.savefig('../../data/trial2/results/XGboost.png')
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
plt.savefig('../../data/trial2/results/KNN.png')
plt.show()
plt.close()



# ML plot summary for all the features and all the models
Models_test = pd.DataFrame()
Models_test["Accuracy"] = accuracy_RF+accuracy_KNN
Models_test["Number of features"] = names_features*2
Models_test["ML Algorithms"] = ["RandomForest"]*14+["KNN"]*14
Models_test


import seaborn as sns
sns.set_theme(style="darkgrid")

sns.set(rc={'figure.figsize':(15,15)})
sns_plot = sns.relplot(x="Number of features",y="Accuracy",hue="ML Algorithms",data=Models_test)
sns_plot.set_xticklabels(rotation=90)
sns_plot.savefig("../../data/trial2/results/Models_test.png")
plt.close()
