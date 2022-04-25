import pandas as pd
import seaborn as sns
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestRegressor
from sklearn.feature_selection import VarianceThreshold
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

#load data set 
df = pd.read_csv('E3_ubiquitin_ligase_04_bioactivity_data_3class_pIC50_pubchem_fp.csv')

#target variable for learning drug potency, aka input feature 
X = df.drop('pIC50', axis=1)

#output feature
Y = df.pIC50

#Remove low variance features
selection = VarianceThreshold(threshold=(.8 * (1 - .8)))    
X = selection.fit_transform(X)

#data split 80/20
X_train, X_test, Y_train, Y_test = train_test_split(X, Y, test_size=0.2)

#build regression model with random forest
np.random.seed(100)
model = RandomForestRegressor(n_estimators=100)
model.fit(X_train, Y_train)
r2 = model.score(X_test, Y_test)
print(r2)
Y_pred = model.predict(X_test)
#Scatter Plot of Experimental vs Predicted pIC50 Values
sns.set(color_codes=True)
sns.set_style("white")

plt.figure(figsize=(5, 10))
sns.set_theme(style="whitegrid")
ax = sns.regplot(Y_test, Y_pred, scatter_kws={'alpha':0.4})
ax.set_xlabel('Experimental pIC50', fontsize='large', fontweight='bold')
ax.set_ylabel('Predicted pIC50', fontsize='large', fontweight='bold')
ax.set_xlim(0, 12)
ax.set_ylim(0, 12)
ax.figure.set_size_inches(5, 5)
plt.savefig('scatter_plot_experimental_vs_predicted.pdf',dpi=300, bbox_inches = "tight")
