import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski
import seaborn as sns
sns.set(style='ticks')
import matplotlib.pyplot as plt

df = pd.read_csv('bioactivity_preprocessed_data.csv')

#Inspired by: https://codeocean.com/explore/capsules?query=tag:data-curation
#this function finds the lipinski values for each drug, bases on ADME A.k.A. pharmocokinetic profile
#lipinski rule:
# -Molecular weight < 500 Dalton
# - Octanol-water partition coefficient (LogP) < 5
# - Hydrogen bond donors < 5
# - Hydrogen bond acceptors < 10

def lipinski(smiles, verbose=False):

    moldata= []
    for elem in smiles:
        mol=Chem.MolFromSmiles(elem) 
        moldata.append(mol)
       
    baseData= np.arange(1,1)
    i=0  
    for mol in moldata:        
       
        desc_MolWt = Descriptors.MolWt(mol)
        desc_MolLogP = Descriptors.MolLogP(mol)
        desc_NumHDonors = Lipinski.NumHDonors(mol)
        desc_NumHAcceptors = Lipinski.NumHAcceptors(mol)
           
        row = np.array([desc_MolWt,
                        desc_MolLogP,
                        desc_NumHDonors,
                        desc_NumHAcceptors])   
    
        if(i==0):
            baseData=row
        else:
            baseData=np.vstack([baseData, row])
        i=i+1      
    
    columnNames=["MW","LogP","NumHDonors","NumHAcceptors"]   
    descriptors = pd.DataFrame(data=baseData,columns=columnNames)
    
    return descriptors

#conver IC50 to PIC50 for data to be more uniformly distributed, 
#convert IC50 using negative logarithmic scale which is -log10(IC50).
#This function pIC50() will accept a DataFrame as input and will: 
#Take the IC50 values from the standard_value column and converts it from nM to M by multiplying the value by 10^-9
#Take the molar value and apply -log10
#Delete the standard_value column and create a new pIC50 column
def pIC50(input):
    pIC50 = []

    for i in input['standard_value_norm']:
        molar = i*(10**-9) # Converts nM to M
        pIC50.append(-np.log10(molar))

    input['pIC50'] = pIC50
    x = input.drop('standard_value_norm', 1)
        
    return x

#This fucntion takes Values greater than 100,000,000 and fixes them at 100,000,000 otherwise the -log() would make the values negative
# this caps the -log() value at 1 for PIC50 values >= to 100,000,000
def norm_value(input):
    norm = []

    for i in input['standard_value']:
        if i > 100000000:
          i = 100000000
        norm.append(i)

    input['standard_value_norm'] = norm
    x = input.drop('standard_value', 1)
        
    return x
#function modified from machinelearningmastery.com
def mannwhitney(descriptor, verbose=False):
  # https://machinelearningmastery.com/nonparametric-statistical-significance-tests-in-python/
  from numpy.random import seed
  from numpy.random import randn
  from scipy.stats import mannwhitneyu

# seed the random number generator
  seed(1)

# actives and inactives
  selection = [descriptor, 'bioactivity_class']
  df = df_2class[selection]
  active = df[df.bioactivity_class == 'active']
  active = active[descriptor]

  selection = [descriptor, 'bioactivity_class']
  df = df_2class[selection]
  inactive = df[df.bioactivity_class == 'inactive']
  inactive = inactive[descriptor]

# compare samples
  stat, p = mannwhitneyu(active, inactive)
  #print('Statistics=%.3f, p=%.3f' % (stat, p))

# interpret
  alpha = 0.05
  if p > alpha:
    interpretation = 'Same distribution (fail to reject H0)'
  else:
    interpretation = 'Different distribution (reject H0)'
  
  results = pd.DataFrame({'Descriptor':descriptor,
                          'Statistics':stat,
                          'p':p,
                          'alpha':alpha,
                          'Interpretation':interpretation}, index=[0])
  filename = 'mannwhitneyu_' + descriptor + '.csv'
  results.to_csv(filename)

  return results


df_lipinski = lipinski(df.canonical_smiles)

#combine the 2 DataFrames
df_combined = pd.concat([df,df_lipinski], axis=1)

#calculate the norm values
df_norm = norm_value(df_combined)

#find the pIC50 values for uniform data distrubution 
df_final = pIC50(df_norm)

#save final dataframe
df_final.to_csv('E3_ubiquitin_ligase_04_bioactivity_data_3class_pIC50.csv')

#remove intermediate drugs
df_2class = df_final[df_final.bioactivity_class != 'intermediate']


#Frequency plot of the 2 bioactivity classes
plt.figure(figsize=(5.5, 5.5))

sns.countplot(x='bioactivity_class', data=df_2class, edgecolor='black')

plt.xlabel('Bioactivity class', fontsize=14, fontweight='bold')
plt.ylabel('Frequency', fontsize=14, fontweight='bold')
#save plot as pdf
plt.savefig('plot_bioactivity_class.pdf')

#Scatter plot of MW versus LogP
plt.figure(figsize=(5.5, 5.5))

sns.scatterplot(x='MW', y='LogP', data=df_2class, hue='bioactivity_class', size='pIC50', edgecolor='black', alpha=0.7)

plt.xlabel('MW', fontsize=14, fontweight='bold')
plt.ylabel('LogP', fontsize=14, fontweight='bold')
plt.legend(bbox_to_anchor=(1.05,1), loc=2, borderaxespad = 0)
#save plot as pdf
plt.savefig('plot_MW_vs_LogP.pdf',dpi=300, bbox_inches = "tight")

#box plots of PCI50 value
plt.figure(figsize=(5.5, 5.5))

sns.boxplot(x = 'bioactivity_class', y = 'pIC50', data = df_2class)

plt.xlabel('Bioactivity class', fontsize=14, fontweight='bold')
plt.ylabel('pIC50 value', fontsize=14, fontweight='bold')
#save plot as pdf
plt.savefig('plot_ic50.pdf',dpi=300, bbox_inches = "tight")

#apply mannwhitney function to pIC50 value to look for statistical siginifigance for pIC50 variable
mannwhitney('pIC50')
plt.figure(figsize=(5.5, 5.5))

sns.boxplot(x = 'bioactivity_class', y = 'MW', data = df_2class)

plt.xlabel('Bioactivity class', fontsize=14, fontweight='bold')
plt.ylabel('MW', fontsize=14, fontweight='bold')
#save plot 
plt.savefig('plot_MW.pdf',dpi=300, bbox_inches = "tight")

#apply mannwhitney to the other lipinski descriptors
mannwhitney('MW')

plt.figure(figsize=(5.5, 5.5))

sns.boxplot(x = 'bioactivity_class', y = 'LogP', data = df_2class)

plt.xlabel('Bioactivity class', fontsize=14, fontweight='bold')
plt.ylabel('LogP', fontsize=14, fontweight='bold')
#save plot 
plt.savefig('plot_LogP.pdf',dpi=300, bbox_inches = "tight")

mannwhitney('LogP')

plt.figure(figsize=(5.5, 5.5))

sns.boxplot(x = 'bioactivity_class', y = 'NumHDonors', data = df_2class)

plt.xlabel('Bioactivity class', fontsize=14, fontweight='bold')
plt.ylabel('NumHDonors', fontsize=14, fontweight='bold')
#save plot
plt.savefig('plot_NumHDonors.pdf',dpi=300, bbox_inches = "tight")

mannwhitney('NumHDonors')

plt.figure(figsize=(5.5, 5.5))

sns.boxplot(x = 'bioactivity_class', y = 'NumHAcceptors', data = df_2class)

plt.xlabel('Bioactivity class', fontsize=14, fontweight='bold')
plt.ylabel('NumHAcceptors', fontsize=14, fontweight='bold')
#save plot
plt.savefig('plot_NumHAcceptors.pdf',dpi=300, bbox_inches = "tight")

mannwhitney('NumHAcceptors')

