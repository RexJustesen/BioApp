# Import necessary libraries
import pandas as pd
from chembl_webresource_client.new_client import new_client

# Target search for coronavirus
target = new_client.target
target_query = target.search('E3 ubiquitin ligase')
targets = pd.DataFrame.from_dict(target_query)
#assign the fifth entry (which corresponds to the target protein, 
#coronavirus 3C-like proteinase) to the selected_target variable
selected_target = targets.target_chembl_id[4]

#retrieve only bioactivity data for coronavirus 3C-like proteinase
#CHEMBL3927 that are reported as IC50 nanomolar values
activity = new_client.activity
res = activity.filter(target_chembl_id=selected_target).filter(standard_type="IC50")
df = pd.DataFrame.from_dict(res)

#save the resulting bioactivity data to a CSV file bioactivity_data.csv.
df.standard_type.unique()
df.to_csv('bioactivity_data.csv', index=False)

#f any compounds has missing value for the standard_value column then drop it
df2 = df[df.standard_value.notna()]

#Compounds having values of less than 1000 nM 
# will be considered to be active while those greater 
# than 10,000 nM will be considered to be inactive. 
# As for those values in between 1,000 and 10,000 nM 
# will be referred to as intermediate.
bioactivity_class = []
for i in df2.standard_value:
  if float(i) >= 10000:
    bioactivity_class.append("inactive")
  elif float(i) <= 1000:
    bioactivity_class.append("active")
  else:
    bioactivity_class.append("intermediate")

#Iterate the molecule_chembl_id to a list
mol_cid = []
for i in df2.molecule_chembl_id:
  mol_cid.append(i)

#Iterate canonical_smiles to a list
canonical_smiles = []
for i in df2.canonical_smiles:
  canonical_smiles.append(i)

#Iterate standard_value to a list
standard_value = []
for i in df2.standard_value:
  standard_value.append(i)

#Combine the 4 lists into a dataframe
data_tuples = list(zip(mol_cid, canonical_smiles, bioactivity_class, standard_value))
df3 = pd.DataFrame( data_tuples,  columns=['molecule_chembl_id', 'canonical_smiles', 'bioactivity_class', 'standard_value'])

#Saves dataframe to CSV file
df3.to_csv('bioactivity_preprocessed_data.csv', index=False)