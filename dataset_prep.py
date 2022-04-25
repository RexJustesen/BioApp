import pandas as pd
import glob
from padelpy import padeldescriptor

##listing and sortingv .xml files
xml_files = glob.glob("*.xml")
xml_files.sort()
xml_files
##Creating a list of present files
FP_list = ['AtomPairs2DCount',
 'AtomPairs2D',
 'EState',
 'CDKextended',
 'CDK',
 'CDKgraphonly',
 'KlekotaRothCount',
 'KlekotaRoth',
 'MACCS',
 'PubChem',
 'SubstructureCount',
 'Substructure']

 #Creating Data Dictionary
fp = dict(zip(FP_list, xml_files))
print(fp['PubChem'])
 #load in data file 
df3 = pd.read_csv('E3_ubiquitin_ligase_04_bioactivity_data_3class_pIC50.csv')

#prepare data set for padel
selection = ['canonical_smiles','molecule_chembl_id']
df3_selection = df3[selection]
df3_selection.to_csv('molecule.smi', sep='\t', index=False, header=False)

#Setting the fingerprint module
fingerprint = 'PubChem'
fingerprint_output_file = 'descriptor_data_'.join([fingerprint,'.csv']) #'PubChem.csv'
fingerprint_descriptortypes = fp[fingerprint]
padeldescriptor(mol_dir='molecule.smi', 
                d_file=fingerprint_output_file, #'PubChem.csv'
                #descriptortypes='PubChemFingerprint.xml', 
                descriptortypes= fingerprint_descriptortypes,
                standardizenitro=True,
                threads=2,
                removesalt=True,
                log=True,
                fingerprints=True)

#Preparing the X and Y Data Matrices
#X data matrix
df3_X = pd.read_csv(fingerprint_output_file)

#drop column names
df3_X = df3_X.drop(columns=['Name'])

#Convert IC50 to pIC50
df3_Y = df3['pIC50']

#Combining X and Y variable
dataset3 = pd.concat([df3_X,df3_Y], axis=1)

dataset3.to_csv('E3_ubiquitin_ligase_04_bioactivity_data_3class_pIC50_pubchem_fp.csv', index=False)
