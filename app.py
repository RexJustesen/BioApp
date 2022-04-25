import streamlit as st
import pandas as pd
from PIL import Image
import subprocess
import os
import base64
import pickle
import glob
from padelpy import padeldescriptor

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

# Molecular descriptor calculator
def desc_calc():
    # Performs the descriptor calculation
    bashCommand = "java -Xms2G -Xmx2G -Djava.awt.headless=true -jar ./PaDEL-Descriptor/PaDEL-Descriptor.jar -removesalt -standardizenitro -fingerprints -descriptortypes ./PaDEL-Descriptor/descriptors.xml/PubchemFingerprinter -dir ./ -file descriptors_output.csv"
    process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
    output, error = process.communicate()
    os.remove('molecule.smi')

# File download
def filedownload(df):
    csv = df.to_csv(index=False)
    b64 = base64.b64encode(csv.encode()).decode()  # strings <-> bytes conversions
    href = f'<a href="data:file/csv;base64,{b64}" download="prediction.csv">Download Predictions</a>'
    return href

# Model building
def build_model(input_data):
    # Reads in saved regression model
    load_model = pickle.load(open('bioactivity_prediction_model.pkl', 'rb'))
    # Apply model to make predictions
    prediction = load_model.predict(input_data)
    st.header('**Prediction output**')
    prediction_output = pd.Series(prediction, name='pIC50')
    molecule_name = pd.Series(load_data[1], name='molecule_name')
    df = pd.concat([molecule_name, prediction_output], axis=1)
    st.write(df)
    st.markdown(filedownload(df), unsafe_allow_html=True)

# Logo image
#image = Image.open('logo.png')

#st.image(image, use_column_width=True)

# Page title
st.markdown("""
# Bioactivity Prediction App (E3 Ubiquitin Ligase)
This app allows you to predict the bioactivity towards inhibting the `E3 Ubiquitin Ligase` enzyme. `E3 Ubiquitin Ligase` is a drug target for Parkinson's disease.


""")

# Sidebar
with st.sidebar.header('1. Upload your CSV data'):
    uploaded_file = st.sidebar.file_uploader("Upload your input file", type=['txt'])
    st.sidebar.markdown("""
[Example input file](example_E3_Ubiquitin_Ligase)
""")

if st.sidebar.button('Predict'):
    load_data = pd.read_table(uploaded_file, sep=' ', header=None)
    load_data.to_csv('molecule.smi', sep='\t', index=False, header=False)
    


    st.header('**Original input data**')
    st.write(load_data)

   
    with st.spinner("Calculating descriptors..."):
        desc_calc()

    # Read in calculated descriptors and display the dataframe
    st.header('**Calculated molecular descriptors**')
    desc = pd.read_csv('descriptor_output.csv')
    st.write(desc)
    st.write(desc.shape)


    # Read descriptor list used in previously built model
    st.header('**Subset of descriptors from previously built models**')
    Xlist = list(pd.read_csv('PubChemdescriptor_data_.csv').columns)
    desc_subset = desc[Xlist]
    st.write(desc_subset)
    st.write(desc_subset.shape)

    # Apply trained model to make prediction on query compounds
    build_model(desc_subset)
else:
    st.info('Upload input data in the sidebar to start!')