import requests
from constants import ONCOKB_API_KEY,HEADER
import pandas as pd

def all_curated_genes():
	response = requests.get('https://www.oncokb.org/api/v1/utils/allCuratedGenes',headers = HEADER)
	return (response.json())

def therapies():
	tabular = pd.read_csv('oncokb_biomarker_drug_associations.tsv',sep="\t")
	return list(tabular['Drugs (for therapeutic implications only)'].dropna().unique())

