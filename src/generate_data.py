import requests
from src.constants import ONCOKB_API_KEY
import pandas as pd

def all_curated_genes():
	HEADER = {
	'Authorization' : f'Bearer {ONCOKB_API_KEY}',
	'accept': 'application/json'
	}
	response = requests.get('https://www.oncokb.org/api/v1/utils/allCuratedGenes',headers = HEADER)
	return (response.json())

def therapies():
	tabular = pd.read_csv('src/oncokb_biomarker_drug_associations.tsv',sep="\t")
	return list(tabular['Drugs (for therapeutic implications only)'].dropna().unique())


def diseases():
	HEADER = {
	'accept': 'application/json'
	}
	response = requests.get('http://oncotree.mskcc.org/api/tumorTypes/tree',headers = HEADER)
	return (response.json())

