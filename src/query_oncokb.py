import requests
from constants import ONCOKB_API_KEY,HEADER


def all_curated_genes():
	response = requests.get('https://www.oncokb.org/api/v1/utils/allCuratedGenes',headers = HEADER)
	return (response.json())