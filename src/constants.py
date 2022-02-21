import os
ONCOKB_API_KEY = "28c66f60-7f3d-48d9-b791-4babf383f3ad" #os.get('ONCOKB_API_KEY')

HEADER = {
	'Authorization' : f'Bearer {ONCOKB_API_KEY}',
	'accept': 'application/json'
	}


