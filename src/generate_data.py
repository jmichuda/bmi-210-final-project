import requests
from src.constants import ONCOKB_API_KEY
import pandas as pd
import types 

def ontology_classes(onto):
	classes = []
	for i in list(onto.classes()):
		classes.append((str(i)[5:]))
	for i in list(onto.individuals()):
		classes.append((str(i)[5:]))
	return classes

def therapy_normalize(diagnosis):
	return diagnosis.replace(" ","_").replace(",","").replace("+","_").replace("__","_").replace("__","_")



## add oncokb curated genes
def all_curated_genes(onto):
	HEADER = {
	'Authorization' : f'Bearer {ONCOKB_API_KEY}',
	'accept': 'application/json'
	}
	response = requests.get('https://www.oncokb.org/api/v1/utils/allCuratedGenes',headers = HEADER)

	# generate gene subclasses
	for i in response.json():
		gene_subclass = onto['Gene'](i['hugoSymbol'])
		gene_subclass.comment = i['background'].replace('','')
		gene_subclass.grch38RefSeq = i['grch38RefSeq']
		gene_subclass.highestResistanceLevel = i['highestResistanceLevel']
		gene_subclass.hugoSymbol = i['hugoSymbol']
		gene_subclass.entrezGeneId = i['entrezGeneId']
		gene_subclass.grch38Isoform = i['grch38Isoform']
		gene_subclass.highestSensitiveLevel = i['highestSensitiveLevel']
		gene_subclass.oncogene = i['oncogene']
		gene_subclass.hasVariant = []
	return onto


## add therapies
def therapies(onto):
	tabular = pd.read_csv('src/oncokb_biomarker_drug_associations.tsv',sep="\t")
	therapies = list(tabular['Drugs (for therapeutic implications only)'].dropna().unique())
	for i in therapies:
		therapy_regimen = therapy_normalize(i)
		therapy_subclass = onto['TherapyRegimen'](therapy_regimen)
	return onto

## Add oncotree cancer types

def oncotree():
	HEADER = {
	'accept': 'application/json'
	}
	response = requests.get('http://oncotree.mskcc.org/api/tumorTypes/tree',headers = HEADER)
	return (response.json())

def add_oncotree(onto):
	tree = oncotree()
	node = tree['TISSUE']
	stack = [node]
	while len(stack) > 0:
		node = stack.pop()
		for key, child in node['children'].items():
			stack.append(child)
		parent = node['parent']
		if parent == "TISSUE":
			parent = "Disease"
		if node['code'] != "TISSUE":
			if node['code'] not in ontology_classes(onto):
				NewClass = types.new_class(node['code'], (onto[parent],))
	return onto

def clean_mutation_effect(mutation_effect):
	mapping = {
		"Unknown": None,
		"Likely Loss-of-function":"LossOfFunction",
		"Gain-of-function":"GainOfFunction",
		"Loss-of-function":"LossOfFunction",
		"Likely Gain-of-function":"GainOfFunction",
		"Likely Neutral":None,
		"Inconclusive":None,
		"Likely Switch-of-function":None,
		"Switch-of-function":None,
		None: None
	}
	return mapping[mutation_effect]



## vaiants and biomarker relationships
def parse_maf(onto, variants_path):
	# generate  variant subclasses and is_biomarker_for object properties
	variants = pd.read_csv(variants_path)


	variants = variants.loc[
		(variants.Hugo_Symbol!='Hugo_Symbol') 
		& variants.Hugo_Symbol.notna()
		* (variants.Hugo_Symbol!='nan') 
	]


	for index, row in variants.iterrows():
		mutation_effect = clean_mutation_effect(row['MUTATION_EFFECT'])
		variant = onto['Mutation'](row['HGVSp_Short'])
		if mutation_effect is not None and onto[row['Hugo_Symbol']] is not None:
			variant.is_a.append(onto[mutation_effect])
			variant.hasGene = onto[row['Hugo_Symbol']]
			onto[row['Hugo_Symbol']].hasVariant.append(onto[row['HGVSp_Short']])




		# gene = types.new_class(row['Hugo_Symbol'], (onto['Gene'],))
		# gene.hasVariant = [variant]
		# regimens = row['LEVEL_1']
		# for i in row['LEVEL_1'].split(","):
		# 	therapy_regimen = therapy_normalize(i)
		# 	therapy = types.new_class(therapy_regimen, (onto['TherapyRegimen'],))
		# 	therapy.hasEvidenceLevel1 =  [variant]

	print(onto['ATM'].hasVariant)
	return onto

