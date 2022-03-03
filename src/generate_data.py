import requests
from src.constants import ONCOKB_API_KEY, MAF_COLUMNS,THERAPEUTIC_COLUMNS,DEVELOPMENT_MODE
import pandas as pd
import types 
import numpy as np

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
		gene_subclass = types.new_class(i['hugoSymbol'],(onto['Gene'],))
		gene_subclass.comment = i['background'].replace('','')
		gene_subclass.grch38RefSeq = i['grch38RefSeq']
		gene_subclass.highestResistanceLevel = i['highestResistanceLevel']
		gene_subclass.hugoSymbol = i['hugoSymbol']
		# gene_subclass.entrezGeneId = i['entrezGeneId']
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
		therapy_subclass = types.new_class(therapy_regimen,(onto['TherapyRegimen'],))
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


def clean_variant_classification(variant_classification):
	mapping = {
	"Missense_Mutation":"Missense",
	"Nonsense_Mutation":"Nonsense",
	"Frame_Shift_Del":"Frameshift",
	"Splice_Site":"Splice",
	"Frame_Shift_Ins":"Frameshift",
	"In_Frame_Del":"INDEL",
	"In_Frame_Ins":"INDEL",
	"Translation_Start_Site":"TranslationStartSite",
	}
	return mapping[variant_classification]

def clean_variant(variant):
	return variant.replace("p.","").replace("*","")



## variants and biomarker relationships
def parse_maf(onto, variants_path):
	# generate  variant subclasses and is_biomarker_for object properties
	variants = pd.read_csv(variants_path, sep ="\t",low_memory=False)
	variants = variants[MAF_COLUMNS]
	variants = variants.loc[variants[THERAPEUTIC_COLUMNS].notna().any(axis=1)]

	if DEVELOPMENT_MODE:
		variants = variants.sample(100)

	for index, row in variants.iterrows():
		mutation_effect = clean_mutation_effect(row['MUTATION_EFFECT'])
		variant_classification = clean_variant_classification(row['Variant_Classification'])
		variant_name = clean_variant(row['HGVSp_Short'])
		cancer_type = row['Cohort']


		biomarker_name = f"{row['Hugo_Symbol']}_{variant_name}_{row['Cohort']}"
		biomarker = types.new_class(biomarker_name, (onto['Biomarker'],))
		
		if variant_classification is not None:
			variant = types.new_class(variant_name,(onto[variant_classification],))
			biomarker.hasVariant = [variant]
			variant.hasBiomarker.append(biomarker)

			gene  = onto[row['Hugo_Symbol']]
			variant.hasGene = [gene]
			gene.hasVariant.append(variant)


			level_mapping = {
				"Level_1":["LEVEL_1"],
				"Level_2":["LEVEL_2"],
				"Level_3":["LEVEL_3A","LEVEL_3B"],
				"Level_4":["LEVEL_4"],
				"Level_R1":["LEVEL_R1"],
				"Level_R2":["LEVEL_R2"],
			}
			for key, value in level_mapping.items():
				level = ','.join(row[value].dropna())
				if level != '':
					for i in level.split(","):
						therapy_regimen = therapy_normalize(i)
						if therapy_regimen in ontology_classes(onto):
							therapy = onto[therapy_regimen]
						else:
							therapy = types.new_class(therapy_regimen,(onto['TherapyRegimen'],))
						if key == "Level_1":
							therapy.hasEvidenceLevel1.append(biomarker)
							biomarker.hasEvidenceLevel1.append(therapy)
						if key == "Level_2":
							therapy.hasEvidenceLevel2.append(biomarker)
							biomarker.hasEvidenceLevel2.append(therapy)
						if key == "Level_3":
							therapy.hasEvidenceLevel3.append(biomarker)
							biomarker.hasEvidenceLevel3.append(therapy)
						if key == "Level_4":
							therapy.hasEvidenceLevel4.append(biomarker)
							biomarker.hasEvidenceLevel4.append(therapy)
						if key == "Level_R1":
							therapy.hasEvidenceLevelR1.append(biomarker)
							biomarker.hasEvidenceLevelR1.append(therapy)
						if key == "Level_R2":
							therapy.hasEvidenceLevelR2.append(biomarker)
							biomarker.hasEvidenceLevelR2.append(therapy)
					

	return onto

