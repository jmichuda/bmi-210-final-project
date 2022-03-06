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
	"Nonstop_Mutation":"Nonsense"
	}
	return mapping[variant_classification]

def clean_variant(variant):
	return variant.replace("p.","").replace("*","")

def add_levels(onto, row, biomarker):
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
				therapy_name = therapy_normalize(i)
				if therapy_name in ontology_classes(onto):
					therapy_regimen = onto[therapy_name]
				else:
					therapy_regimen = types.new_class(therapy_name,(onto['TherapyRegimen'],))
				onto, therapy_regimen, biomarker = add_evidence(onto, therapy_regimen, biomarker, key)
	return onto,biomarker

def add_evidence(onto, therapy, biomarker, key):
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
	return onto, therapy, biomarker

## variants and biomarker relationships
def parse_maf(onto, variants_path):
	# generate  variant subclasses and is_biomarker_for object properties
	variants = pd.read_csv(variants_path, sep ="\t",low_memory=False)
	variants = variants[MAF_COLUMNS]
	variants = variants.loc[variants['HIGHEST_LEVEL'].notna()]

	if DEVELOPMENT_MODE:
		variants = variants.sample(100)

	for index, row in variants.iterrows():
		mutation_effect = clean_mutation_effect(row['MUTATION_EFFECT'])
		variant_classification = clean_variant_classification(row['Variant_Classification'])
		variant_name = clean_variant(row['HGVSp_Short'])
		cancer_type = onto[row['Cohort']]


		biomarker_name = f"{row['Hugo_Symbol']}_{variant_name}_{row['Cohort']}"
		biomarker = types.new_class(biomarker_name, (onto['Biomarker'],))
		gene  = onto[row['Hugo_Symbol']]
		biomarker.hasGene.append(gene)
		gene.hasBiomarker.append(biomarker)
		biomarker.hasDisease.append(cancer_type)
		biomarker.evidenceSource.append("oncokb")
		if variant_classification is not None:
			variant = types.new_class(variant_name,(onto[variant_classification],))
			biomarker.hasVariant = [variant]
			variant.hasBiomarker.append(biomarker)
			variant.hasGene = [gene]
			gene.hasVariant.append(variant)
			onto, biomarker = add_levels(onto, row,biomarker)

	return onto


def add_fusions(onto,fusion_path):
	fusions = pd.read_csv(fusion_path,sep="\t",low_memory=False)
	fusions = fusions.loc[fusions['HIGHEST_LEVEL'].notna()]
	for index, row in fusions.iterrows():
		gene_1 = row['Fusion'].split("-")[0]
		gene_2 = row['Fusion'].split("-")[1]
		biomarker_name = f"fusion_{gene_1}_{gene_2}"
		fusion_name = f"{gene_1}_{gene_2}"
		if fusion_name not in ontology_classes(onto):
			fusion = types.new_class(fusion_name, (onto['GeneFusion'],))
			biomarker = types.new_class(biomarker_name, (onto['Biomarker'],))
		else:
			fusion = onto[fusion_name]
			biomarker = onto[biomarker_name]

		if gene_1 in ontology_classes(onto):
			fusion.hasGene = [onto[gene_1]]
		if gene_2 in ontology_classes(onto):
			fusion.hasGene.append(onto[gene_2])
		biomarker.hasVariant = [fusion]
		fusion.hasBiomarker.append(biomarker)
		onto, biomarker = add_levels(onto, row, biomarker)
		biomarker.evidenceSource.append("oncokb")
	return onto

def add_cnas(onto, cna_path):
	cnas = pd.read_csv(cna_path, sep = "\t",low_memory = False).drop_duplicates(subset=['CANCER_TYPE','HUGO_SYMBOL','ALTERATION'])
	cnas = cnas.loc[cnas['HIGHEST_LEVEL'].notna()]
	for index, row in cnas.iterrows():
		gene_name = row['HUGO_SYMBOL']
		alteration = row['ALTERATION']
		cna_name = f"{gene_name}_{alteration}"
		cancer_type = row['CANCER_TYPE']
		biomarker_name = f"{gene_name}_{alteration}_{cancer_type}"
		disease = onto[cancer_type]
		if gene_name not in ontology_classes(onto):
			gene = types.new_class(gene_name, (onto['Gene']))
		else:
			gene = onto[gene_name]

		cna = types.new_class(cna_name, (onto[alteration],))
		biomarker = types.new_class(biomarker_name,(onto['Biomarker'],))

		cna.hasGene.append(gene)
		cna.hasBiomarker.append(biomarker)
		biomarker.hasVariant.append(cna)
		biomarker.hasDisease.append(disease)
		gene.hasVariant.append(cna)
		biomarker.evidenceSource.append("oncokb")
		onto, biomarker = add_levels(onto, row, biomarker)
	return onto


def map_civic_evidence(clin_sig, evidence_level):
	mapping = {
		"A":{"Sensitivity/Response":"Level_2", "Resistance":"Level_R1", "Adverse Response":"Level_R1", "Reduced Sensitivity":"Level_R1",},
		"B":{"Sensitivity/Response":"Level_3", "Resistance":"Level_R1", "Adverse Response":"Level_R1", "Reduced Sensitivity":"Level_R1",},
		"C":{"Sensitivity/Response":"Level_3", "Resistance":"Level_R2", "Adverse Response":"Level_R2", "Reduced Sensitivity":"Level_R2",},
		"D":{"Sensitivity/Response":"Level_3", "Resistance":"Level_R2", "Adverse Response":"Level_R2", "Reduced Sensitivity":"Level_R2",},
		"E":{"Sensitivity/Response":"Level_4", "Resistance":"Level_R2", "Adverse Response":"Level_R2", "Reduced Sensitivity":"Level_R2",},
	}
	return mapping[evidence_level][clin_sig]

def add_civic_variants(onto, civic_path):
	civic_evidence = pd.read_csv(civic_path)
	# subset to cleanly formatted civic variants
	civic_evidence= civic_evidence.loc[civic_evidence.variant.str.contains("^[A-Z][0-9]*[A-Z]$",na=False)]
	civic_evidence = civic_evidence[['Gene','variant','TherapyRegimen','oncotree','ClinicalSignificance','EvidenceLevel']].dropna().drop_duplicates()
	for index, row in civic_evidence.iterrows():
		therapy_name = therapy_normalize(row['TherapyRegimen'])
		gene_name = row['Gene']
		mutation_name = clean_variant(row['variant'])
		disease_name = row['oncotree']
		evidence_level = map_civic_evidence(row['ClinicalSignificance'], row['EvidenceLevel'])
		biomarker_name = f"{gene_name}_{mutation_name}_{disease_name}"
		if disease_name in ontology_classes(onto):
			disease = onto[disease_name]
		else:
			disease = types.new_class(disease_name, (onto['Disease'],))

		if gene_name in ontology_classes(onto):
			gene = onto[gene_name]
			
		else:
			gene = types.new_class(gene_name, (onto['Gene'],))

		if therapy_name in ontology_classes(onto):
			therapy_regimen = onto[therapy_name]
		else:
			therapy_regimen = types.new_class(therapy_name, (onto['TherapyRegimen'],))

		if mutation_name in ontology_classes(onto):
			mutation = onto[mutation_name]
		else:
			mutation = types.new_class(mutation_name, (onto['Missense'],))

		if biomarker_name in ontology_classes(onto):
			biomarker = onto[biomarker_name]
		else:
			biomarker = types.new_class(biomarker_name, (onto['Biomarker'],))

		biomarker.evidenceSource.append("civic")
		biomarker.hasGene.append(gene)
		gene.hasBiomarker.append(biomarker)
		biomarker.hasDisease.append(disease)
		biomarker.hasVariant.append(mutation)
		mutation.hasBiomarker.append(biomarker)

		onto, therapy_regimen, biomarker = add_evidence(onto, therapy_regimen, biomarker, evidence_level)



	return onto

def add_civic_cnas(onto, civic_path):
	civic_evidence = pd.read_csv(civic_path)
	# subset to cleanly formatted civic variants
	civic_cnas= civic_evidence.loc[civic_evidence.variant.str.contains("amplification|deletion",na=False, case=False)]
	civic_cnas = civic_cnas[['Gene','variant','TherapyRegimen','oncotree','ClinicalSignificance','EvidenceLevel']].dropna().drop_duplicates()
	for index, row in civic_cnas.iterrows():
		therapy_name = therapy_normalize(row['TherapyRegimen'])
		gene_name = row['Gene']
		alteration = "Amplification" if "amplification" in row['variant'] else "Deletion"
		cna_name = f"{gene_name}_{alteration}"
		disease_name = row['oncotree']
		evidence_level = map_civic_evidence(row['ClinicalSignificance'], row['EvidenceLevel'])
		biomarker_name = f"{gene_name}_{alteration}_{disease_name}"
		
		if disease_name in ontology_classes(onto):
			disease = onto[disease_name]
		else:
			disease = types.new_class(disease_name, (onto['Disease'],))

		if gene_name in ontology_classes(onto):
			gene = onto[gene_name]
			
		else:
			gene = types.new_class(gene_name, (onto['Gene'],))

		if therapy_name in ontology_classes(onto):
			therapy_regimen = onto[therapy_name]
		else:
			therapy_regimen = types.new_class(therapy_name, (onto['TherapyRegimen'],))

		if cna_name in ontology_classes(onto):
			cna = onto[mutation_name]
		else:
			cna = types.new_class(cna_name, (onto[alteration],))

		if biomarker_name in ontology_classes(onto):
			biomarker = onto[biomarker_name]
		else:
			biomarker = types.new_class(biomarker_name, (onto['Biomarker'],))
			


		biomarker.hasVariant.append(cna)
		biomarker.evidenceSource.append("oncokb")
		gene.hasBiomarker.append(biomarker)
		biomarker.hasDisease.append(disease)
		biomarker.hasVariant.append(cna)
		cna.hasBiomarker.append(biomarker)


		onto, therapy_regimen, biomarker = add_evidence(onto, therapy_regimen, biomarker, evidence_level)

	return onto


