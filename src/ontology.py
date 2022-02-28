from owlready2 import *
from src.generate_data import all_curated_genes, therapies, diseases
import re
import types 
import defopt
import pandas as pd
import os

def add_oncotree(onto):
	oncotree = diseases()
	node = oncotree['TISSUE']
	stack = [node]
	while len(stack) > 0:
		node = stack.pop()
		print(node['code'],node['parent'])
		for key, child in node['children'].items():
			stack.append(child)

		parent = node['parent']
		if parent == "TISSUE":
			parent = "Disease"
		if node['code'] != "TISSUE":
			NewClass = types.new_class(node['code'], (onto[parent],))
	return onto

def main(variants_path: str, output_path: str):
	onto = get_ontology("ontology/base.owl").load()

	with onto:
		# generate gene subclasses
		for i in all_curated_genes():
			gene_subclass = types.new_class(i['hugoSymbol'], (onto['Gene'],))
			gene_subclass.comment = i['background'].replace('','')

		# generate therapy regimen subclassses
		for i in therapies():
			therapy_regimen = i
			therapy_regimen = i.replace(" ","_").replace(",","").replace("+","").replace("__","_")
			therapy_subclass = types.new_class(therapy_regimen, (onto['TherapyRegimen'],))

		onto = add_oncotree(onto)

		
		# generate  variant subclasses and is_biomarker_for object properties
		level_1 = pd.read_csv(variants_path).dropna()
		level_1 = level_1.loc[level_1.Hugo_Symbol!='Hugo_Symbol']
		for index, row in level_1.iterrows():
			variant_subclass = types.new_class(row['HGVSp_Short'], (onto[row['Hugo_Symbol']],))
			variant_subclass.is_a.append(onto['Variant'])
			for regimen in row['LEVEL_1'].split(",")[0]:
				therapy_regimen = regimen.replace(" ","_").replace(",","").replace("+","").replace("__","_")
			if onto[therapy_regimen] is not None:
				variant_relationship = onto[row['HGVSp_Short']].hasEvidenceLevel1 = [onto[therapy_regimen]]




	# save ontology to disk
	onto.save(file = output_path)

if __name__ == "__main__":
	defopt.run(main)

