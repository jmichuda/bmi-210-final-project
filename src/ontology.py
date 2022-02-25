from owlready2 import *
from src.generate_data import all_curated_genes, therapies
import re
import types 
import defopt
import pandas as pd
import os

def main(variants_path: str, output_path: str):
	print(os.path.exists("ontology/base.owl"))
	onto = get_ontology("ontology/base.owl").load()
	print(list(onto.classes()))

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

		# generate  variant subclasses and is_biomarker_for object properties
		level_1 = pd.read_csv(variants_path).dropna()
		level_1 = level_1.loc[level_1.Hugo_Symbol!='Hugo_Symbol']
		for index, row in level_1.iterrows():
			variant_subclass = types.new_class(row['HGVSp_Short'], (onto[row['Hugo_Symbol']],))
			variant_subclass.is_a.append(onto['Variant'])
			for regimen in row['LEVEL_1'].split(",")[0]:
				therapy_regimen = regimen.replace(" ","_").replace(",","").replace("+","").replace("__","_")
			if onto[therapy_regimen] is not None:
				variant_relationship = onto[row['HGVSp_Short']].is_biomarker_for = [onto[therapy_regimen]]

	# save ontology to disk
	onto.save(file = output_path)

if __name__ == "__main__":
	defopt.run(main)

