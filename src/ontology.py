from owlready2 import *
from generate_data import all_curated_genes, therapies, variants
import re
import types 
import defopt

def main(variants_path: str, output_path: str):
	onto = get_ontology("http://test.org/onto.owl") 

	with onto:
		class Gene(Thing):
			pass

		class Oncogene(DataProperty):
			domain = [Gene]
			range = [bool]

		for i in all_curated_genes():
			gene_subclass = types.new_class(i['hugoSymbol'], (Gene,))
			gene_subclass.comment = i['background'].replace('','')
			gene_subclass.Oncogene = [i['oncogene']]

		class TherapyRegimen(Thing):
			pass

		for i in therapies():
			therapy_regimen = i
			therapy_regimen = i.replace(" ","_").replace(",","").replace("+","").replace("__","_")
			therapy_subclass = types.new_class(therapy_regimen, (TherapyRegimen,))


		class Variants():
			pass
		class LevelOneRelationship():
			pass

	onto.save(file = output_path)

