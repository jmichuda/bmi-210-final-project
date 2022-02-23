from owlready2 import *
from src.generate_data import all_curated_genes, therapies
import re
import types 
import defopt
import pandas as pd

def main(variants_path: str, output_path: str):
	onto = get_ontology("http://test.org/onto.owl") 

	with onto:
		class Gene(Thing):
			pass


		for i in all_curated_genes():
			gene_subclass = types.new_class(i['hugoSymbol'], (Gene,))
			gene_subclass.comment = i['background'].replace('','')

		class TherapyRegimen(Thing):
			pass

		for i in therapies():
			therapy_regimen = i
			therapy_regimen = i.replace(" ","_").replace(",","").replace("+","").replace("__","_")
			therapy_subclass = types.new_class(therapy_regimen, (TherapyRegimen,))


		class is_biomarker_for(Gene >> TherapyRegimen):
			pass

		class has_biomarker(TherapyRegimen >> Gene):
			pass

		level_1 = pd.read_csv(variants_path).dropna()
		level_1 = level_1.loc[level_1.Hugo_Symbol!='Hugo_Symbol']

		for index, row in level_1.iterrows():
			variant_subclass = types.new_class(row['HGVSp_Short'], (onto[row['Hugo_Symbol']],))
			therapy_regimen = row['LEVEL_1'].split(",")[0].replace(" ","_").replace(",","").replace("+","").replace("__","_")
			if onto[therapy_regimen] is not None:
				variant_relationship = onto[row['HGVSp_Short']].is_biomarker_for = [onto[therapy_regimen]]

	onto.save(file = output_path)

if __name__ == "__main__":
	defopt.run(main)

