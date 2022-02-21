from owlready2 import *
from query_oncokb import all_curated_genes, therapies
import re
import types 
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


onto.save(file = "oncokb.owl")

