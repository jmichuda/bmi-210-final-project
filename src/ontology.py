from owlready2 import *
from query_oncokb import all_curated_genes
import re
onto = get_ontology("http://test.org/onto.owl") 

with onto:
	class Gene(Thing):
		pass

	class grch37Isoform(DataProperty):
		domain = [Gene]
		range = [str]



for i in all_curated_genes():
	print(i.keys())
	gene = Gene(i['hugoSymbol'])
	gene.comment = i['background'].replace('','')

	gene.grch37Isoform = [i['grch37Isoform']]
	break


onto.save(file = "oncokb.owl")

