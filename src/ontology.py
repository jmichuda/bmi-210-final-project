from owlready2 import *
from src.generate_data import therapies, all_curated_genes,add_oncotree, parse_maf, add_fusions, add_cnas
import re
import types 
import defopt
import pandas as pd
import os





def main(maf_path: str, fusion_path: str, cna_path: str, output_path: str):
	onto = get_ontology("ontology/base.owl").load()
	with onto:
		onto = therapies(onto)
		onto = all_curated_genes(onto)
		onto = add_oncotree(onto)
		onto = add_oncotree(onto)
		onto = parse_maf(onto, maf_path)
		onto = add_fusions(onto,fusion_path)
		onto = add_cnas(onto, cna_path)
		onto.save(file = output_path)

if __name__ == "__main__":
	defopt.run(main)

