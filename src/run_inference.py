import pandas as pd
from owlready2 import *
import defopt
from src.generate_data import clean_variant
import multiprocessing
import functools
import itertools
import tqdm
import tempfile


onto = get_ontology("ontology/oncokb_civic.owl").load()
fd, path = tempfile.mkstemp()
print(path)
default_world.set_backend(filename = path)

def get_therapy_regimens(onto, gene, variant, disease, evidence_level):
    with onto:
        treatments = list(
            default_world.sparql(
                f"""
                    SELECT distinct ?regimen ?source
                    {{
                        ?biomarker rdfs:subClassOf oncokb_civic:Biomarker.
                        ?biomarker rdfs:subClassOf ?R1.
                            ?R1 owl:onProperty oncokb_civic:hasGene.
                            ?R1 owl:someValuesFrom oncokb_civic:{gene}.
                        ?biomarker rdfs:subClassOf ?R2.
                            ?R2 owl:onProperty oncokb_civic:hasVariant.
                            ?R2 owl:someValuesFrom oncokb_civic:{variant}.
                        ?biomarker rdfs:subClassOf ?R3.
                            ?R3 owl:onProperty oncokb_civic:hasDisease.
                            ?R3 owl:someValuesFrom oncokb_civic:{disease}.
                        ?biomarker rdfs:subClassOf ?R4.
                            ?R4 owl:onProperty oncokb_civic:evidenceSource.
                            ?R4 owl:hasValue  ?source.
                       
                           
                        ?regimen rdfs:subClassOf oncokb_civic:TherapyRegimen.
                        ?regimen rdfs:subClassOf ?R5.
                            ?R5 owl:onProperty oncokb_civic:hasEvidenceLevel{evidence_level}.
                            ?R5 owl:someValuesFrom ?biomarker.
                    }}
                """
            , error_on_undefined_entities=False)
        )
    

    return [(str(tx[0]), str(tx[1]))  for tx in treatments]

def infer_row(row):
	gene = row['Gene']
	variant = clean_variant(row['Description'])
	disease = row['TCGA_Cohort']
	patient = row['Participant_ID']
	patient_regimens=[]
	for level in ["1","2","3","4","R1","R2"]:
		# try:
		regimens = get_therapy_regimens(onto, gene,variant,disease,level)
		for reg in regimens:
			regimen, source = reg
			patient_regimens.append(pd.Series([patient, gene, variant, disease,level,source,regimen]))
		# except:
		# 	continue
	return patient_regimens


def main(ontology: str, tcga_variants: str, output: str, threads:int):
	tcga_variants = pd.read_csv(tcga_variants,sep="\t",low_memory=False)
	tcga_variants = tcga_variants.loc[tcga_variants['Variant_Type']!='Copy_Number_Alteration']
	tcga_variants = tcga_variants.loc[tcga_variants['TCGA_Cohort'] =='LUAD']
	rows = [row for index,row in tcga_variants.iterrows()]
	with multiprocessing.Pool(threads) as p:
		patient_regimens = list(tqdm.tqdm(p.imap(infer_row,  rows), total = len(rows)))
	flattened = itertools.chain.from_iterable(patient_regimens)

	all_patient_regimens = pd.concat(flattened,axis=1).T
	all_patient_regimens.columns = ['patient','gene','variant','disease','level','therapy']
	all_patient_regimens.to_csv(output,index=False)
	return

if __name__=="__main__":
	defopt.run(main)


