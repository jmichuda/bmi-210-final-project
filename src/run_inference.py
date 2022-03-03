import pandas as pd
from owlready2 import *
import defopt
from src.generate_data import clean_variant

def get_therapy_regimens(onto, gene, variant, disease, evidence_level):
	with onto:
		treatments = list(
			default_world.sparql(
				f"""
					SELECT distinct ?regimen
					{{
						?biomarker rdfs:subClassOf oncokb:Biomarker.
						?biomarker rdfs:subClassOf ?R1.
							?R1 owl:onProperty oncokb:hasGene.
							?R1 owl:someValuesFrom oncokb:{gene}.
						?biomarker rdfs:subClassOf ?R2.
							?R2 owl:onProperty oncokb:hasVariant.
							?R2 owl:someValuesFrom oncokb:{variant}.
						?biomarker rdfs:subClassOf ?R3.
							?R3 owl:onProperty oncokb:hasDisease.
							?R3 owl:someValuesFrom oncokb:{disease}.
						?regimen rdfs:subClassOf oncokb:TherapyRegimen.
						?regimen rdfs:subClassOf ?R4.
							?R4 owl:onProperty oncokb:hasEvidenceLevel{evidence_level}.
							?R4 owl:someValuesFrom ?biomarker.
					}}
				"""
			,error_on_undefined_entities=False)
		)
	return [str(tx[0]) for tx in treatments]

def main(ontology: str, tcga_variants: str, output: str):
	onto = get_ontology(ontology).load()
	tcga_variants = pd.read_csv(tcga_variants,sep="\t",low_memory=False)
	tcga_variants = tcga_variants.loc[tcga_variants['Variant_Type']!='Copy_Number_Alteration']
	tcga_variants = tcga_variants.loc[tcga_variants.Gene == 'ATM'].sample(300)

	patient_regimens = []
	for index, row in tcga_variants.iterrows():
		gene = row['Gene']
		variant = clean_variant(row['Description'])
		disease = row['TCGA_Cohort']
		patient = row['Participant_ID']
		
		regimens = get_therapy_regimens(onto, gene,variant,disease,"3")
		print(gene,variant,disease,regimens)
		if len(regimens) > 0:
			patient_regimens.append(pd.Series([patient,gene,variant,disease,regimens]))
	all_patient_regimens = pd.concat(patient_regimens,axis=1).T
	all_patient_regimens.to_csv(output)
	return

if __name__=="__main__":
	defopt.run(main)