"""Query for `TherapyRegimen` associated with mutated genes.

Corresponding Authors:  Joseph Wakim, Vinodh Rajapakse
Affiliation:            Stanford
Date:                   March 3, 2022

USAGE: python query_therapy_regimen.py <ONTOLOGY_PATH> <MUTATION_TABLE_PATH> <GENE_LIST_PATH> <PATIENT_NAME> <REGIMEN_SAVE_PATH>
"""
import csv
import sys
import os
import re
from typing import List, Optional

import numpy as np
import pandas as pd
import owlready2 as or2

from sample_patients import get_patient_mutation_data


def load_ongology(path: str) -> or2.namespace.Ontology:
    """Load ontology from specified file path.

    Parameters
    ----------
    path : str
        Path to the ontology in the local directory

    Returns
    -------
    or2.namespace.Ontology
        Object representing the ontology in owlready2
    """
    return or2.get_ontology(path).load()


def get_mutation_data(
        csv_path: str, gene_path: str,
        df_mutations: Optional[pd.DataFrame] = None,
        df_genes: Optional[pd.DataFrame] = None
) -> pd.DataFrame:
    """Get cBioPortal mutation data and add gene labels, if available.

    Parameters
    ----------
    csv_path : str
        Path to csv file containing mutation data (if `df_mutation` is not None,
        this path has no effect)
    gene_path : str
        Path to table of genes associated with cancer and their entrez IDs (if
        `df_genes` is not None, this path has no effect)
    df_mutations: Optional[pd.DataFrame]
        If a mutation table is loaded as a Pandas DataFrame, you can enter the
        DataFrame object under this argument to avoid reloading the mutation
        data (default = None)
    df_genes : Optional[pd.DataFrame]
        If a gene table is already loaded as a Pandas DataFrame, you can enter
        the DataFrame object under this argument to avoid reloading the gene
        data (default = None)

    Returns
    -------
    pd.DataFrame
        Table of patients with mutations and associated genes
    """
    columns_of_interest = ["patientId", "proteinChange", "entrezGeneId"]
    if df_mutations is None:
        df_mutations = pd.read_csv(csv_path)
    df_mutations = df_mutations.loc[:, columns_of_interest]
    if df_genes is None:
        df_genes = pd.read_csv(gene_path, sep="\t", usecols=[0, 1])
    df_genes["Entrez_Id"] = \
        df_genes["Entrez_Id"].dropna().astype("int").astype("str")
    df = df_mutations.merge(
        df_genes, left_on="entrezGeneId", right_on="Entrez_Id", how="left"
    ).drop("Entrez_Id", axis=1)
    return df


def get_therapy_given_gene(
    ontology: or2.namespace.Ontology, gene: str, evidence_level: str
) -> List[str]:
    """Query ontology for therapy regimen associated with mutation in gene.

    Parameters
    ----------
    ontology : or2.namespace.Ontology
        OWL ontology from which to query therapy regimens matching mutated gene
    gene : str
        Gene for which to query therapy regimen
    evidence_level : str
        Evidence level of therapy regemen to query for in ontology

    Returns
    -------
    List[str]
        List of therapy regimens associated with gene of interest
    """
    with ontology:
        therapy_list = list(
            or2.default_world.sparql(
                f"""
                SELECT distinct ?regimen
                {{
                    ?biomarker rdfs:subClassOf oncokb:Biomarker.
                    ?biomarker rdfs:subClassOf ?r1.
                    ?r1 owl:onProperty oncokb:hasGene.
                    ?r1 owl:someValuesFrom oncokb:{re.escape(gene)}.
                    
                    ?regimen rdfs:subClassOf oncokb:TherapyRegimen.
                    ?regimen rdfs:subClassOf ?r2.
                    ?r2 owl:onProperty oncokb:hasEvidenceLevel{evidence_level}.
                    ?r2 owl:someValuesFrom ?biomarker.
                }}
                """
            )
        )
    therapy_list = [str(gene[0]) for gene in therapy_list]
    return therapy_list


def get_therapy_given_gene_variant(
    ontology: or2.namespace.Ontology, gene: str, variant: str,
    evidence_level: str
) -> List[str]:
    """Query ontology for therapy matching gene and variant.

    Given a patient's gene and variant, query the ontology for all therapy
    regimens.

    Parameters
    ----------
    ontology : or2.namespace.Ontology
        OWL ontology from which to query therapy regimens matching mutated gene
    gene : str
        Gene for which to query therapy regimen
    variant : str
        Specific variant of gene for which to query therapy regiment
    evidence_level : str
        Evidence level of therapy regemen to query for in ontology

    Returns
    -------
    List[str]
        List of therapy regimens associated with the gene and variant
    """
    with ontology:
        treatments = list(
            or2.default_world.sparql(
                f"""
                SELECT distinct ?regimen
                {{
                    ?biomarker rdfs:subClassOf oncokb:Biomarker.
                    ?biomarker rdfs:subClassOf ?R1.
                    ?R1 owl:onProperty oncokb:hasGene.
                    ?R1 owl:someValuesFrom oncokb:{re.escape(gene)}.
                    
                    ?biomarker rdfs:subClassOf ?R2.
                    ?R2 owl:onProperty oncokb:hasVariant.
                    ?R2 owl:someValuesFrom oncokb:{re.escape(variant)}.
                    
                    ?regimen rdfs:subClassOf oncokb:TherapyRegimen.
                    ?regimen rdfs:subClassOf ?R3.
                    ?R3 owl:onProperty oncokb:hasEvidenceLevel{evidence_level}.
                    ?R3 owl:someValuesFrom ?biomarker.
                }}
                """
            )
        )
    return [str(tx[0]) for tx in treatments]


def get_therapy_given_gene_variant_disease(
    ontology: or2.namespace.Ontology, gene: str, variant: str, disease: str,
    evidence_level: str
) -> List[str]:
    """Query ontology for therapy matching gene, variant, disease, and evidence.

    Given a patient's gene, variant, and disease, query the ontology for
    therapy regimens with a specified evidence level.

    Parameters
    ----------
    ontology : or2.namespace.Ontology
        OWL ontology from which to query therapy regimens matching mutated gene
    gene : str
        Gene for which to query therapy regimen
    variant : str
        Specific variant of gene for which to query therapy regiment
    disease : str
        Specific disease associated with gene and variant for which to query the
        ontology
    evidence_level : str
        Evidence level of therapy regemen to query for in ontology

    Returns
    -------
    List[str]
        List of therapy regimens associated with the gene, variant, and disease
        matching records in the ontology
    """
    with ontology:
        treatments = list(
            or2.default_world.sparql(
                f"""
                SELECT distinct ?regimen
                {{
                    ?biomarker rdfs:subClassOf oncokb:Biomarker.
                    ?biomarker rdfs:subClassOf ?R1.
                    ?R1 owl:onProperty oncokb:hasGene.
                    ?R1 owl:someValuesFrom oncokb:{re.escape(gene)}.
                    
                    ?biomarker rdfs:subClassOf ?R2.
                    ?R2 owl:onProperty oncokb:hasVariant.
                    ?R2 owl:someValuesFrom oncokb:{re.escape(variant)}.
                    
                    ?biomarker rdfs:subClassOf ?R3.
                    ?R3 owl:onProperty oncokb:hasDisease.
                    ?R3 owl:someValuesFrom oncokb:{re.escape(disease)}.
 
                    ?regimen rdfs:subClassOf oncokb:TherapyRegimen.
                    ?regimen rdfs:subClassOf ?R4.
                    ?R4 owl:onProperty oncokb:hasEvidenceLevel{evidence_level}.
                    ?R4 owl:someValuesFrom ?biomarker.
                }}
                """
            )
        )
    return [str(tx[0]) for tx in treatments]


def main(
    ontology_path: str,
    mutation_path: str,
    gene_list_path: str,
    patient_id: str,
    regimen_save_path: str
):
    """Save potential therapy regimen associated with mutated gene in patient.

    Parameters
    ----------
    ontology_path : str
        Path to our custom ontology combining biomedical data sources
    mutation_path : str
        Path to patient mutation data extracted from cBioPortal
    gene_list_path : str
        Path to table of genes associated with cancer and their entrez IDs
    patient_id : str
        Patient identifier matching record in cBioPortal mutation table
    regimen_save_path : str
        Path at which to save list of regimens associated with each affected
        gene in the patient
    """
    ontology = load_ongology(ontology_path)
    df_mutation = get_patient_mutation_data(
        mutation_path, patient_id
    )
    mutation_data = get_mutation_data(
        mutation_path, gene_list_path, df_mutations=df_mutation
    )
    genes = mutation_data["Gene_Symbol"].to_numpy()
    variants = mutation_data["proteinChange"].to_numpy()
    evidence_levels = ["1", "2", "3", "4", "R1", "R2"]

    # Clear file
    if os.path.isfile(regimen_save_path):
        file_temp = open(regimen_save_path, 'r+')
        file_temp.truncate(0)
        file_temp.close()

    # Write therapy regimen to file
    with open(regimen_save_path, "a") as f:
        writer = csv.writer(f)
        for i, gene in enumerate(genes):
            variant = re.escape(variants[i])
            for evidence in evidence_levels:
                try:
                    therapy_list = get_therapy_given_gene_variant(
                        ontology, gene, variant, evidence
                    )
                    writer.writerow(
                        [gene, variant, f"level{evidence}"] + therapy_list
                    )
                except Exception:
                    continue


if __name__ == "__main__":
    """Run the module as a command line tool.
    
    Creates a file with therapy regimen for specified user.
    """
    onto_path = sys.argv[1]
    mut_tbl_path = sys.argv[2]
    gene_lst_path = sys.argv[3]
    pat_id = sys.argv[4]
    therapy_save_path = sys.argv[5]
    main(
        onto_path, mut_tbl_path, gene_lst_path, pat_id, therapy_save_path
    )
