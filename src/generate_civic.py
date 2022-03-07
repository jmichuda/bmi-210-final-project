import pandas as pd
import defopt
import src.api_tools as apit
from src.generate_data import oncotree
import json
from thefuzz import process
from thefuzz import fuzz


def generate_oncotree_mapping():
    tree = oncotree()
    node = tree['TISSUE']
    stack = [node]
    mapping = []
    while len(stack) > 0:
        node = stack.pop()
        for key, child in node['children'].items():
            stack.append(child)
        parent = node['parent']
        if parent == "TISSUE":
            parent = "Disease"
        if node['code'] != "TISSUE":
            if "NCI" in node['externalReferences'].keys():
                for umls in node['externalReferences']['NCI']:
                    mapping.append([umls, node['code'], node['name']])
    return pd.DataFrame(mapping, columns = ['NCI','oncotree','name'])


def add_oncotree_code(df_variant_evidence_filtered, mapping_path):
    if mapping_path is None:
        oncotree = generate_oncotree_mapping()
        oncotree_mapping = [ ]
        for disease in df_variant_evidence_filtered["Disease"].drop_duplicates():
            oncotree_name, score, score2 = process.extractOne(disease, oncotree.name, scorer=fuzz.token_sort_ratio)
            oncotree_code = oncotree.set_index('name').loc[oncotree_name,'oncotree']
            oncotree_mapping.append([disease,oncotree_name,score,score2,oncotree_code])
            
        full_mapping = pd.DataFrame(oncotree_mapping, columns = ['disease','oncotree_name','score','score2','oncotree'])
        mapper = full_mapping.set_index('disease')['oncotree']
    else:
        full_mapping = pd.read_csv(mapping_path)
        mapper = full_mapping.set_index('DO disease')['New']
    df_variant_evidence_filtered['oncotree'] = df_variant_evidence_filtered.Disease.map(mapper)
    return df_variant_evidence_filtered


def main(output:str, mapping_path: str):
    civic_gene = apit.Endpoint(url="https://civicdb.org/api/genes?count=238")
    df_civic_gene = civic_gene.data_as_pandas("records")

    civic_variants = apit.Endpoint(url="https://civicdb.org/api/variants?count=3056")
    df_civic_variants = civic_variants.data_as_pandas("records")

    civic_evidence = apit.Endpoint(url="https://civicdb.org/api/evidence_items?count=8579")
    df_civic_evidence = civic_evidence.data_as_pandas("records")

    df_civic_evidence_filtered = df_civic_evidence[df_civic_evidence["evidence_type"] == "Predictive"]
    df_civic_evidence_filtered = df_civic_evidence_filtered[df_civic_evidence_filtered["evidence_direction"] == "Supports"]
    df_civic_evidence_filtered = df_civic_evidence_filtered[df_civic_evidence_filtered["clinical_significance"] != "N/A"]

    df_variant_evidence = df_civic_evidence_filtered.merge(df_civic_variants, left_on="variant_id", right_on="id", how="left")


    columns_of_interest = ["name_y", "drugs", "evidence_level", "clinical_significance", "disease", "entrez_id"]
    df_variant_evidence_filtered = df_variant_evidence.loc[:, columns_of_interest]

    for i in range(df_variant_evidence_filtered.shape[0]):
        drug_list = [drug["name"] for drug in df_variant_evidence_filtered.loc[i, "drugs"]]
        therapy_regimen = "+".join(drug_list)
        df_variant_evidence_filtered.loc[i, "therapy_regimen"] = therapy_regimen
        df_variant_evidence_filtered.loc[i, "disease"] = df_variant_evidence_filtered.loc[i, "disease"]["name"]


    gene_list = pd.read_csv("source_data/CancerGeneList.tsv", sep="\t", header=0, usecols=[0, 1])
    df_variant_evidence_filtered = df_variant_evidence_filtered.merge(gene_list, left_on="entrez_id", right_on="Entrez_Id", how="left")
    df_variant_evidence_filtered = df_variant_evidence_filtered.drop(columns=["Entrez_Id", "drugs"])

    column_mappings = {
        "name_y": "variant",
        "therapy_regimen": "TherapyRegimen",
        "evidence_level": "EvidenceLevel",
        "clinical_significance": "ClinicalSignificance",
        "disease": "Disease",
        "Gene_Symbol": "Gene"
    }
    df_variant_evidence_filtered.rename(columns=column_mappings, inplace=True)

    df_variant_evidence_filtered = add_oncotree_code(df_variant_evidence_filtered, mapping_path)
    df_variant_evidence_filtered.to_csv(output)

if __name__=="__main__":
    defopt.run(main)
