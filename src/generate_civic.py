import pandas as pd
import defopt
import api_tools as apit

def generate_evidence(output:str):
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

    gene_list = pd.read_csv("CancerGeneList.tsv", sep="\t", header=0, usecols=[0, 1])
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

    df_variant_evidence_filtered.to_csv(output)

if __name__=="__main__":
    defopt.run(main)
