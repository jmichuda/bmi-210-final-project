"""Command line tool for loading CIViC data.

USAGE: python load_CIViC_data.py <OUT_PATH>

Corresponding Author:   Joseph Wakim
Affiliation:            Stanford
Date:                   March 1, 2022
"""
from typing import Optional
import sys

import pandas as pd

import api_tools as apit


def save_df(df, save_path: Optional[str]=None):
    """Save table to CSV file if a `save_path` is provided.

    Parameters
    ----------
    df : pd.DataFrame
        Table to save to CSV file
    save_path : Optional[str]
        Path at which to save data from API endpoint, relative to current
        working directory (default = None); if None, table is not saved
    """
    if save_path is not None:
        df.to_csv(save_path)


def load_data_from_civic_endpoint(url: str, save_path: Optional[str]=None):
    """Load/save data from API endpoint as a csv file

    Parameters
    ----------
    url : str
        URL to the API endpoint for plotting
    save_path : Optional[str]
        Path at which to save data from API endpoint, relative to current
        working directory (default = None); if None, table is not saved

    Returns
    -------
    pd.DataFrame
        Table containing data from API endpoint
    """
    data = apit.Endpoint(url=url)
    df = data.data_as_pandas("records")
    save_df(df, save_path)
    return df


def filter_civic_evidence_items(df_evidence: pd.DataFrame) -> pd.DataFrame:
    """Filter records in CIViC evidence items table.

    Include only "Predictive" evidence type, "Supports" evidence direction. Drop
    all records with "N/A" clinical significance.

    Parameters
    ----------
    df_evidence : pd.DataFrame
        Table of evidence items loaded from CIViC database.

    Returns
    -------
    pd.DataFrame
        Filtered table from CIViC database
    """
    df_evidence = \
        df_evidence[df_evidence["evidence_type"] == "Predictive"]
    df_evidence = \
        df_evidence[df_evidence["evidence_direction"] == "Supports"]
    df_evidence = \
        df_evidence[df_evidence["clinical_significance"] != "N/A"]
    return df_evidence


def merge_civic_variant_evidence(
    df_variant: pd.DataFrame,
    df_evidence: pd.DataFrame,
    save_path: Optional[str]=None
) -> pd.DataFrame:
    """Merge variant and evidence data from CIViC.

    Notes
    -----
    Saves only columns relevant to our ontology.

    Parameters
    ----------
    df_variant : pd.DataFrame
        Table of variants from CIViC database
    df_evidence : pd.DataFrame
        Filtered table of evidence items from CIViC database
    save_path : Optional[str]
        Path at which to save joined variant/evidence item table (default =
        None); if None, table is not saved

    Returns
    -------
    pd.DataFrame
        Table of combined variant and evidence data
    """
    df_joined = df_evidence.merge(
        df_variant, left_on="variant_id", right_on="id", how="left"
    )
    columns_of_interest = [
        "name_y", "drugs", "evidence_level", "clinical_significance", "disease"
    ]
    df_joined = df_joined.loc[:, columns_of_interest]

    for i in range(df_joined.shape[0]):
        drug_list = [drug["name"] for drug in df_joined.loc[i, "drugs"]]
        therapy_regimen = "+".join(drug_list)
        df_joined.loc[i, "therapy_regimen"] = therapy_regimen
        df_joined.loc[i, "disease"] = df_joined.loc[i, "disease"]["name"]

    df_joined.drop(columns=["drugs"], inplace=True)
    column_mappings = {
        "name_y": "variant",
        "therapy_regimen": "TherapyRegimen",
        "evidence_level": "EvidenceLevel",
        "clinical_significance": "ClinicalSignificance",
        "disease": "Disease"
    }
    df_joined.rename(columns=column_mappings, inplace=True)

    save_df(df_joined, save_path)


def main(save_path: Optional[str]=None):
    """Load variant and evidence item data from CIViC and save table to CSV.

    Parameters
    ----------
    save_path : Optional[str]
        Path at which to save joined variant/evidence item table (default =
        None); if None, table is not saved
    """
    variant_api_url = "https://civicdb.org/api/variants?count=3056"
    evidence_item_api_url = "https://civicdb.org/api/evidence_items?count=8579"

    df_variant = load_data_from_civic_endpoint(url=variant_api_url)
    df_evidence = filter_civic_evidence_items(
        load_data_from_civic_endpoint(url=evidence_item_api_url)
    )
    merge_civic_variant_evidence(df_variant, df_evidence, save_path)


if __name__ == "__main__":
    """Load data from CIViC.
    """
    main(save_path=sys.argv[1])
