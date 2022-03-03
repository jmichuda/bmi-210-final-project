"""Access cBioPortal in a format compatible with Snakefile.

This module accesses tools in cBioPortal in a format that works with our
pipeline for automatic ontology generation.

Corresponding Author:   Joseph Wakim
Affiliation:            Stanford
Date:                   February 27, 2022
"""
from typing import Optional, Sequence

import pandas as pd
import numpy as np
from bravado.client import SwaggerClient

import api_tools as apit


endpoint_set = {
    "Molecular Profiles", "Studies", "Sample Lists", "Cancer Types",
    "Patients", "Samples", "Clinical Data", "Mutations"
}

simple_api_methods = {
    "Molecular Profiles": "getAllMolecularProfilesUsingGET",
    "Studies": "getAllStudiesUsingGET",
    "Sample Lists": "getAllSampleListsUsingGET",
    "Cancer Types": "getAllCancerTypesUsingGET",
    "Patients": "getAllPatientsUsingGET"
}


def get_all_from_bioportal_endpoint(
    client: SwaggerClient, endpoint: str, out_path: Optional[str] = None,
    **kwargs
) -> pd.DataFrame:
    """Load all data from an endpoint of cBioPortal

    Wrapper for cBioPortal convenience functions. Returns a Pandas DataFrame,
    which we can use to populate our ontology.

    Parameters
    ----------
    client : SwaggerClient
        Swagger client accessing the cBioPortal database
    endpoint : str
        Name of the cBioPortal endpoint from which to load data
    out_path : Optional[str]
        Path relative to current working directory at which to save output of
        mutations data.

    Returns
    -------
    pd.DataFrame
        Table of data at specified endpoint.
    """
    if endpoint not in endpoint_set:
        raise ValueError(
            f"Specified endpoint `{endpoint}` not found in cBioPortal."
        )
    obj = getattr(client, endpoint)

    if endpoint in simple_api_methods.keys():
        all_records = getattr(obj, simple_api_methods[endpoint])().result()
        df = apit.swagger_endpoint_data_to_df(all_records)
        if out_path is not None:
            df.to_csv(out_path)
        return df

    else:
        if endpoint == "Samples":
            df = get_all_samples(client)
            if out_path is not None:
                df.to_csv(out_path)
            return df
        elif endpoint == "Clinical Data":
            df = get_patient_cancer_types(client)
            if out_path is not None:
                df.to_csv(out_path)
            return df
        elif endpoint == "Mutations":
            if out_path is None:
                raise ValueError(
                    "Specify an output directory to save mutation data."
                )
            get_mutations(client, out_path=out_path, **kwargs)


def get_unique_study_ids(client: SwaggerClient):
    """Get unique study IDs from cBioPortal

    Parameters
    ----------
    client : SwaggerClient
        Swagger Client accessing the  cBioPortal database

    Returns
    -------
    Sequence[str]
        Array of unique study IDs
    """
    studies = get_all_from_bioportal_endpoint(client, "Studies")
    return np.unique(studies.studyId.to_numpy())


def get_all_samples(
    client: SwaggerClient, unique_study_ids: Optional[Sequence[str]] = None
) -> pd.DataFrame:
    """Get all samples from all studies.

    Parameters
    ----------
    client : SwaggerClient
        Swagger client accessing the cBioPortal database
    unique_study_ids : Optional[Sequence[str]]
        List of unique study IDs (default = None). If argument is not specified,
        then a list of all unique study IDs is determined.

    Returns
    -------
    pd.DataFrame
        Table of sample data
    """
    obj = client.Samples

    if unique_study_ids is None:
        unique_study_ids = get_unique_study_ids(client)

    df_samples = pd.DataFrame()
    n_studies = len(unique_study_ids)

    for i, study in enumerate(unique_study_ids):
        if (i+1) % 100 == 0:
            print(f"Parsing Study {i+1} of {n_studies}")

        all_samples = obj.getAllSamplesInStudyUsingGET(studyId=study).result()

        if i == 0:
            headers = dir(all_samples[0])

        df_temp = apit.swagger_endpoint_data_to_df(all_samples, headers)
        df_samples = pd.concat([df_samples, df_temp])

    return df_samples


def get_patient_cancer_types(
    client: SwaggerClient, unique_study_ids: Optional[Sequence[str]] = None
) -> pd.DataFrame:
    """Get all samples from all studies.

    Parameters
    ----------
    client : SwaggerClient
        Swagger client accessing the cBioPortal database
    unique_study_ids : Optional[Sequence[str]]
        List of unique study IDs (default = None). If argument is not specified,
        then a list of all unique study IDs is determined.

    Returns
    -------
    pd.DataFrame
        Table of patient cancer types
    """
    obj = client.Clinical_Data

    if unique_study_ids is None:
        unique_study_ids = get_unique_study_ids(client)

    df_clinical_data = pd.DataFrame()
    n_studies = len(unique_study_ids)

    for i, study in enumerate(unique_study_ids):
        if (i+1) % 100 == 0:
            print(f"Parsing Study {i+1} of {n_studies}")

        all_clinical_data = obj.getAllClinicalDataInStudyUsingGET(
            studyId = study
        ).result()

        if i == 0:
            headers = dir(all_clinical_data[0])

        df_temp = apit.swagger_endpoint_data_to_df(all_clinical_data, headers)
        df_clinical_data = pd.concat([df_clinical_data, df_temp])

    df_clinical_data = df_clinical_data[
        df_clinical_data["clinicalAttributeId"] == "CANCER_TYPE"
    ].drop_duplicates()

    relevant_columns = ["patientId", "sampleId", "studyId", "value"]
    df_patient_cancer_types = df_clinical_data.loc[:, relevant_columns]
    df_patient_cancer_types.rename(columns={'value':'name'}, inplace=True)

    df_cancer_types = get_all_from_bioportal_endpoint(client, "Cancer Types")
    df_patient_cancer_types = df_patient_cancer_types.merge(
        df_cancer_types, on="name", how="left"
    )

    return df_patient_cancer_types


def get_mutations(
    client: SwaggerClient,
    out_path: str,
    *,
    gene_list: Optional[Sequence[str]]=None,
    df_molecular_profiles: Optional[pd.DataFrame]=None,
    df_sample_lists: Optional[pd.DataFrame]=None
):
    """Get all mutations in cBioPortal corresponding to a gene in gene_list.

    Parameters
    ----------
    client : SwaggerClient
        Swagger client accessing the cBioPortal database
    out_path : str
        Path relative to current working directory at which to save CSV file of
        variants
    gene_list : Optional[Sequence[str]]
        List of genes to which variants will be filtered (default = None); if
        the argument is not specified, no variants will be filtered
    df_molecular_profiles : Optional[pd.DataFrame]
        Table of molecular profiles
    df_sample_lists : Optional[pd.DataFrame]
        Table of sample lists
    """
    mutations = client.Mutations
    df_temp = pd.DataFrame(list())
    df_temp.to_csv(out_path)
    empty = True

    if df_molecular_profiles is None:
        df_molecular_profiles = get_all_from_bioportal_endpoint(
            client, "Molecular Profiles"
        )

    if df_sample_lists is None:
        df_sample_lists = get_all_from_bioportal_endpoint(
            client, "Sample Lists"
        )

    df_joined = df_molecular_profiles.merge(
        df_sample_lists,
        on="studyId", how="outer"
    ).loc[:, ["molecularProfileId", "sampleListId"]]
    n_iters = df_joined.shape[0]

    for i in range(n_iters):

        if (i+1) % 100 == 0:
            print(f"Getting mutation list {i+1} of {n_iters}")

        try:
            all_mutations =\
                mutations.getMutationsInMolecularProfileBySampleListIdUsingGET(
                    molecularProfileId = df_joined.iloc[i, 0],
                    sampleListId = df_joined.iloc[i, 1]
                ).result()
        except Exception:
            continue

        if len(all_mutations) == 0:
            continue

        df_temp = apit.swagger_endpoint_data_to_df(all_mutations)
        if gene_list is not None:
            df_temp = df_temp.query(
                "proteinChange in @gene_list | gene in @gene_list | entrezGeneId in @gene_list"
            )

        if empty:
            headers = list(df_temp.columns)
            df_temp.to_csv(out_path, mode="a", index=False, header=headers)
            empty = False

        else:
            df_temp.to_csv(out_path, mode="a", index=False, header=False)
