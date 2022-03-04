"""Load patient variant data from cBioPortal.

USAGE: python load_cBioPortal_data.py <SAVE_PATH> <OPTIONAL_CANCER_GENE_PATH>

Corresponding Author:   Joseph Wakim
Affiliation:            Stanford
Date:                   March 1, 2022
"""
import sys
from typing import Optional

import pandas as pd
from bravado.client import SwaggerClient

import api_tools as apit
import access_cBioPortal as cbio


def load_mutations(
    client: SwaggerClient,
    *,
    save_path: str,
    gene_path: Optional[str] = None
):
    """Load patient mutation data from cBioPortal.

    Parameters
    ----------
    client : SwaggerClient
        Swagger API client accessing cBioPortal
    save_path : str
        Path at which to save table of mutation data
    gene_path : Optional[str]
        Path to TSV file containing gene set (name in column 0, entrez gene ID
        in column 1) to which mutation data will be filtered (default = None);
        if None, mutation data will not be filtered and all variants will be
        maintained
    """
    if gene_path is not None:
        risk_genes = pd.read_csv(gene_path, sep="\t", header=0, usecols=[0, 1])
        gene_list = risk_genes.iloc[:, 1].to_numpy().astype(int)
    else:
        gene_list = None

    cbio.get_all_from_bioportal_endpoint(
        client, "Mutations", out_path=save_path, gene_list=gene_list
    )


def main(save_path: str, gene_path: Optional[str] = None):
    """Save patient mutation data from cBioPortal.

    Parameters
    ----------
    save_path : str
        Path at which to save table of mutation data
    gene_path : Optional[str]
        Path to TSV file containing gene set (name in column 0, entrez gene ID
        in column 1) to which mutation data will be filtered (default = None);
        if None, mutation data will not be filtered and all variants will be
        maintained
    """
    cBioPortal_api_url = "https://www.cbioportal.org/api/api-docs"
    client = apit.get_swagger_api_client(url=cBioPortal_api_url)
    load_mutations(client=client, save_path=save_path, gene_path=gene_path)


if __name__ == "__main__":
    sys_args = sys.argv
    if len(sys_args) == 2:
        main(save_path=sys_args[1])
    elif len(sys_args) == 3:
        main(save_path=sys_args[1], gene_path=sys_args[2])
    else:
        raise ValueError(
            "Please specify an output path at which to save mutation data. "
            "Optionally, you may also specify a path to a gene list to include."
        )
