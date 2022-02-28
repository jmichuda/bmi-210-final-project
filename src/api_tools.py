"""Tools for accessing API clients using Swagger.

Corresponding Author:   Joseph Wakim
Affiliation:            Stanford University
Date:                   February 26, 2022
"""

from typing import Optional, List, Mapping, Sequence

import pandas as pd
import requests
from bravado.client import SwaggerClient


class Endpoint:
    """Interface with API endpoint using Python.

    Serves as a wrapper for methods in `requests` package.
    """
    def __init__(self, url: str):
        """Initialize the object representation of the API endpoint.

        Parameters
        ----------
        url : str
            URL to the API endpoint
        """
        self.url = url
        self.requests_obj = requests.get(url)

        if self.requests_obj.status_code == 404:
            raise ValueError("The specified API endpoint url was not found.")

        self.headers = self.get_headers()

    def get_headers(self):
        """Return headers of the API endpoint.

        Returns
        -------
        list
            Headers at the API endpoint
        """
        return set(self.requests_obj.json().keys())

    def check_header(self, header: str):
        """Check that a requested header exists in the API endpoint.

         Parameters
        ----------
        header : str
            Name of the header to check
        """
        if header not in self.headers:
            raise ValueError(f"The header `{header}` was not found.")

    def access_data_at_header(self, header: str) -> List[dict]:
        """Access data at the API endpoint under a specified header.

        Parameters
        ----------
        header : str
            Name of the header to access

        Returns
        -------
        List[dict]
            The data associated with the endpoint under the specified header,
            which takes the form of a list of dictionaries
        """
        self.check_header(header)
        return self.requests_obj.json()[header]

    def data_as_pandas(
        self,
        headers: str or Sequence[str],
        search_keys: Optional[Sequence[str]] = [],
        search_vals: Optional[Sequence[str]] = []
    ) -> pd.DataFrame:
        """Tabulate data at specified level of endpoint's JSON data structure.

        Notes
        -----
        The JSON data structure returned from the GET call to the API Endpoint
        is represented by a list of dictionaries. Each dictionary in this list
        can have a key whose corresponding value is also a list of dictionaries.
        This nested "list-of-dicts" format makes it difficult to parse our data.

        This method locates one of the nested lists-of-dicts and returns its
        values as a pandas DataFrame. To do this, the user must specify an
        ordered sequence of `headers`, each of which corresponds to a key whose
        value is a list of dictionaries.

        Each time a header is located, a list of dictionaries (all of which have
        common keys) is identified. We may be interested in a dictionary in this
        list which has a key whose value is also a list-of-dicts (hense the
        nested list-of-dicts pattern). The `search_keys` and `search_vals`
        arguments represent ordered sequences of keys and corrresponding values
        with which we locate the next-layer list-of-dicts that we are interested
        in tabulating.

        Pseudocode
        ----------
        dict = JSON data from API endpoint as list-of-dicts
        for i in range(number_of_headers):
            list_of_dicts = dict[headers[i]]
            if i == number_of_headers - 1:
                return tabular representation of list_of_dicts
            else:
                for row in list-of-dict:
                    if row[search_keys[i]] == search_values[i]:
                    dict = row

        Parameters
        ----------
        headers : str or Sequence[str]
            Sequence of nested headers pointing to data for tabulation, in the
            order by which they appear in the JSON data structure; if entered
            as a string, then the list of dictionaries for that string is
            returned as a data frame.
        search_keys : Optional[Sequence[str]]
            Sequence of keys on which to search for a specific dictionary in a
            list of dictionaries (default = [])
        search_vals : Optional[Sequence[str]]
            Sequence of values for respective search key on which  to isolate a
            specific dictionary in a list of dictionaries (default = [])

        Returns
        -------
        pd.DataFrame
            Table representing data at some position of the endpoint.
        """
        if isinstance(headers, str):
            list_of_dicts = self.access_data_at_header(headers)

        else:
            n_headers = len(headers)

            if n_headers == 0:
                return pd.DataFrame()

            list_of_dicts = self.access_data_at_header(headers[0])
            for i in range(n_headers-1):
                for row in list_of_dicts:
                    if row[search_keys[i]] == search_vals[i]:
                        list_of_dicts = row[headers[i+1]]
                        break

        if len(list_of_dicts) == 0:
            raise ValueError(f"No data stored at endpoint header sequence.")

        keys = list(list_of_dicts[0].keys())
        return list_dict_to_pandas(keys, list_of_dicts)


def list_dict_to_pandas(
        keys: List[str], dictionaries: List[dict]
) -> pd.DataFrame:
    """Return pandas DataFrame representation of a list of dictionaries.

    Notes
    -----
    All dictionaries should represent a single record and should have a common
    set of keys. The dictionaries represent JSON-formatted structures. This
    function is relevant for representing data at an API endpoint as a pandas
    DataFrame.

    Parameters
    ----------
    keys : List[str]
        List of keys common to each dicitonary in `dictionaries`
    dictionaries : dict
        List of dictionaries representing individual records to be tabulated
    """
    records = [
        {
            key: dictionary[key] for key in keys
        } for dictionary in dictionaries
    ]
    return pd.DataFrame.from_records(records)


def get_swagger_api_client(
        url: str,
        *,
        val_requests: Optional[bool] = False,
        val_responses: Optional[bool] = False,
        val_swagger_spec: Optional[bool] = False
) -> SwaggerClient:
    """Generate a Swagger client for the API from a URL.

    This function wraps the `bravo.client.SwaggerClient.from_url` for
    our particular use case.

    Parameters
    ----------
    url : str
        Path to the API endpoint
    val_requests : Optional[bool]
        Check if requests sent to the REST API meet validation rules
        set forth by the API (default = False)
    val_responses : Optional[bool]
        Check if responses obtained from the REST API meet validation
        rules set forth by the API (default = False)
    val_swagger_spec : Optional[bool]
        Check if the responses obtained from the REST API match the
        schema (or other structure) defined in the metadata of the API
        (default = False)

    Returns
    -------
    SwaggerClient
        Swagger client object for the specified API
    """
    return SwaggerClient.from_url(
        spec_url = url,
        config = {
            "validate_requests": val_requests,
            "validate_responses": val_responses,
            "validate_swagger_spec": val_swagger_spec
        }
    )


def swagger_endpoint_data_to_df(
        data: List[Mapping],
        headers: Optional[List[str]] = None
) -> pd.DataFrame:
    """Load results from cBioPortal API endpoints to pandas DataFrame.

    Parameters
    ----------
    data : List[Mapping]
        List of dictionaries at an endpoint of the API
    headers : Optional[List[str]]
        Headers of the data contained at the endpoint; if None, headers are read
        from the first endpoint in the list

    returns
    -------
    pd.DataFrame
        Table representation of data at endpoint
    """
    if headers is None:
        headers = dir(data[0])
    return list_dict_to_pandas(headers, data)
