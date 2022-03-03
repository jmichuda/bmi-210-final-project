"""Sample patients from the CSV file of cBioPortal gene variants.

USAGE: python sample_patients.py <MUTATION_CSV_FILE> <OUTPUT_CSV_FILE> <OPTIONAL_NUMBER_OF_PATIENTS>

Corresponding Author:   Joseph Wakim
Affiliation:            Stanford
Date:                   March 1, 2022
"""
import csv
import sys
from typing import Optional

import numpy as np
import pandas as pd


def sample_patient_records(
    csv_path: str, out_path: str, patient_column: Optional[int] = 24,
    n_patients: Optional[int] = 1000, n_row_skip: Optional[int] = 2
):
    """Load a subset of patient records from a CSV file

    Parameters
    ----------
    csv_path : str
        Path to CSV file containing variant data
    out_path : str
        Path at which to save CSV file containing variant data from a subset of
        patients
    patient_column : Optional[int]
        Column index identifying patient identifiers in CSV file of variant data
        (default = 24)
    n_patients : Optional[int]
        Number of patients to sample from variant data CSV file (default = 1000)
    n_row_skip : Optional[int]
        Number of rows corresponding to header information in variant CSV file
        (default = 2)
    """
    count = 0
    patient_set = set()

    with open(csv_path, "r") as read_file:
        with open(out_path, "w") as write_file:

            reader = csv.reader(read_file)
            writer = csv.writer(write_file, delimiter=',')

            i = 0
            for row in reader:
                i += 1
                if i < n_row_skip:
                    continue
                elif i == n_row_skip:
                    writer.writerow(row)
                else:
                    if row[patient_column] in patient_set:
                        writer.writerow(row)
                    else:
                        patient_set.add(row[patient_column])
                        count += 1
                        if count > n_patients:
                            print("Done")
                            return
                        writer.writerow(row)


def get_patient_mutation_data(
        csv_path: str, patient_id: str, patient_column: int = 24,
        n_row_skip: Optional[int] = 1, ind_header_row: Optional[int] = 1
) -> pd.DataFrame:
    """Load a table of mutation data associated with a single patient.

    Parameters
    ----------
    csv_path : str
        Path to CSV file containing variant data
    patient_id : str
        Patient ID for which the mutation data is requested
    patient_column : Optional[int]
        Column index identifying patient identifiers in CSV file of variant data
        (default = 24)
    n_row_skip : Optional[int]
        Number of rows corresponding to header information in variant CSV file
        (default = 2)
    ind_header_row : Optional[int]
        Index of the header row, starting from 0 (default = 1)

    Returns
    -------
    pd.DataFrame
        Table of mutation data associated with the specific patient of interest
    """
    with open(csv_path, "r") as read_file:
        reader = csv.reader(read_file)
        searching = True
        i = 0
        for row in reader:
            i += 1
            if i <= n_row_skip:
                continue

            if (i-1) == ind_header_row:
                headers = np.array(row)
                df = pd.DataFrame(columns=headers)

            if row[patient_column] == patient_id:
                if searching:
                    print(f"Patient '{patient_id}' found!")
                    searching = False
                new_row = pd.DataFrame.from_records(
                    [{headers[ind]: row[ind] for ind in range(len(row))}]
                )
                df = pd.concat([df, new_row])
            else:
                if not searching:
                    return df
    return df


if __name__ == "__main__":
    """Extract the records of interest.
    """
    csv_file = sys.argv[1]
    out_file = sys.argv[2]
    if len(sys.argv) > 3:
        n_patients = int(sys.argv[3])

    sample_patient_records(
        csv_path=csv_file, out_path=out_file, n_patients=n_patients
    )
