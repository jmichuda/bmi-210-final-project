"""Sample patients from the CSV file of cBioPortal gene variants.

USAGE: python sample_patients.py <MUTATION_CSV_FILE> <OUTPUT_CSV_FILE> <OPTIONAL_NUMBER_OF_PATIENTS>

Corresponding Author:   Joseph Wakim
Affiliation:            Stanford
Date:                   March 1, 2022
"""
import csv
import sys


def sample_patient_records(
    csv_path: str, out_path: str, patient_column: int = 24,
    n_patients: int = 1000, n_row_skip: int = 2
):
    """Load a subset of patient records from a CSV file

    Parameters
    ----------
    csv_path : str
        Path to CSV file containing variant data
    out_path : str
        Path at which to save CSV file containing variant data from a subset of
        patients
    patient_column : int
        Column index identifying patient identifiers in CSV file of variant data
    n_patients : int
        Number of patients to sample from variant data CSV file
    n_row_skip : int
        Number of rows corresponding to header information in variant CSV file
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
                if i < (n_row_skip):
                    continue
                elif i == (n_row_skip):
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
