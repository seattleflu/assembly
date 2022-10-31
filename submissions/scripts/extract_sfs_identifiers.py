"""
Pull out NWGC ID and SFS sample barcodes from a provided NWGC metadata Excel
file and output them in a two column CSV file.

Will do a simple check to verify SFS sample barcodes are 8 characters made up
of 0-9 and a-f. Prints barcodes that do not follow this format to stdout.
"""
import argparse
import sys
import pandas as pd
from typing import Optional


def read_all_identifiers(metadata_file: str) -> pd.DataFrame:
    """
    Read all sheets from *metadata_file* and return identifiers from the
    preferred sheet. Refer to the SHEET_PREFERENCE for the order of sheet
    preference.
    """
    SHEET_PREFERENCE = ['Samplify FC data', 'Metadata']
    NWGC_COLUMN_MAP = {
        'Samplify FC data' : {
            'Sample ID': 'nwgc_id',
            'Project': 'project',
            'Investigator\'s sample ID': 'sfs_sample_barcode',
            'Origin': 'sample_origin',
        },
        'Metadata' : {
            'LIMS': 'nwgc_id',
            'Project': 'project',
            'lab_accession_id': 'sfs_sample_barcode',
            'county': 'sample_origin'
        },
    }

    all_sheets = pd.read_excel(metadata_file, dtype='string',
                               engine='openpyxl', sheet_name=None)

    for sheet in SHEET_PREFERENCE:
        if sheet in all_sheets.keys():
            identifiers = all_sheets[sheet]

            # The "Samplify FC data" sheet has an extra row on top
            # Rename the columns based on values in the first row and delete the row
            if sheet == 'Samplify FC data':
                identifiers.columns = identifiers.iloc[0]
                identifiers.drop(identifiers.index[0], inplace=True)

            identifiers.rename(columns = NWGC_COLUMN_MAP[sheet], inplace=True)

            return identifiers[NWGC_COLUMN_MAP[sheet].values()]

    sys.exit(f"No known sheet names found in metadata Excel file. It contains the following sheet names: {all_sheets.keys()}")


def find_sfs_identifiers(identifiers: pd.DataFrame, sample_origin_filter=None) -> Optional[pd.DataFrame]:
    """
    Find SFS identifiers within provided *identifiers* based on regex matches
    of the 'project' column.

    Outputs incorrectly formatted barcodes to stdout.
    """
    SFS_PROJECT_REGEX = r'(SeattleChildrensDirect_|starita_bbi_)'
    SFS_BARCODE_REGEX = r'^[a-f0-9]{8}$'
    SFS_OUTPUT_COLUMNS = ['nwgc_id', 'sfs_sample_barcode', 'sample_origin']

    sfs_samples = identifiers.loc[identifiers['project'].str.match(SFS_PROJECT_REGEX)]

    if len(sfs_samples.index) == 0:
        print("No Seattle Flu Study samples found!")
        return None

    bad_barcodes = sfs_samples.loc[~sfs_samples['sfs_sample_barcode'].str.match(SFS_BARCODE_REGEX)]

    if len(bad_barcodes.index) > 0:
        print(f"Found the following erroneous barcodes: {list(bad_barcodes['sfs_sample_barcode'].values)}")

    if sample_origin_filter:
        return sfs_samples[sfs_samples['sample_origin'].str.contains(sample_origin_filter)][SFS_OUTPUT_COLUMNS]
    else:
        return sfs_samples[SFS_OUTPUT_COLUMNS]



if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description = __doc__,
        formatter_class = argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--metadata", type=str, required=True,
        help = "File path to the metadata Excel file provided by NWGC")
    parser.add_argument("--output", type=str, required=True,
        help = "File path to the output CSV file")

    args = parser.parse_args()

    identifiers = read_all_identifiers(args.metadata)
    sfs_identifiers = find_sfs_identifiers(identifiers)

    if sfs_identifiers is not None:
        sfs_identifiers.to_csv(args.output, index=False)
