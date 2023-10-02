"""
Merge identifiers with GISAID and GenBank accessions to keep track of all
sequence identifiers.
"""
import argparse
import pandas as pd
from typing import List
from create_submissions_sars_cov_2 import IDENTIFIER_COLUMNS
from create_submissions_rsv import IDENTIFIER_COLUMNS as RSV_IDENTIFIER_COLUMNS


def parse_identifiers(filename: str, pathogen: str) -> pd.DataFrame:
    """
    Parse the identifiers from the provided *filename*.
    File is expected to be a TSV with all IDENTIFIER_COLUMNS.

    Creates a new `seq_id` column for easy matching with GenBank Sequence IDs
    """
    identifiers = pd.read_csv(filename, sep='\t', dtype='string',
                                keep_default_na=False, usecols=RSV_IDENTIFIER_COLUMNS if pathogen.startswith('rsv-') else IDENTIFIER_COLUMNS )
    # Pull out the seq id from the strain name USA/<seq_id>/<year>
    identifiers['seq_id'] = identifiers['strain_name'].apply(
        lambda x: x.split('/')[1] if not x == 'N/A' else x)

    return identifiers


def parse_gisaid(filename: str) -> pd.DataFrame:
    """
    Parse the 'Accession ID' and 'Virus name' columns from the provided *filename*.
    Pulls out the `seq_id` from the `Virus name` column for easy merging with identifiers.
    """
    gisaid_columns = {
        'Accession ID': 'new_gisaid_accession',
        'Virus name': 'seq_id'
    }
    gisaid_accessions = pd.read_csv(filename, sep='\t', dtype='string', usecols=gisaid_columns.keys())
    gisaid_accessions.rename(columns=gisaid_columns, inplace=True)
    # Pull out the seq id from the GISAID virus name hCoV-19/USA/<seq_id>/<year>
    gisaid_accessions['seq_id'] = gisaid_accessions['seq_id'].apply(lambda x: x.split('/')[2])

    return gisaid_accessions


def parse_genbank(filenames: List[str]) -> pd.DataFrame:
    """
    Parse the '#Accession' and 'Sequence ID' columns from the provided *filenames*.
    Combines the data from multiple files into one DataFrame for easy merging.
    """
    genbank_columns = {
        '#Accession': 'new_genbank_accession',
        'Sequence ID': 'seq_id'
    }
    genbank_accessions = pd.DataFrame()
    for filename in filenames:
        genbank_accessions = genbank_accessions.append(
            pd.read_csv(filename, sep='\t', dtype='string', index_col=False, usecols=genbank_columns.keys()))

    return genbank_accessions.rename(columns=genbank_columns)


def merge_accessions(identifiers: pd.DataFrame,
                     accessions: pd.DataFrame,
                     accession_column: str) -> pd.DataFrame:
    """
    Merge the provided *identifiers* and *accessions* on the column `seq_id`.

    Expects *identifiers* to have the provided *accession_column*
    and the *accessions* to have a new_*accession_column*.

    Prints a warning if the new accession does not match the existing accession.
    """
    def fill_accession(row: pd.Series) -> str:
        """
        Given a *row* of identifiers, fills in the *accession_column*
        with the value from the new_*accession_column*.

        Prints a warning if row already has a value in the *accession_column*
        that does not match the new value.
        """
        new_accession_column = f'new_{accession_column}'

        old_accession = row[accession_column]
        new_accession = row[new_accession_column]

        if old_accession != 'N/A' and old_accession != new_accession:
            print(f"Warning: replacing the old accession {old_accession} with the new accession {new_accession}")

        return new_accession

    identifier_columns = identifiers.columns

    identifiers = identifiers.merge(accessions, on=['seq_id'], how='left')
    identifiers[accession_column] = identifiers.apply(fill_accession, axis=1)

    return identifiers[identifier_columns]



if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description = __doc__,
        formatter_class = argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--identifiers", type=str, required=True,
        help = "File path to identifiers TSV")
    parser.add_argument("--gisaid-accessions", type=str, required=False,
        help = "File path to GISAID accessions TSV")
    parser.add_argument("--genbank-accessions", type=str, required=False, nargs='*',
        help = "File path to GenBank accessions TSV. Can list multiple files.")
    parser.add_argument("--output", type=str, required=True,
        help = "File path for final output TSV.")
    parser.add_argument("--pathogen", type=str, required=False,
        choices=['sars-cov-2', 'rsv-a', 'rsv-b', 'flu-a', 'flu-b'],
        default='sars-cov-2')
    args = parser.parse_args()

    identifiers = parse_identifiers(args.identifiers, pathogen=args.pathogen)

    if args.gisaid_accessions:
        gisaid = parse_gisaid(args.gisaid_accessions)
        identifiers = merge_accessions(identifiers, gisaid, 'gisaid_accession')


    if args.genbank_accessions:
        genbank = parse_genbank(args.genbank_accessions)
        identifiers = merge_accessions(identifiers, genbank, 'genbank_accession')

    # For consistency, adding ' (Omicron)' to clade name to match previous format from NextClade
    if args.pathogen == 'sars-cov-2':
        identifiers.loc[identifiers.clade > '21J', 'clade'] = identifiers.clade + ' (Omicron)'
        identifiers.loc[identifiers.clade.isin(['21A', '21I', '21J']), 'clade'] = identifiers.clade + ' (Delta)'
        identifiers.loc[identifiers.clade == '20J', 'clade'] = '20J (Gamma, V3)'
        identifiers.loc[identifiers.clade == '20H', 'clade'] = '20H (Beta, V2)'
        identifiers.loc[identifiers.clade == '20I', 'clade'] = '20I (Alpha, V1)'

        identifiers[IDENTIFIER_COLUMNS].fillna('N/A').to_csv(args.output, sep='\t', index=False)
    elif args.pathogen.startswith('rsv-'):
        identifiers[RSV_IDENTIFIER_COLUMNS].fillna('N/A').to_csv(args.output, sep='\t', index=False)
