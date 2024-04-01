"""
Creates genbank sequence accession tsv for flu tracking,
which maps strain name and sequence id to GenBank accession id.
"""

import pandas as pd
import argparse


def parse_genbank_metadata(filename: str, pathogen: str) -> pd.DataFrame:
    genbank_columns = {
        'flu-a': {
            'Sequence_ID': 'sequence_id',
            'isolate': 'strain_name',
            'serotype':  'subtype'
        },
        'flu-b': {
            'Sequence_ID': 'sequence_id',
            'isolate': 'strain_name'
        }
    }

    genbank_metadata = pd.read_csv(filename, sep='\t', dtype='string',
                                keep_default_na=False, usecols=genbank_columns[pathogen].keys())

    genbank_metadata.rename(columns=genbank_columns[pathogen], inplace=True)

    if pathogen == 'flu-a':
        genbank_metadata['subtype'] = genbank_metadata['subtype'].apply(lambda x: x.lower())
    elif pathogen == 'flu-b':
        # add subtype column and set to N/A. This is because we have a single tracker for flu-a and flu-b.
        genbank_metadata['subtype'] = None

    return genbank_metadata


def parse_genbank_accessions(filename: str) -> pd.DataFrame:
    accessions_columns = {
        '#Accession': 'genbank_accession',
        'Sequence ID': 'sequence_id'
    }

    genbank_accessions = pd.read_csv(filename, sep='\t', dtype='string', 
                                     index_col=False, usecols=accessions_columns.keys())

    return genbank_accessions.rename(columns=accessions_columns)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description = __doc__,
        formatter_class = argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--genbank-metadata", type=str, required=True,
        help = "File path to genbank metadata TSV")
    parser.add_argument("--genbank-accessions", type=str, required=True,
        help = "File path to GenBank accessions TSV")
    parser.add_argument("--batch-date", type=str, required=True,
        help = "Batch date, formatted YYYYMMDD")
    parser.add_argument("--output", type=str, required=True,
        help = "File path for final output TSV")
    parser.add_argument("--pathogen", type=str, required=True,
        choices=['flu-a', 'flu-b'])
    args = parser.parse_args()


    genbank_metadata = parse_genbank_metadata(args.genbank_metadata, args.pathogen)

    accessions = parse_genbank_accessions(args.genbank_accessions)

    # merge metadata to accessions
    full_df = genbank_metadata.merge(accessions, on='sequence_id')

    # add batch name and pathogen
    full_df['pathogen'] = args.pathogen
    full_df['batch'] = f'{args.batch_date}_fastq'

    # reorder columns
    cols = [
        'sequence_id',
        'batch',
        'strain_name',
        'genbank_accession',
        'pathogen',
        'subtype'
    ]

    full_df = full_df[cols]

    # write out to tsv
    full_df.fillna('N/A').to_csv(args.output, sep='\t', index=False)

