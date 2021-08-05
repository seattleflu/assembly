"""
Adds BioSample accessions to GenBank metadata TSVs for submissions.

Outputs the new TSVs to the same directory as the provided GenBank TSVs
with the same filename + `_with_biosample` suffix.
"""
import argparse
import pandas as pd
from pathlib import Path


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description = __doc__,
        formatter_class = argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--biosample", type=str, required=True,
        help = "File path to the BioSample accessions TSV")
    parser.add_argument("--genbank", type=str, required=True, nargs='*',
        help = "File path to the GenBank metadata TSV file. Can list multiple files")

    args = parser.parse_args()

    biosample_columns = ['accession', 'isolate', 'bioproject_accession']
    biosample = pd.read_csv(args.biosample, sep='\t', dtype='string', usecols=biosample_columns)
    biosample.drop_duplicates(inplace=True, keep='first')
    # Rename BioSample columns to fit the GenBank submissions template
    biosample.rename(inplace=True, columns={
        'accession': 'BioSample',
        'bioproject_accession': 'BioProject'
    })

    for filepath in args.genbank:
        genbank = pd.read_csv(filepath, sep='\t', dtype='string')
        genbank = genbank.merge(biosample, on=['isolate'], how='left')

        input_file = Path(filepath)
        output_file = input_file.with_name(f'{input_file.stem}_with_biosample.tsv')
        genbank.to_csv(output_file, sep='\t', index=False)

