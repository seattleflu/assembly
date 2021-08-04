"""
Creates BioSample submission TSV for previously submitted GenBank sequences
by pulling originating lab name for the same sequences from GISAID.

This should be a one-time script since we are submitting to BioSample
going forward.
"""
import argparse
import math
import numpy as np
import pandas as pd


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description = __doc__,
        formatter_class = argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--gisaid", type=str, required=True,
        help = "File path to the GISAID TSV file, expected to contain the originating lab name")
    parser.add_argument("--genbank", type=str, required=True, nargs='*',
        help = "File path to the GenBank metadata TSV file. Can list multiple files.")
    parser.add_argument("--output-dir", type=str, required=True,
        help = "Path to directory for the final output TSV files")

    args = parser.parse_args()

    gisaid = pd.read_csv(args.gisaid, sep='\t', dtype='string', usecols=['strain', 'originating_lab'])
    gisaid['Sequence_ID'] = gisaid['strain'].apply(lambda x: x.split('/')[1])

    genbank_columns = ['Sequence_ID', 'isolate', 'country', 'host', 'collection-date']
    genbank = pd.DataFrame()
    for file_path in args.genbank:
        genbank = genbank.append(pd.read_csv(file_path, sep='\t', dtype='string', usecols=genbank_columns))

    biosample = genbank.merge(gisaid, on=['Sequence_ID'], how='left')

    biosample_columns = [
        'sample_name',
        'bioproject_accession',
        'organism',
        'collected_by',
        'collection_date',
        'geo_loc_name',
        'host',
        'host_disease',
        'isolate',
        'isolation_source'
    ]

    biosample.rename(inplace=True, columns={
        'strain': 'sample_name',
        'country': 'geo_loc_name',
        'collection-date': 'collection_date',
        'originating_lab': 'collected_by'
    })

    # Hard-coded values
    biosample['bioproject_accession'] = 'PRJNA746979'
    biosample['organism'] = 'Severe acute respiratory syndrome coronavirus 2'
    biosample['host_disease'] = 'COVID-19'
    biosample['isolation_source'] = 'clinical'

    # Maximum number of records per submission to BioSample
    max_biosample = 1000
    number_of_files = math.ceil(len(biosample) / max_biosample)

    for index, chunk in enumerate(np.array_split(biosample, number_of_files), start=1):
        output_file = f'{args.output_dir}/biosample-{index}.tsv'
        chunk[biosample_columns].to_csv(output_file, sep='\t', index=False)
