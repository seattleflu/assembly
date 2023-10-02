"""
Create files for submission of RSV sequences to GenBank.

Creates the following files:
- <batch_name>_<submission_group>_genbank_metadata.tsv: TSV to upload to GenBank's submission portal
- <batch_name>_<submission_group>_genbank.fasta: FASTA to upload to GenBank's submission portal
"""
import argparse
import sys
import pandas as pd
from datetime import datetime
from typing import List, Set, Optional
from pathlib import Path
from Bio import SeqIO
from utils import create_biosample_submission, create_genbank_submission

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description = __doc__,
        formatter_class = argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--batch-name", type=str, required=True,
        help = "The name for this batch of sequences, usually the date of sequence release, e.g. `20210701`")
    parser.add_argument("--biosample-metadata", type=str, required=True,
        help = "File path to the metadata CSV file")
    parser.add_argument("--biosample", type=str, required=True,
        help = "File path to the biosample accessions TSV file")
    parser.add_argument("--fasta", type=str, required=True,
        help = "File path to the FASTA file")
    parser.add_argument("--output-dir", type=str, required=True,
        help = "Path to the output directory for all output files")
    parser.add_argument("--pathogen", type=str, required=True,
        choices=['rsv-a', 'rsv-b', 'flu-a', 'flu-b'])

    args = parser.parse_args()

    biosample_metadata = pd.read_csv(args.biosample_metadata)
    biosample_accessions = pd.read_csv(args.biosample, sep='\t')

    biosample_metadata_with_accessions = pd.merge(biosample_metadata, biosample_accessions, on='sample_name', suffixes=(None,'_2'))

    create_genbank_submission(biosample_metadata_with_accessions, args.fasta, Path(args.output_dir), args.batch_name, args.pathogen)
