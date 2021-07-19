"""
Pull out a subset of sequences for re-submission provided the metadata CSV/TSV,
the FASTA, and a list of the subset of sequence names.
"""
import os
import argparse
import pandas as pd
from typing import List
from Bio import SeqIO
from create_submissions import text_to_list


def create_output_filename_and_file_path(input_file_path: str, output_dir:str) -> List[str]:
    """
    Given the *input_file_path* and *output_dir* create the output filename,
    which is the original filename with a `-subset` added, and the full
    output file path.
    """
    full_filename = os.path.basename(input_file_path)
    filename, file_extension = os.path.splitext(full_filename)

    output_filename = filename + '-subset' + file_extension
    output_file_path = os.path.join(output_dir, output_filename)

    return [output_filename, output_file_path]


def extract_fasta_sequences(subset: set, input_fasta: str, output_dir: str) -> str:
    """
    Extract sequences from *input_fasta* that are listed in provided *subset*
    and create a new FASTA file in the *output_dir*.

    Returns the output FASTA filename to be used in downstream
    metadata processing.
    """
    subset_sequences = []
    for record in SeqIO.parse(input_fasta, 'fasta'):
        if record.id in subset:
            subset_sequences.append(record)

    output_fname,output_path = create_output_filename_and_file_path(input_fasta, output_dir)

    SeqIO.write(subset_sequences, output_path, 'fasta-2line')

    return output_fname


def extract_metadata(subset: set, input_metadata: str, metadata_format: str,
                     output_fasta: str, output_dir: str):
    """
    Extract the metadata for sequences in the provided *subset*.
    If the *metadata_format* is GISAID, then also replace the FASTA filename
    with the provided *output_fasta*.

    The extracted metadata is output to a new file in the *output_dir*.
    """
    sequence_column = output_fasta_column = None
    if metadata_format == 'GISAID':
        sequence_column = 'covv_virus_name'
        output_fasta_column = 'fn'
    else:
        sequence_column = 'Sequence_ID'

    separator = '\t' if input_metadata.endswith('tsv') else ','
    metadata = pd.read_csv(input_metadata, sep=separator)

    if output_fasta_column and output_fasta:
        metadata.loc[metadata[sequence_column].isin(subset), output_fasta_column] = output_fasta

    subset_metadata = metadata.loc[metadata[sequence_column].isin(subset)]

    output_metadata_fname,output_metadata_path = create_output_filename_and_file_path(input_metadata, output_dir)
    subset_metadata.to_csv(output_metadata_path, sep=separator, index=False)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description = __doc__,
        formatter_class = argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--subset", type=str,
        help = "File path to the text file that lists a subset of sequence names in a single column")
    parser.add_argument("--metadata", type=str,
        help = "File path to the metadata CSV/TSV file")
    parser.add_argument("--metadata-format", choices=['GISAID', 'GenBank'],
        help = "The format of the metdata file provided. Choices are GISAID or GenBank")
    parser.add_argument("--fasta", type=str,
        help = "File path to the FASTA file.")
    parser.add_argument("--output-dir", type=str,
        help = "File path to the output directory")

    args = parser.parse_args()

    subset = set(text_to_list(args.subset))

    output_fasta = None
    if args.fasta:
        output_fasta = extract_fasta_sequences(subset, args.fasta, args.output_dir)

    extract_metadata(subset, args.metadata, args.metadata_format, output_fasta, args.output_dir)
