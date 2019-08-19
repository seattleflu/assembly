"""
Converts a batch FASTA file into newline-delimited JSON file, where each line
is a record of one sequence.

Sequences are output to stdout, you will likely want to redirect
stdout to a file.
"""
import json
import argparse
from Bio import SeqIO


def fasta_to_ndjson(fasta_file):
    """
    Prints each sequence as a newline-delimited JSON record.
    """
    consensus_genomes = list(SeqIO.parse(fasta_file, "fasta"))

    for segment in consensus_genomes:
        header = segment.id.split("|")
        sample_identifier = header[0]
        sequence_identifier = header[1]
        organism = header[2]
        sequence_segment = header[3]
        sequence = str(segment.seq)
        genome_data = {
            "sample_identifier": sample_identifier,
            "sequence_identifier": sequence_identifier,
            "organism": organism,
            "sequence_segment": sequence_segment,
            "genomic_sequence": sequence
        }
        print(json.dumps(genome_data))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description = __doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--fasta-file",
        type=argparse.FileType('r'),
        metavar = "<batch.fasta>",
        help="File path to FASTA file containing a batch of consensus genomes")

    args = parser.parse_args()
    fasta_to_ndjson(args.fasta_file)
