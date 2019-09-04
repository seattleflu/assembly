"""
Converts a fasta file input to machine-readable, Json format. The Json document
is printed to stdout, so the user will likely want to redirect this to a new
file.
"""
import re
import json
import argparse
from Bio import SeqIO


def fasta_to_json(fasta_file):
    consensus_genomes = list(SeqIO.parse(fasta_file, "fasta"))

    organism = None
    sample_identifier = None
    converted_genomes = []

    for segment in consensus_genomes:
        header = segment.id.split("|")

        sample_identifier = check_sample_identifier(header, sample_identifier)
        sequence_identifier = header[1]
        organism = check_organism(header, organism)
        sequence_segment = header[3]

        sequence = str(segment.seq)

        genome_data = {
            "sequence_identifier": sequence_identifier,
            "sequence_segment": sequence_segment,
            "genomic_sequence": sequence
        }
        converted_genomes.append(genome_data)

    print(json.dumps(
        {
            'sample_identifier': sample_identifier,
            'reference_organism': organism,
            'masked_consensus': converted_genomes
        }, indent=4))


def check_sample_identifier(header: str, sample_id: str) -> str:
    """"""
    if sample_id and header[0] != sample_id:
        raise AssertionError("There are multiple sample identifiers in the given fasta file.")

    return header[0]


def check_organism(header: str, organism_id: str) -> str:
    """"""
    if organism_id and header[2] != organism_id:
        raise AssertionError("There are multiple target organisms in the given fasta file.")

    return header[2]


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description = __doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("fasta_file", type=str, nargs="?",
        metavar = "<batch.fasta>",
        help="File path to FASTA file containing a batch of consensus genomes")

    args = parser.parse_args()
    fasta_to_json(args.fasta_file)
