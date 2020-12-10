"""
Calculate the percent identity for a consensus genome to the reference.

Assumes the consensus genome is generated through reference alignment and
therefore has the same length as the reference sequence.
"""
import argparse
import json
import numpy as np
from Bio import SeqIO


def get_segment_sequence(fasta_file: str, segment: str) -> np.ndarray:
    """
    Get the sequence for the *segment* in the *fasta_file*.
    Expects the *segment* to be in the FASTA sequence header divided by
    vertical pipes (|), e.g. 'HA' is in '>KX058884.1|VIC|HA'
    """
    for record in SeqIO.parse(fasta_file, "fasta"):
        sequence_ids = record.id.split("|")
        if segment in sequence_ids:
            return np.array(record.seq.lower())


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description = __doc__,
        formatter_class = argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--reference", type=str, required=True,
        help = "Reference FASTA file")
    parser.add_argument("--consensus", type=str, required=True,
        help = "Consensus genome FASTA file")
    parser.add_argument("--segment", type=str, required=True,
        help = "Segment of the genome to check percent identity, e.g. HA")
    parser.add_argument("--output", type=str, required=True,
        help = "The output of the checkpoint in Snakemake")

    args = parser.parse_args()

    reference_seq = get_segment_sequence(args.reference, args.segment)
    consensus_seq = get_segment_sequence(args.consensus, args.segment)

    # Because the consensus genome is generated from aligning mapped reads to
    # the reference genome, these two sequences should always have the same
    # length!
    assert len(reference_seq) == len(consensus_seq), \
        "Reference sequence and consensus genome sequence have different lengths"

    matches = np.sum(reference_seq == consensus_seq)
    percent_identity = (matches * 100) / len(reference_seq)

    with open(args.output, "w+") as output:
        json.dump({"percent_identity": percent_identity}, output)
