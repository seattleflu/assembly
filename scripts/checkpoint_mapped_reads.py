"""
Calculate minimum number of reads necessary to meet minimum depth for a
reference and determine the number of mapped reads from bowtie2 alignment
summary.

Print both to the <output>.
"""
import argparse
from Bio import SeqIO

def determine_ref_genome_length(reference):
    """
    Calculate the genome length of a given *reference*
    """
    genome_length = 0
    for segment in SeqIO.parse(reference, 'fasta'):
        seq = segment.seq
        seq_length = len(seq)
        genome_length += seq_length
    return genome_length


def determine_number_of_mapped_reads(bowtie2):
    """
    Use the given *bowtie2* alignment summary to determine the total number of
    reads that aligned with the reference
    """
    mapped_reads = 0
    with open(bowtie2, "r") as summary:
        for line in summary.readlines():
            if "1 time" not in line:
                continue
            aligned = int(line.split()[0])
            mapped_reads += aligned
    return mapped_reads


def print_results(minimum_reads, mapped_reads, output_file):
    """
    Print *minimum_reads* and *mapped_reads* in provided *output_file*
    """
    with open(output_file, "w") as output:
        output.write(f"Minimum reads required: {str(minimum_reads)}\n" + \
            f"Mapped reads: {str(mapped_reads)}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description = __doc__,
        formatter_class = argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--bowtie2", type=str, nargs="?", required=True,
        help = "Summary log from bowtie2")
    parser.add_argument("--min-cov", type=int, required=True,
        help = "Minimum coverage set for varscan")
    parser.add_argument("--raw-read-length", type=int, required=True,
        help = "Length of raw reads in FASTQ files")
    parser.add_argument("--reference", type=str, nargs="?", required=True,
        help = "Reference FASTA file")
    parser.add_argument("--output", type=str, nargs="?", required=True,
        help = "The output of the checkpoint in Snakemake")

    args = parser.parse_args()

    genome_length = determine_ref_genome_length(args.reference)
    minimum_reads = (genome_length * args.min_cov) / args.raw_read_length
    mapped_reads = determine_number_of_mapped_reads(args.bowtie2)
    print_results(minimum_reads, mapped_reads, args.output)
