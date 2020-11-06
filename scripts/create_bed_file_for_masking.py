"""
Create a bed file from the pileup file that contains all the sites where
coverage is below a set minimum.
"""
import argparse
import re


def create_bed_file(pileup, min_coverage, min_freq, bed_file):
    """
    """
    number_of_low_coverage = 0
    with open(pileup, 'r') as pile:
        with open(bed_file, 'w') as bed:
            for line in pile:
                line = line.split()
                sequence = line[0]
                end_position = line[1]
                start_position = str(int(end_position) - 1)
                coverage_depth = line[3]
                bed_line = sequence + "\t" + start_position + "\t" + end_position + "\t"
                if int(coverage_depth) < min_coverage:
                    bed.write(bed_line + "coverage depth: " + coverage_depth + "\n")
                    number_of_low_coverage += 1
            bed.close()
        pile.close()
    print(f"Number of low coverage (<{min_coverage}) positions: {number_of_low_coverage}")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description = __doc__,
        formatter_class = argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--pileup", type=str, nargs="?",
        help = "Pileup file generated from samtools mpileup")
    parser.add_argument("--min-cov", type=int, nargs="?",
        help = "The minimum coverage depth needed to not mask a base")
    parser.add_argument("--min-freq", type=float, nargs="?",
        help = "The minimum frequency of a SNP to be included.")
    parser.add_argument("--bed-file", type=str, nargs="?",
        help = "File path to the output bedfile")

    args = parser.parse_args()
    create_bed_file(args.pileup, args.min_cov, args.min_freq, args.bed_file)
