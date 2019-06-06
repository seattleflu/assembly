"""
Check to ensure R1 and R2 FASTQ files have the same reads in them.

This is to check that read-pair splitting is not occuring during
demultiplexing. 
"""
import sys
import gzip
import argparse

def check_r1_r2_reads(r1, r2):
    """
    Verify that reads in R1 and R2 are the same.
    """
    r1_read_ids = find_read_ids(r1)
    r2_read_ids = find_read_ids(r2)
    different_read_ids = r1_read_ids.symmetric_difference(r2_read_ids)
    percent_different = len(different_read_ids) / len(r1_read_ids) * 100
    if r1_read_ids == r2_read_ids:
        print("100% of R1 and R2 reads are the paired.")
        return 0
    print(f"R1 and R2 reads differ by {percent_different}%")
    print(different_read_ids)
    return percent_different
    

def find_read_ids(r1):
    """
    Pull out the read IDs FASTQ file
    """
    read_ids = set()
    with gzip.open(r1, "rt") as r1_file:
        for line in r1_file:
            if line.startswith("@"):
                id = line.split()[0]
                read_ids.add(id)
    return read_ids


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="""
        Check merged R1 and R2 FASTQ files have the same reads. 
        """,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("R1", type=str, nargs="?",
        metavar="R1",
        help="File path to merged R1 FASTQ file")
    parser.add_argument("R2", type=str, nargs="?",
        metavar="R2",
        help="File patht to merged F2 FASTQ file")
    
    
    try:
        args = parser.parse_args()
        check_r1_r2_reads(args.R1, args.R2)
    except:
        parser.print_help()
        sys.exit(0)