"""
Creates a directory of symlinks that point to FASTQs for recaptures and their
original FASTQs.
"""
import os
import argparse
import glob
import pandas as pd
from urllib.parse import urljoin

def create_symlink(sample: str, src_file: str, dst_directory: str) -> None:
    """
    Create a symlink in *dst_directory* that refers to the *src_file*.
    """
    src_path_split = src_file.split('/')
    src_directory = src_path_split[-2]
    src_filename = src_path_split[-1]
    new_file = '_'.join([str(sample), src_directory, src_filename])
    dst = urljoin(args.dst_directory, new_file)
    print(f"Creating a symlink «{dst}»")
    os.symlink(src_file, dst)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description = __doc__,
        formatter_class = argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument("--recaptures", type=str,
        help = "CSV file that lists samples in recapture run")

    parser.add_argument("--sample-column", type=str,
        help = "Name of sample column within recap CSV file")

    parser.add_argument("--recap-fastqs", type=str,
        help = "Directory of recaptured fastq.gz files, " +
               "must be full explicit path.")

    parser.add_argument("--original-fastqs", type=str, nargs="+",
        help = "Directories that contain original fastq.gz files, " +
               "must be full explicit path.")

    parser.add_argument("--dst-directory", type=str,
        help = "File path to directory that will contain symlinks")

    args = parser.parse_args()

    samples = pd.read_csv(
        args.recaptures,
        usecols=[args.sample_column],
        squeeze=True).tolist()

    for sample in samples:
        # Create symlinks for all recap FASTQs
        for f in glob.glob(f'{args.recap_fastqs}/{sample}*.fastq.gz'):
            create_symlink(sample, f, args.dst_directory)
        # Find original FASTQs for sample
        for directory in args.original_fastqs:
            # Skip directory if it is the recap directory
            if directory == args.recap_fastqs: continue
            print(f"Searching for sample «{sample}» in directory «{directory}»")
            # Create symlinks for original FASTQs
            for f in glob.glob(f'{directory}/{sample}*.fastq.gz'):
                create_symlink(sample, f, args.dst_directory)
