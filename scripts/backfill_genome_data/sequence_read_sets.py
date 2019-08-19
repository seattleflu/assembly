"""
Find all fastq.gz files within a provided directory and group them by samples
into sequence read sets.

All sequence read sets are output to stdout as newline-delimited JSON
records.  You will likely want to redirect stdout to a file.
"""
import json
import os
import re
import argparse
from collections import defaultdict
from urllib.parse import urljoin
from pathlib import Path

def find_fasta_urls(fastq_directory, filename_pattern, url_prefix) -> dict:
    """
    Find all *.fastq.gz files within provided *fastq_directory*.

    The provided --filename-pattern regular expression is used to extract the
    sample ID from each FASTQ filename.  The regex should contain a capture
    group named "sample".  Each set of files with the same sample ID are
    grouped into a single sequence read set.

    Return as a dict with keys as NWGC sample IDs and values as an array of
    file paths.
    """
    sequence_read_sets = defaultdict(list)
    filename_pattern = re.compile(filename_pattern)

    for filepath in list(Path(fastq_directory).glob("*.fastq.gz")):
        # Skip the Undetermined FASTQ files
        if filepath.name.startswith("Undetermined"):
            continue
        filename = filepath.name
        # Check the filename matches provided filename pattern
        filename_match = filename_pattern.match(filename)
        assert filename_match, f"Filename {filename} doesn't match provided --filename-pattern"

        # Extract the sample from the filename_match
        try:
            sample = filename_match.group("sample")
        except IndexError:
            print(f"Filename {filename} matched provided --filename-pattern, but didn't extract a «sample» capture group")
            raise

        sequence_read_sets[sample].append(urljoin(url_prefix, str(filepath)))

    return sequence_read_sets


def dump_ndjson(sequence_read_set_dict):
    """
    Prints sequence read sets as newline-delimited JSON records.
    """
    for sample in sequence_read_set_dict:
        print(json.dumps({"sample": sample, "urls": sequence_read_set_dict[sample]}))


def dir_path(directory: str):
    """
    Confirm that the path for given *directory* is valid.
    """
    if os.path.isdir(directory):
        return directory
    else:
        raise NotADirectoryError(directory)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description = __doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--fastq-directory",
        metavar="<FASTQ directory>",
        type=dir_path,
        required=True,
        help="Absolute file path to directory that holds *.fastq.gz files")
    parser.add_argument("--filename-pattern",
        metavar="<regex>",
        default=r'^(?P<sample>\d+)_',
        help="Regex pattern to match sample in expected filename")
    parser.add_argument("--url-prefix",
        metavar="<url>",
        default = "file://rhino.fhcrc.org",
        help="Base for fully-qualifying sequence read set URLs")

    args = parser.parse_args()
    sequence_read_sets = find_fasta_urls(args.fastq_directory,
                                         args.filename_pattern,
                                         args.url_prefix)
    dump_ndjson(sequence_read_sets)
