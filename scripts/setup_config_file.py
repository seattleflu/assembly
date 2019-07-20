"""
Edits the config file for Snakemake by filling in two fields:
"ignored_samples" and "sample_reference_pairs".

Currently need to download Excel file from Metabase that contains the
NWGC sample ID and the target identifiers where present is True for the sample.
"""
import glob
import sys
import json
import argparse
import gzip
import pandas as pd
from typing import Tuple

def edit_config_file(file_dir: str,
                     config_file:str,
                     presence_file: str,
                     sample: str,
                     target:str,
                     max_difference: int):
    """
    For all the samples listed in *file_dir*, create a dict that maps each
    sample to a list of references.
    """
    with open("references/target_reference_map.json", 'r') as f:
        reference_map = json.load(f)
    sample_target_map = create_sample_target_map(presence_file, sample, target)
    all_samples = sample_ids(file_dir)

    (sample_ref_pairs,
     ignored_samples) = generate_config_values(all_samples,
                                               file_dir,
                                               max_difference,
                                               sample_target_map,
                                               reference_map)

    with open(config_file, "r+") as f:
        data = json.load(f)
        data["sample_reference_pairs"] = sample_ref_pairs
        data["ignored_samples"] = dict.fromkeys(ignored_samples, {})
        f.seek(0)
        f.truncate()
        json.dump(data, f, indent=4)


def create_sample_target_map(presence_file: str,
                             sample: str, target: str) -> dict:
    """
    Create a map of samples and targets from the given *presence_file*,
    with *sample* as column name of samples and *target* as column
    name of targets.

    In the returned dict, the keys are the sample IDs and the values are a
    list of targets.
    """
    present_refs = pd.read_excel(presence_file, usecols=[sample, target])
    present_refs[sample] = present_refs[sample].astype(str)
    present_refs[target] = present_refs[target].apply(lambda x: x.split(","))
    return present_refs.set_index(sample)[target].to_dict()


def sample_ids(file_dir: str) -> set:
    """
    Create a set of all sample Ids from the provided *file_dir*
    """
    samples = set()
    for f in glob.glob(file_dir + "/*.fastq.gz"):
        sample_id = f.split("/")[-1].split("_", 2)[0]
        samples.add(sample_id)
    return samples


def generate_config_values(all_samples: set,
                           file_dir: str,
                           max_diff: int,
                           sample_target_map: dict,
                           ref_map: dict) -> Tuple[dict, set]:
    """
    Generates the sample_reference_pairs dict and ignored_samples set.

    Samples are added to ignored_samples if read difference between R1 and
    R2 is greater than *max_diff* or if the sample cannot be found in the
    *sample_target_map*.

    The sample_reference_pairs is made using the *sample_target_map* and
    the *ref_map* for the remaining samples.
    """
    sample_reference_pairs = {}
    ignored_samples = set()

    for sample in all_samples:
        if not sample_reads_pair(file_dir, sample, max_diff):
            print(f"{sample} has too many unpaired reads."
                  + "Adding to ignored samples")
            ignored_samples.add(sample)
        elif sample_target_map.get(sample) is None:
            print(f"{sample} not found. Adding to ignored samples")
            ignored_samples.add(sample)
        else:
            target_list = sample_target_map[sample]
            config_targets = set()
            for target in target_list:
                config_targets.update(ref_map[target])
            sample_reference_pairs[sample] = list(config_targets)

    return (sample_reference_pairs, ignored_samples)



def sample_reads_pair(file_dir: str, sample: str, max_diff: int) -> bool:
    """
    Check all *sample* files in *file_dir* to ensure reads are paired.

    Returns True if the difference between R1 and R2 is
    less than or equal to the *max_diff*
    """
    print(f"Sample: {sample}")
    for i in range(1, 5):
        print(f"Lane: {i}")
        r1 = glob.glob(file_dir + f"/{sample}*L00{i}*R1*.fastq.gz")[0]
        r2 = glob.glob(file_dir + f"/{sample}*L00{i}*R2*.fastq.gz")[0]
        difference = check_r1_r2_reads(r1, r2)
        if difference > max_diff:
            return False
    return True


def check_r1_r2_reads(r1: str , r2: str) -> int:
    """
    Checks the difference between the reads in R1 and R2.

    Returns the percent different as an int.
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


def find_read_ids(fastq_file: str) -> set:
    """
    Finds all read IDs in FASTQ file and returns them as a set.
    """
    read_ids = set()
    with gzip.open(fastq_file, "rt") as r1_file:
        for line in r1_file:
            if line.startswith("@"):
                id = line.split()[0]
                read_ids.add(id)
    return read_ids


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="""
        Edits the Snakemake config file by adding sample reference pairs and
        ignored samples.
        """,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("directory", type=str, nargs="?",
        metavar="FASTQ directory",
        help="File path to directory containing fastq files")
    parser.add_argument("config_file", type=str, nargs="?",
        metavar="Config file",
        help="File path to config file to be used with Snakemake")
    parser.add_argument("sample_target_map", type=str, nargs="?",
        metavar="Sample/target map",
        help="File path to Excel file containing samples and present targets")
    parser.add_argument("--sample", type=str, nargs="?",
        metavar="Sample column",
        default="sample",
        help="Column name of column containg NWGC sample IDs")
    parser.add_argument("--target", type=str, nargs="?",
        metavar="Target column",
        default="target",
        help="Column name of column containing target identifiers")
    parser.add_argument("--max_difference", type=int, nargs="?",
        metavar="Max percent difference",
        default=20,
        help="The maximum difference acceptable between R1 and R2 reads")


    try:
        args = parser.parse_args()
        edit_config_file(args.directory, args.config_file,
                         args.sample_target_map,
                         args.sample, args.target,
                         args.max_difference)
    except:
        parser.print_help()
        sys.exit(0)
