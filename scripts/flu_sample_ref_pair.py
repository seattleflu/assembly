"""
Generate sample/ref pairing based on presence/absence test results.

Currently need to download Excel file from Metabase that contains the 
NWGC sample ID and the target identifiers where present is True for the sample.
"""
import os
import sys
import json
import argparse
import pandas as pd

def generate_sample_ref_pairing(file_dir: str, presence_file: str,
                                sample: str, target:str, config_file:str):
    """
    For all the samples listed in *file_dir*, create a dict that maps each
    sample to a list of references.
    """
    reference_map = {
        "Flu_A_H3": {"h3n2_Texas_50_2012"},
        "Flu_A_H1": {"h1n1pdm_Michigan_45_2015"},
        "Flu_A_pan": {"h3n2_Texas_50_2012", "h1n1pdm_Michigan_45_2015"},
        "Flu_B_pan": {"vic_Brisbane_60_2008", "yam_Wisconsin_01_2010"}
    }

    present_refs = pd.read_excel(presence_file, usecols=[sample, target])
    present_refs[sample] = present_refs[sample].astype(str)
    present_refs[target] = present_refs[target].apply(lambda x: x.split(","))
    sample_mapper = present_refs.set_index(sample)[target].to_dict()
    
    config_map = {}
    ignore_samples = {}
    for filename in os.listdir(file_dir):
        if filename.endswith("fastq.gz"):
            sample_id = filename.split("_", 2)[0]
            try:
                target_list = sample_mapper[sample_id]
            except KeyError:
                print(f"{sample_id} not found. Adding to ignored samples")
                ignore_samples[sample_id] = {}
            else:
                config_targets = set()
                for target in target_list:
                    config_targets.update(reference_map[target])
                config_map[sample_id] = list(config_targets)

    with open(config_file, "r+") as f:
        data = json.load(f)
        data["sample_reference_pairs"] = config_map
        data["ignored_samples"] = ignore_samples
        f.seek(0)
        f.truncate()
        json.dump(data, f, indent=4)
        

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="""
        Generate sample/reference pairs for Flu only.
        Samples not found in «Sample/target map» are added to ignored samples.
        """,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("directory", type=str, nargs="?",
        metavar="FASTQ directory",
        help="File path to directory containing fastq files")
    parser.add_argument("sample-target-map", type=str, nargs="?",
        metavar="Sample/target map",
        help="File path to Excel file containing samples and present targets")
    parser.add_argument("sample", type=str, nargs="?",
        metavar="Sample column",
        help="Column name of column containg NWGC sample IDs")
    parser.add_argument("target", type=str, nargs="?",
        metavar="Target column",
        help="Column name of column containing target identifiers")
    parser.add_argument("config-file", type=str, nargs="?",
        metavar="Config file",
        help="File path to config file to be used with Snakemake")
    
    try:
        args = parser.parse_args()
        generate_sample_ref_pairing(args.file_dir, args.present, 
                                args.sample, args.target, args.config_file)
    except:
        parser.print_help()
        sys.exit(0)
