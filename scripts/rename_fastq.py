"""
Some fastq files come in with a different identifier that are not the 
NWGC sample ID. This renames the files so that they contain the corresponding
NWGC sample ID so it can be pushed through the assembly pipeline. 
"""
import os
import sys
import argparse
import pandas as pd

def rename_fastq_files(name_map: str, barcode: str, NWGC: str, file_dir: str):
    """
    Rename all the fastq files within *file_dir* according the the *name_map*.
    """
    mapping_df = pd.read_csv(name_map)
    mapper = mapping_df.set_index(barcode)[NWGC].to_dict()
    
    for filename in os.listdir(file_dir):
        if filename.endswith("fastq.gz"):
            split_file = filename.split("_", 2)
            file_barcode = "_".join(split_file[:2])
            try:
                file_id = str(mapper[file_barcode])
            except KeyError:
                print(f"Unknown barcode: {file_barcode}")
                continue
            else:
                new_filename = file_id + "_" + split_file[-1]
                os.rename(file_dir + "/" + filename, file_dir + "/" + new_filename)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Rename fastq files to contain NWGC sample ID",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("map", type=str, nargs="?",
        metavar="Sample-ID-Map",
        help="File path to csv file that matches random barcode to NWGC ID")
    parser.add_argument("barcode", type=str, nargs="?",
        metavar="Barcode-column",
        help="Column name of column containing random barcodes")
    parser.add_argument("NWGC", type=str, nargs="?",
        metavar="ID-column",
        help="Column name of column containing NWGC sample IDs")
    parser.add_argument("directory", type=str, nargs="?",
        metavar="FASTQ-directory",
        help="File path to directory of FASTQ files that need to be renamed")
    
    try:
        args = parser.parse_args()
        rename_fastq_files(args.map, args.barcode, 
                           args.NWGC, args.directory)
    except:
        parser.print_help()
        sys.exit(0)
    